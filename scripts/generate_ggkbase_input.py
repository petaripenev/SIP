#!/usr/bin/env python3
import sys, csv, argparse

from waterYearSamplesInfrastructure import FractionatedSIPSample, DNAextractionSample

SAMPLE_METADATA = '../all_samples.csv'
#SEQUENCING_FILE_NAMES = '../Sequencing/Final_report_unfract/X202SC23061317-Z01-F001.csv'
SEQUENCING_FILE_NAMES = '../Sequencing/Final_report_fract/md5_check.csv'
ECOSYSTEM = {'H': "HREC", 'A': "Angelo"}
FRACTIONS = True

def build_dna_sample_objects(metadata_file, seq_file_names, use_fractionated_samples=False):
    seq_ids = [x[1:5] for x in seq_file_names]
    if use_fractionated_samples:
        seq_ids = [x.split('_')[0] for x in seq_file_names]
    dna_samples = list()
    with open(metadata_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            sample_identifier = row['id']
            if use_fractionated_samples:
                sample_identifier = row['Extraction ID']
            if sample_identifier not in seq_ids:
                continue
            date_collected = row['date_collected']
            if row['sample_id'][2] == '8':
                date_collected = '2022-09-15'
            current_sample = DNAextractionSample(
                        id = row['id'],
                        sampling_site = row['location'][0],
                        replicate = row['location'][1],
                        sampling_week = row['sample_id'][2],
                        depth = tuple(row['depth'].split('-')),
                        depth_str = row['sample_id'][4:6],
                        date = date_collected,
                        metadata_file_path='../all_samples.csv',
                        experiment = row['sample_id'][7:])
            current_sample.extraction_id = row['Extraction ID']
            current_sample.seq_file_names = [x for x in seq_file_names if current_sample.id == x[1:5]]
            if use_fractionated_samples:
                current_sample.seq_file_names = [x for x in seq_file_names if current_sample.extraction_id == x.split('_')[0]]
            dna_samples.append(current_sample)
    return dna_samples

def main(commandline_arguments):
    seq_file_names = list()
    with open(SEQUENCING_FILE_NAMES) as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if '.fq.gz' not in row[0]:
                continue
            seq_file_names.append(row[0].split('/')[-1])
    
    dna_samples = build_dna_sample_objects(SAMPLE_METADATA, seq_file_names, use_fractionated_samples=FRACTIONS)
    print(f'raw_forward_reads_name\traw_reverse_reads_name\tggkbase_project_name\tnotes')
    
    if FRACTIONS:
        fraction_dict = dict()
        for sample in dna_samples:
            for file in sample.seq_file_names:
                extraction_id = f'{file.split("_")[0]}_{file.split("_")[1]}'
                if extraction_id not in fraction_dict.keys():
                    fraction_dict[extraction_id] = list()
                fraction_num = file.split('_')[1]
                note = f'Fraction {fraction_num} sequencing from {ECOSYSTEM[sample.sampling_site]} plot {sample.replicate}, sampled on {sample.date} at {sample.depth[0]}-{sample.depth[1]}cm depth.'
                project_name = f'WaterYear_{sample}_{fraction_num}'
                fraction_dict[extraction_id].append([file, project_name, note])
        for extraction_id in fraction_dict.keys():
            for fw, rev in zip(fraction_dict[extraction_id][0::2], fraction_dict[extraction_id][1::2]):
                if fw[0][-7] != '1':
                    raise ValueError(f'Forward file {fw[0]} does not end in 1!')
                if rev[0][-7] != '2':
                    raise ValueError(f'Reverse file {rev[0]} does not end in 2!')
                print(f'{fw[0]}\t{rev[0]}\t{fw[1]}\t{fw[2]}')
    else:
        for sample in dna_samples:
            for fw, rev in zip(sample.seq_file_names[0::2], sample.seq_file_names[1::2]):
                print(f'{fw}\t{rev}\tWaterYear_{sample}\tUnfractionated sequencing from {ECOSYSTEM[sample.sampling_site]} plot {sample.replicate}, sampled on {sample.date} at {sample.depth[0]}-{sample.depth[1]}cm depth.')


if __name__ == '__main__':
    main(sys.argv[1:])