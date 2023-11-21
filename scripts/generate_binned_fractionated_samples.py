#!/usr/bin/env python3
import sys, csv
import pandas as pd

from waterYearSamplesInfrastructure import FractionatedSIPSample, DNAextractionSample

NOVOGENE_DATA = '../Sequencing/QC_fract/QC_data.csv'
UCB_VOLUME_DATA = 'binned_expanded_wells_for_volume_input.csv'
UCB_DNA_DATA = 'binned_fractions_concentrations.csv'
DENSITY_DATA = 'binned_fractions_average_density.csv'
SAMPLE_METADATA = '../all_samples.csv'

def build_dna_sample_objects(metadata_file, seq_ids):
    dna_samples = dict()
    with open(metadata_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row['Extraction ID'] not in seq_ids:
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
            dna_samples[row['Extraction ID']] = current_sample
    return dna_samples

def main():
    # Combine novogene, UCB volume, and UCB DNA concentration data
    novogene_data = pd.read_csv(NOVOGENE_DATA)
    ucb_volume_data = pd.read_csv(UCB_VOLUME_DATA)
    ucb_dna_data = pd.read_csv(UCB_DNA_DATA)
    density_data = pd.read_csv(DENSITY_DATA)
    combined_data = pd.merge(novogene_data, ucb_dna_data, on='sample_id')
    combined_data = pd.merge(combined_data, ucb_volume_data, on='sample_id')
    combined_data = pd.merge(combined_data, density_data, on='sample_id')
    #Group by sample_id
    grouped_data = combined_data.groupby('sample_id', as_index=False).agg({'volume': 'sum', 
                                                                           'well': lambda x: ','.join(x), 
                                                                           'plate': 'min', 
                                                                           'batch':'min', 
                                                                           'Concentration': 'mean', 
                                                                           'volume':'sum', 
                                                                           'Concentration(ng/ul)':'mean', 
                                                                           'Volume(ul)':'mean',
                                                                           'average_density': 'mean' })
    #Rename concentration and volume columns
    grouped_data.rename(columns={'Concentration':'conc UCB', 
                                 'volume':'vol UCB', 
                                 'Concentration(ng/ul)':'conc NG', 
                                 'Volume(ul)':'vol NG'}, inplace=True)
    grouped_data['conc UCB'] = grouped_data['conc UCB'].round(2)
    grouped_data['conc NG'] = grouped_data['conc NG'].round(2)
    grouped_data['average_density'] = grouped_data['average_density'].round(4)
    dna_samples = build_dna_sample_objects(SAMPLE_METADATA, [x[:-2] for x in grouped_data['sample_id'].to_list()])

    with open("binned_fractionated_samples_QR_create.csv", 'w') as csvfile:
        csvfile.write("sample_id,core_id,date_collected,viromics,host,incubation_length\n")
        for i, row in grouped_data.iterrows():
            sample = dna_samples[row['sample_id'][:-2]]
            csvfile.write(f'{str(sample)}{row["sample_id"][-2:]},{sample.core_id},2023-9-1,0,LLNL-UCB,{sample.experiment[-1]}\n')

    with open("binned_fractionated_samples_QR_update.csv", 'w') as csvfile:
        csvfile.write("sample_id,Extraction ID,density,batch,plate,wells,volume-UCB,conc-UCB,volume-NG,conc-NG\n")
        for i, row in grouped_data.iterrows():
            sample = dna_samples[row['sample_id'][:-2]]
            csvfile.write(f'{str(sample)}{row["sample_id"][-2:]},{row["sample_id"]},{row["average_density"]},{row["batch"]},{row["plate"]},{row["well"].replace(",","-")},{row["vol UCB"]},{row["conc UCB"]},{row["vol NG"]},{row["conc NG"]}\n')

if __name__ == '__main__':
    main()