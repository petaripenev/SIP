#!/usr/bin/env python3
from glob import glob
import re, csv, os
import pandas as pd 
from dataclasses import dataclass, field

TUBES_TO_PLATES = {'A' : 1,
                   'B' : 1,
                   'C' : 1,
                   'D' : 1,
                   'E' : 2,
                   'F' : 2,
                   'G' : 2,
                   'H' : 2,
                   'I' : 3,
                   'J' : 3,
                   'K' : 3,
                   'L' : 3}

@dataclass
class SIPSample:
    sampling_site: str
    replicate: str
    sampling_week: str
    depth: "tuple[int]"
    isotope: str
    incubation_length: int
    sample_weight: float
    grams_soil_extracted: float
    concentration_DNA_extracted: float
    total_ul_DNA_extracted: float
    ul_DNA_SIP_loaded: float
    h20_added: float
    initial_GWC: float
    plate: str
    tube: str
    dna_yield: float
    concentrations: "list[float]" = field(default_factory=list, repr=False)
    densities: "list[float]" = field(default_factory=list, repr=False)
    wells: "list[str]" = field(default_factory=list, repr=False)

    def __post_init__(self):
        self.sample_id = self.__repr__()
        self.oneword_id = self.sample_id.replace('-','').replace('_','')

    def __repr__(self):
        return f"{self.sampling_site}{self.replicate}{self.sampling_week}_{self.depth[0]}-{self.depth[1]}_{self.isotope}O-{self.incubation_length}"

def getGWC(sampleID, filePath):
    '''Get the GWC of a sample from the metadata file'''
    with open(filePath) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row['sample_id'] == sampleID:
                return float(row['GWC %'])

def parseSIPoutput(df):
    '''Parse the SIP output from Mike's script into a dictionary of SIPSamples.
    Skips wells excluded from the analysis (first one and last two)'''
    sipData = dict()
    for col in df.columns:
        if col == 'Sample ID':
            continue
        if re.match(r'^Unnamed:', str(col)):
            continue
        if re.search(r'-density$', str(col)) or re.search(r'-conc$', str(col)):
            continue
        if len(str(col)) < 4:
            continue

        sipData[col] = {'concentrations':list(df[f'{col}-conc'][3:24]),
                    'densities'     :list(df[f'{col}-density'][3:24]),
                    'wells'         :list(df[col][3:24]),
                    'plate'         :TUBES_TO_PLATES[df[col][0]],
                    'tube'          :df[col][0],
                    'dna_yield'     :df[f'{col}-conc'][24]}
                        
    return sipData

def createSIPSampleClasses(SIPSample, sipData, filePath):
    sipSamples = list()
    with open(filePath) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if int(row['id']) not in sipData.keys():
                continue
            gwc_before_water_addition = float(row['Lab dried GWC (%)'])
            if gwc_before_water_addition == '':
                gwc_before_water_addition = getGWC(f"{row['sample_id'][0:6]}_GSC", filePath)
            sample_weight = float(row['Weight dried (gr)'])
            if sample_weight == '':
                sample_weight = float(row['Weight (gr)'])
            sipData[int(row['id'])]['sample_id'] = row['sample_id']
            newSample = SIPSample(row['sample_id'][0],
                              row['sample_id'][1],
                              row['sample_id'][2],
                              (row['depth'].split('-')[0],row['depth'].split('-')[1]),
                              row['sample_id'][7:9],
                              row['sample_id'].split('-')[1],
                              sample_weight,
                              float(row['Grams soil DNA extr1']),
                              float(row['DNA concentration extr1 (ng/ul)']),
                              float(row['Volume extr1 (ul)']),
                              float(row['DNA loaded for SIP (ul)']),
                              float(row['H2O amount (ml)']),
                              gwc_before_water_addition,
                              sipData[int(row['id'])]['plate'],
                              sipData[int(row['id'])]['tube'],
                              sipData[int(row['id'])]['dna_yield'],
                              sipData[int(row['id'])]['concentrations'],
                              sipData[int(row['id'])]['densities'],
                              sipData[int(row['id'])]['wells'])
            sipSamples.append(newSample)
    return sipSamples

def main():
    #Read in the SIP output
    df = pd.read_excel("./metadata/Angelo_W4_SIP.xlsx", sheet_name='Summary', engine='openpyxl')

    sipData = parseSIPoutput(df)
    sipSamples = createSIPSampleClasses(SIPSample, sipData, './metadata/all_samples_afterW10.csv')

    # Iterate over files in folder ./fastq
    fastq = glob('./fastq/*.fastq')
    print("file_name,new_name,sample_id")
    for file in fastq:
        # Get the part of filename before L001
        sampleID = re.search(r'(./fastq/)(.*).fastq', file).group(2).split('_')[0]
        readPair = re.search(r'(./fastq/)(.*).fastq', file).group(2).split('_')[3]
        if sampleID == 'negative' or sampleID[:-1] == 'negative-control':
            continue
        if sampleID == 'positive-control' or sampleID == 'Undetermined':
            continue
        redo=''
        split_sampleID = sampleID.split('-')
        if len(split_sampleID) == 3:
            redo = f"-{split_sampleID[2]}"
        plate = int(sampleID.split('-')[0][-1])
        well = sampleID.split('-')[1]
        for sample in sipSamples:
            if sample.plate == plate and well in sample.wells:
                print(f"{file},{sample.oneword_id}{plate}{well}{redo}_{readPair}_.fastq,{sample.sample_id}_{plate}_{well}_{readPair}{redo}")
                os.rename(file, f"./fastq/{sample.oneword_id}{plate}{well}{redo}_{readPair}_.fastq")

    # #Print out the samples
    # for sample in sipSamples:
    #     print(sample)
    

if __name__ == "__main__":
    main()