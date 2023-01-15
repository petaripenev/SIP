#!/usr/bin/env python3
import csv
from pandas import read_excel
from rename_MiSeq_outputs import parseSIPoutput, createSIPSampleClasses, SIPSample

ABC_TO_NUMBERS = {'A':1, 'B':2, 'C':3}
FINAL_SIP_VOLUME = 36.0

def laodTaxonomyTable(taxonomyTablePath):
    '''Load the taxonomy table from the sequencing pipeline output in a list of zOTUs and a list of their taxonomies'''

    zOTUs, taxonomies = list(), list()
    with open(taxonomyTablePath, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for i, row in enumerate(reader):
            if i == 0:
                continue
            zOTUs.append(row[0].replace('ASV','X'))
            taxonomies.append(row[1].replace(' ',''))#.replace('d__','').replace('p__','').replace('c__','').replace('o__','').replace('f__','').replace('g__','').replace('s__',''))
    return zOTUs, taxonomies

def loadZOTUtable(zOTUtablePath):
    '''Load the zOTU table from the sequencing pipeline output in a dict of dicts with the zOTU as key and the sample as key and the abundance as value'''

    zOTUtable, sampleIDsWithWells = dict(),list()
    with open(zOTUtablePath, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for i, row in enumerate(reader):
            renamedASV = row[0].replace('ASV','X')
            if i == 0:
                sampleIDsWithWells = row[1:]
                continue
            zOTUtable[renamedASV] = dict()
            for j, abundance in enumerate(row[1:]):
                zOTUtable[renamedASV][sampleIDsWithWells[j]] = abundance

    return zOTUtable

def getASVsPerSample(zOTUs, zOTUtable, sipSample, well):
    asvAbundances = list()
    for asv in zOTUs:
        if asv not in zOTUtable.keys():
            continue
        longID = sipSample.oneword_id+str(sipSample.plate)+well
        
        if longID+'redo' in zOTUtable[asv].keys() and longID in zOTUtable[asv].keys():
            asvAbundances.append(zOTUtable[asv][longID+'redo'])
            continue
        elif longID in zOTUtable[asv].keys():
            asvAbundances.append(zOTUtable[asv][longID])
    return asvAbundances

def main():
    df = read_excel("./data/fractionation/220817_Batch_122_Water_Yr_Summary_mod.xlsx", sheet_name='Summary', engine='openpyxl')

    sipData = parseSIPoutput(df)
    sipSamples = createSIPSampleClasses(SIPSample, sipData, './data/metadata/all_samples.csv')

    zOTUs, taxonomies = laodTaxonomyTable('./data/sequencing/amplicon_2022/taxonomy.tsv')
    zOTUtable = loadZOTUtable('./data/sequencing/amplicon_2022/zotutab.txt')

    with open('./processing/taxonomy_legend.csv', 'w') as f:
        for zOTU, taxonomy in zip(zOTUs, taxonomies):
            f.write(f'{taxonomy}\t{zOTU}\n')
    
    with open('./processing/zotu_table.csv', 'w') as f:
        f.write(f'ID\ttube\tFraction\tHabitat\tReplicate\tTime\tIsotope\tTrt.code\tDensity\tDNA.ng.fraction\tcopies.ul\t{chr(0x9).join(zOTUs)}\n')
        i = 1
        for sipSample in sipSamples:
            for j, (density, conc, well) in enumerate(zip(sipSample.densities, sipSample.concentrations, sipSample.wells)):
                if conc < 0:
                    conc = 0
                asvAbundances = getASVsPerSample(zOTUs, zOTUtable, sipSample, well)
                if len(asvAbundances) == 0:
                    continue
                totalAbundance = sum([float(x) for x in asvAbundances])
                relative_abundances = [str(float(x)/totalAbundance) for x in asvAbundances]
                tube = sipSample.sample_id.replace('50-80','D5').replace('-','_')
                f.write(f'{i}\t{tube}\t{j+1}\tD5\t{ABC_TO_NUMBERS[sipSample.replicate]}\t{sipSample.incubation_length}\t{sipSample.isotope}O\t{tube}_{j+1}\t{density}\t{conc*FINAL_SIP_VOLUME}\t{conc*FINAL_SIP_VOLUME}\t{chr(0x9).join(relative_abundances)}\n')
                i += 1

    with open('./processing/copies_per_g_soil.csv', 'w') as f:
        f.write(f'unique.tube\tg.wet.soil.extracted\ttotal.DNA.extracted.ng.ul\ttotal.DNA.extracted.ng\ttotal.DNA.extracted.ug\tug.DNA.g.DM.soil\tdry.wet\tg.dry.soil.extracted\tug.DNA.added.to.tube\tg.dry.soil.tube\n')
        for sipSample in sipSamples:
            tube = sipSample.sample_id.replace('50-80','D5').replace('-','_')
            totalDNAextracted = sipSample.concentration_DNA_extracted*sipSample.total_ul_DNA_extracted
            totalDNAextracted_ug = totalDNAextracted/1000
            totalDNAloaded_ug = sipSample.ul_DNA_SIP_loaded*sipSample.concentration_DNA_extracted/1000
            initial_water = sipSample.initial_GWC*sipSample.sample_weight
            final_wet_soil = sipSample.sample_weight + sipSample.h20_added
            final_GWC = ((initial_water + sipSample.h20_added)/(final_wet_soil))/100
            final_dry_soil = final_wet_soil - final_GWC*final_wet_soil
            final_dry_soil_extracted = sipSample.grams_soil_extracted - sipSample.grams_soil_extracted*final_GWC
            dry_soil_in_tube = final_dry_soil_extracted/(totalDNAextracted_ug/totalDNAloaded_ug)

            f.write(f'{tube}\t{sipSample.grams_soil_extracted}\t{sipSample.concentration_DNA_extracted}\t{totalDNAextracted}\t{totalDNAextracted_ug}\t{totalDNAextracted_ug/final_dry_soil_extracted}\t{final_GWC}\t{final_dry_soil_extracted}\t{totalDNAloaded_ug}\t{dry_soil_in_tube}\n')
if __name__ == "__main__":
    main()