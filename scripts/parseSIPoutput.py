#!/usr/bin/env python3

import sys, csv, argparse, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from os import listdir
from itertools import combinations
from math import ceil
from scipy.stats import linregress

#So that tests run from the root directory
sys.path.append('./scripts/')

from waterYearSamplesInfrastructure import FractionatedSIPSample

EXCEL_LOCATION = './data/fractionation/test_fract_adjust/'
PICO_LOCATION = './data/fractionation/Quant/'
SAMPLE_METADATA = '../all_samples.csv'

def parse_args(argument_list):
    parser = argparse.ArgumentParser(description='Parse SIP fractionation files, plot enrichment curves, and calculate fractionation bins.')
    parser.add_argument('-f', '--fractionation_files', help='Folder with fractionation output files in xlsx formtat.', default=EXCEL_LOCATION)
    parser.add_argument('-p', '--pico_quantification_file', help='Folder with picogreen results in xlsx formtat.', default=PICO_LOCATION)
    parser.add_argument('-m', '--metadata_file', help='CSV file with sample metadata.', default=SAMPLE_METADATA)
    parser.add_argument('-o', '--output', help='Output folder for plots.', default='./figures/')
    parser.add_argument('-d', '--target_dna', help='Target DNA amount in ng for binned fractions.', default=100, type=int)
    parser.add_argument('-b', '--bins', help='Calculate fractionation bins and save to CSV file.', action='store_true')
    parser.add_argument('-lc', '--light_cutoff', help='Light density cutoff for fractionation wells.', default=1.62, type=float)
    parser.add_argument('-hc', '--heavy_cutoff', help='Heavy density cutoff for fractionation wells.', default=1.76, type=float)
    parser.add_argument('-dc', '--dna_cutoff', help='Total DNA amount cutoff for fractionation wells.', default=0, type=float)
    parser.add_argument('-l', '--linear_fix', help='Fix densities using linear regression.', action='store_true')
    return parser.parse_args(argument_list)

def extract_fractionation_samples_from_excel_files(EXCEL_LOCATION):
    samples = dict()
    duplicate_samples = list()
    for file in listdir(EXCEL_LOCATION):
        if not file.endswith(".xlsx"):
            continue
        data = pd.read_excel(EXCEL_LOCATION+file, 'Summary')
        toc = pd.read_excel(EXCEL_LOCATION+file, 'Table of Contents', header=3)

        #Get the column name from toc.columns starting with 'Final Volume'
        res = [key for key in toc.columns if 'Final Volume' in key]
        volName = res[0]
        for dbID, plate, tube, dna_yield, wells_volume in zip(toc['Sample ID'], toc['Plate Label'], 
            toc['Tube Letter'], toc['Percent DNA Recovered'], toc[volName]):
            if np.isnan(dbID):
                break
            if dbID not in samples:
                samples[dbID] = (plate, tube[-1], dna_yield, list(data[dbID][3:22]), list(data[f'{dbID}-density'][3:22]), list(data[f'{dbID}-conc'][3:22]), [wells_volume]*19, file.split(' ')[2])
            else:
                duplicate_samples.append((dbID, file, samples[dbID][-1]))
                #raise ValueError(f"Duplicate dbID {dbID} found in {file} and batch {samples[dbID][-1]}!")
    if len(duplicate_samples) > 0:
        print(f"Duplicate samples found:")
        for sample in duplicate_samples:
            print(f"{sample[0]} in {sample[1]} and batch {sample[2]}")
        raise ValueError("Duplicate samples found!")
    return samples

def load_picogreen_quantification_data(PICO_LOCATION):
    samples = dict()
    for file in listdir(PICO_LOCATION):
        if not file.endswith(".xlsx"):
            continue
        data = pd.read_excel(PICO_LOCATION+file, 'Consolidated')
        batch, plate = file.split()[2].replace('batch',''), file.split()[3].replace('.xlsx','')
        wells, ave_yields, stdevs = list(data['Fractionation well']), list(data['Average yield']), list(data['stdev'])
        dna_wells = [wells[i * 22:(i + 1) * 22] for i in range((len(wells) + 22 - 1) // 22 )] 
        dna_yields = [ave_yields[i * 22:(i + 1) * 22] for i in range((len(ave_yields) + 22 - 1) // 22 )] 
        dna_stdevs = [stdevs[i * 22:(i + 1) * 22] for i in range((len(stdevs) + 22 - 1) // 22 )] 
        if batch not in samples:
            samples[batch] = {}
        samples[batch] = {**samples[batch], **{plate: (dna_wells, dna_yields, dna_stdevs)}}
    return samples

def build_sample_objects(metadata_file, samples, picogreen_data):
    fractionation_samples = list()
    tube_to_ix = {key: value for key, value in zip(['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P'],[0,1,2,3]*4)}
    with open(metadata_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            db_id = int(row['id'])
            if db_id not in samples.keys():
                continue
            if row['sample_id'][-2:] == '30':
                continue
            if picogreen_data[samples[db_id][7]][samples[db_id][0]][0][tube_to_ix[samples[db_id][1]]][1:20] != samples[db_id][3]:
                print(f"Mismatched wells for batch {samples[db_id][7]}, plate {samples[db_id][0]}, sample {db_id}!")
            current_sample = FractionatedSIPSample(
                id = db_id,
                sampling_site = row['location'][0],
                replicate = row['location'][1],
                sampling_week = row['sample_id'][2],
                depth = tuple(row['depth'].split('-')),
                depth_str = row['sample_id'][4:6],
                date = row['date_collected'],
                metadata_file_path='../all_samples.csv',
                experiment = row['sample_id'][7:],
                ul_DNA_SIP_loaded = row['DNA loaded for SIP (ul)'],
                plate = samples[db_id][0],
                tube = samples[db_id][1],
                dna_yield = samples[db_id][2],
                wells = samples[db_id][3],
                densities = samples[db_id][4],
                concentrations = [(abs(x)+x)/2 for x in samples[db_id][5]],
                fraction_volumes = samples[db_id][6],
            )
            current_sample.conc_stdevs = picogreen_data[samples[db_id][7]][samples[db_id][0]][2][tube_to_ix[samples[db_id][1]]][1:20]
            current_sample.batch = samples[db_id][7]
            current_sample.extraction_id = row['Extraction ID']
            fractionation_samples.append(current_sample)
    return fractionation_samples

def findIsotopePairs(fractionation_samples):
    isotope_pairs = list()
    for sample1, sample2 in combinations(fractionation_samples, 2):
        if str(sample1)[0:6] == str(sample2)[0:6]:
            if '16O' in str(sample1) and '18O' in str(sample2):
                isotope_pairs.append((sample1, sample2))
            elif '16O' in str(sample2) and '18O' in str(sample1):
                isotope_pairs.append((sample2, sample1))
            else:
                raise ValueError(f"Sample {sample1} and {sample2} are not isotope pairs!")
    return isotope_pairs

def calcAtomPCT(wm16, isotopeShift):
    gcContent = (wm16-1.66)/0.098
    m16 = 307.691+0.496*gcContent
    m18 = ((isotopeShift/wm16)+1)*m16
    delM = m18-m16
    m18Max = 12.07747+m16
    atomPCTconst = 1-0.002000429
    return 100*((delM)/(m18Max-m16))*atomPCTconst

def plotGraph(sample18, ax, title, sample16=None):

    density18 = sample18.densities 
    conc18 = sample18.concentrations 
    mean18 = sample18.weighted_mean_density
    ax.set_title(title)
    if sample16:
        density16 = sample16.densities
        conc16 = sample16.concentrations
        mean16 = sample16.weighted_mean_density
        ax.plot(density16, conc16, 'b-')
        ax.axvline(x = mean16, color = 'b',  alpha=0.5, label = '16O')
        ax.axvline(x = mean18, color = 'r', alpha=0.5, label = '18O')
        ax.plot(density18, conc18, 'r-')
        if not hasattr(sample18, 'chunked_densities'):
            return True
    cm = plt.get_cmap('Dark2')
    for i,(x,y) in enumerate(zip(sample18.chunked_densities,sample18.chunked_concentrations)):
        if i+1 < len(sample18.chunked_densities):
            x.append(sample18.chunked_densities[i+1][0])
            #x.append((sample18.chunked_densities[i+1][0]+sample18.chunked_densities[i][-1])/2)
            y.append(sample18.chunked_concentrations[i+1][0])
            #y.append((sample18.chunked_concentrations[i+1][0]+sample18.chunked_concentrations[i][-1])/2)
        ax.fill_between(x=x,
                        y1=y,
                        color=cm.colors[i],
                        alpha=0.5)
    ax.plot(density18, conc18, 'r-')
    return True

def lindexsplit(some_list, *args):
    # Checks to see if any extra arguments were passed. If so,
    # prepend the 0th index and append the final index of the 
    # passed list. This saves from having to check for the beginning
    # and end of args in the for-loop. Also, increment each value in 
    # args to get the desired behavior.
    if args:
        args = (0,) + tuple(data+1 for data in args) + (len(some_list)+1,)

    # For a little more brevity, here is the list comprehension of the following
    # statements:
    #    return [some_list[start:end] for start, end in zip(args, args[1:])]
    my_list = []
    for start, end in zip(args, args[1:]):
        my_list.append(some_list[start:end])
    return my_list

def splitsum(L,S):
    '''Splits a list of numbers into sublists whose sums are just about bigger than S'''
    result,t = [[]],0
    for n in L:
        if t>S:
            r,v,t = (result,[n],n)
        else:
            r,v,t = (result[-1],n,t+n)
        r.append(v)
    return result

def extract_heavy_samples(samples, isotope_string_to_match='_18O-'):
    heavy_samples = list()
    for sample in samples:
        if isotope_string_to_match in str(sample) and sample.dna_yield > 10:
            heavy_samples.append(sample)
    return heavy_samples

def mergeFractions(sampleList, targetDNA=100):
    for sample in sampleList:
        heavy_chunk = splitsum(sample.remainingDNAperFraction, targetDNA)[0]
        light_chunk = splitsum(sample.remainingDNAperFraction[::-1], targetDNA)[0]
        mid_chunks = splitsum(sample.remainingDNAperFraction[len(heavy_chunk):-len(light_chunk)], targetDNA)
        sample.chunked_fractions = [heavy_chunk,*mid_chunks,light_chunk[::-1]]
        chunk_ixes = [len(x) for x in sample.chunked_fractions]
        sample = chunkProperties(sample, chunk_ixes)
        # chunks = splitsum(sample.remainingDNAperFraction, targetDNA)
        # sample.chunked_fractions = chunks
    return sampleList

def chunkProperties(sample, chunk_ixes):
    '''Split the dnaFractions, wells, and densities in the same chunks as the provided indexes'''
    chunk_ixes[0]=chunk_ixes[0]-1
    list_ixes=list(np.array(chunk_ixes).cumsum())
    chunked_fractions = lindexsplit(sample.remainingDNAperFraction, *list_ixes)[:-1]
    chunked_concentrations = lindexsplit(sample.concentrations, *list_ixes)[:-1]
    chunked_wells = lindexsplit(sample.wells, *list_ixes)[:-1]
    chunked_densities = lindexsplit(sample.densities, *list_ixes)[:-1]
    sample.chunked_densities = chunked_densities
    sample.chunked_concentrations = chunked_concentrations
    sample.chunked_fractions = chunked_fractions
    sample.chunked_wells = chunked_wells
    return sample

def improve_SIP_bins(sampleList, passThrough=0):
    if passThrough+1 >= max([len(x.chunked_densities) for x in sampleList]):
        return sampleList
    modifiedChunks = False
    lowestFractionFromSet = min([min(x.chunked_densities[passThrough]) for x in sampleList])
    for sample in sampleList:
        remove_indices = list()
        if passThrough+1 >= len(sample.chunked_densities):
            continue
        for i,density in enumerate(sample.chunked_densities[passThrough+1]):
            if density >= lowestFractionFromSet:
                modifiedChunks = True
                sample.chunked_densities[passThrough].append(density)
                remove_indices.append(i)
        if len(remove_indices) > 0:
            sample.chunked_densities[passThrough+1] = list(np.delete(np.array(sample.chunked_densities[passThrough+1]),remove_indices))
            sample.chunked_densities = [ele for ele in sample.chunked_densities if ele != []]
    if modifiedChunks:
        improve_SIP_bins(sampleList, passThrough=passThrough)
    else:
        improve_SIP_bins(sampleList, passThrough=passThrough+1)

def merge_lightest_bins(sample):
    '''Merge the two lightest chunks into one'''
    sample.chunked_densities[-2] = sample.chunked_densities[-2] + sample.chunked_densities[-1]
    sample.chunked_densities = sample.chunked_densities[:-1]
    sample.chunked_concentrations[-2] = sample.chunked_concentrations[-2] + sample.chunked_concentrations[-1]
    sample.chunked_concentrations = sample.chunked_concentrations[:-1]
    sample.chunked_fractions[-2] = sample.chunked_fractions[-2] + sample.chunked_fractions[-1]
    sample.chunked_fractions = sample.chunked_fractions[:-1]
    sample.chunked_wells[-2] = sample.chunked_wells[-2] + sample.chunked_wells[-1]
    sample.chunked_wells = sample.chunked_wells[:-1]
    return sample

def identify_density_thresholds_between_bins(sampleList):
    '''Identify the density thresholds between bins'''
    max_min_densities, threshold_densities = list(), list()
    for sample in sampleList:
        max_min_densities.append([[max(x),min(x)] for x in sample.chunked_densities])
    
    #Identify the minimum and maximum densities across all samples
    max_min_densities = list(zip(*max_min_densities))
    range_densities = [([max(list(zip(*x))[0]),min(list(zip(*x))[1])]) for x in max_min_densities]
    for i,x in enumerate(range_densities):
        if i+1 < len(range_densities):
            threshold_densities.append(np.average([x[1], range_densities[i+1][0]]))
    return threshold_densities

def plotFractionationDataByIsotope(fr_samples16O, fr_samples18O=None):
    #Plot the fractionation data by isotope pair
    if fr_samples18O is not None:
        fig, axs = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(20,20))
    else:
        fig, axs = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(20,20))
        axs = [axs]
    cmap160 = cm.get_cmap('tab20', len(set([x.batch for x in fr_samples16O])))
    #make the colormap a named dictionary with keys batch and plate and values the color
    named_cmap16O = dict()
    for color, batch_plate in zip(cmap160.colors, set([x.batch for x in fr_samples16O])):
        named_cmap16O[batch_plate] = color

    if fr_samples18O is not None:
        cmap18O = cm.get_cmap('tab20', len(set([x.batch for x in fr_samples18O])))
        named_cmap18O = dict()
        for color, batch_plate in zip(cmap18O.colors, set([x.batch for x in fr_samples18O])):
            named_cmap18O[batch_plate] = color
        for i, fr_sample in enumerate(fr_samples18O):
            axs[1].plot(fr_sample.densities, fr_sample.concentrations, 
                        color=named_cmap18O[fr_sample.batch], label=f'{fr_sample.id}-{fr_sample.batch}: {fr_sample.weighted_mean_density:.3f}')
        lgd2 = axs[1].legend(loc='center left', bbox_to_anchor=(1, 0.5))
        axs[1].set_title('18O')

    for i, fr_sample in enumerate(fr_samples16O):
        axs[0].plot(fr_sample.densities, fr_sample.concentrations, 
                    color=named_cmap16O[fr_sample.batch], label=f'{fr_sample.id}-{fr_sample.batch}: {fr_sample.weighted_mean_density:.3f}')
    lgd = axs[0].legend(loc='center left', bbox_to_anchor=(1, 0.5))
    axs[0].set_title('16O')
   
    for ax in axs:
        ax.set_xlim(1.65, 1.77)
        ax.set(xlabel='density (g/ml)', ylabel='DNA (ng/ul)')
    plt.tight_layout()
    plt.savefig('figures/fractionation_data_both_fixed.png', dpi=300)
    plt.cla()
    plt.clf()
    return True

def plotAllSampleCurves(isotopePairs, loc='figures/all_samples_v6_fixed.svg'):
    fig, axs = plt.subplots(7, int(ceil(len(isotopePairs)/7)), sharex=True, sharey=True, figsize=(40,30))
    for isotope_pair, axis in zip(isotopePairs, axs.flatten()):
        iso16O = isotope_pair[0]
        iso18O = isotope_pair[1]
        isotope_shift = iso18O.weighted_mean_density-iso16O.weighted_mean_density
        atomPCTenrichment = calcAtomPCT(iso16O.weighted_mean_density, isotope_shift)
        #print(f'{str(iso16O)} ape: {atomPCTenrichment:.2f}%')
        plotGraph(iso18O, 
                  axis, 
                  f'{str(iso16O)[:6]}({iso16O.id}:{iso16O.dna_yield:.0f},{iso18O.id}:{iso18O.dna_yield:.0f}){atomPCTenrichment:.2f}',
                  sample16=iso16O)

    for ax in axs.flat:
        ax.set_xlim(1.65, 1.77)
        ax.set(xlabel='density (g/ml)', ylabel='DNA (ng/ul)')

    for ax in axs.flat:
       ax.label_outer()

    plt.tight_layout()
    plt.plot()
    plt.savefig(loc, dpi=300)
    return True

def calculate_density_linear_regression(densities, plot=False):
    '''Fit a line through the lists of densities, each list being an independent linear regression'''
    data = np.array(densities)
    y = data.flatten()
    x_axis = np.linspace(1, len(densities[0]), len(densities[0]))
    x = [x_axis] * len(densities)
    x = np.array(x).flatten()
    slope, intercept, r_value, p_value, std_err = linregress(x, y)

    if plot:
        colors = np.repeat(np.arange(len(densities)), len(densities[0]))
        plt.scatter(x, y, c=colors, cmap='viridis')
        plt.plot(x, slope*x + intercept, color='black')
        plt.xlabel('Well Number')
        plt.ylabel('Density')
        plt.title('Linear Regression of Multiple Sets of Linear Data')
        plt.savefig('figures/linear_regression.png', dpi=300)
        plt.cla()
        plt.clf()

    return slope, intercept, r_value, p_value, std_err

def binLightSample(light_sample:FractionatedSIPSample, thresholds:list):
    '''Chunk the properties of the light sample based on density thresholds'''
    chunked_densities = list(np.split(light_sample.densities,np.searchsorted(-np.array(light_sample.densities), -np.array(thresholds))))
    chunk_ixes = [len(x) for x in chunked_densities]
    light_sample = chunkProperties(light_sample, chunk_ixes)
    return light_sample

def print_wells_for_binning(isotopePairs):
    '''Prints batch, plate and wells of the bins for each sample
    and the name of the sample based on extraction id and fraction number'''
    sampleList = list()
    for isotope_pair in isotopePairs:
        sampleList.append(isotope_pair[0])
        sampleList.append(isotope_pair[1])
    #Order samples by batch and plate
    sampleList.sort(key=lambda x: (x.batch, x.plate))
    dict_for_print = dict()
    for sample in sampleList:
        if sample.batch not in dict_for_print.keys():
            dict_for_print[sample.batch] = dict()
        if sample.plate not in dict_for_print[sample.batch].keys():
            dict_for_print[sample.batch][sample.plate] = list()
        dict_for_print[sample.batch][sample.plate].append(sample)
    
    for batch in sorted(dict_for_print.keys()):
        print(f'Batch {batch}:')
        for plate in sorted(dict_for_print[batch].keys()):
            print(f'\tPlate {plate}:')
            for sample in dict_for_print[batch][plate]:
                print(f'\t\tOld ID {sample.id}')
                for i, wells in enumerate(sample.chunked_wells):
                    print(f'\t\t\t{sample.extraction_id}_{i+1}: {" ".join(wells)}')

def main(commandline_arguments):

    args = parse_args(commandline_arguments)

    samples = extract_fractionation_samples_from_excel_files(args.fractionation_files)

    pico_quantification_stdevs = load_picogreen_quantification_data(args.pico_quantification_file)

    fractionation_samples = build_sample_objects(args.metadata_file, samples, pico_quantification_stdevs)

    #Print samples with high stdevs for the concentration measurements
    #print('Samples with high stdevs for the concentration measurements:')
    for sample in fractionation_samples:
        well_conc_stdev = list(zip(sample.wells, sample.concentrations, sample.conc_stdevs))
        well_conc_stdev_warn = [item for item in well_conc_stdev if item[2]*2 > item[1] and item[1] > 0.1]
        well_conc_stdev_nill = [item for item in well_conc_stdev if item[2]*2 > item[1] and 0.1 > item[1] > 0]
        for item in well_conc_stdev_warn:
            print(f'{sample.id}-{sample.batch}-{sample.plate}-{item[0]}: {item[1]:.3f}+-{item[2]:.3f}')
        for item in well_conc_stdev_nill:
            #output warnings to a file
            #with open('warnings.txt', 'a') as f:
            #    f.write(f'Setting to 0: {sample.id}-{sample.batch}-{sample.plate}-{item[0]}: {item[1]:.3f}+-{item[2]:.3f}\n')
            sample.concentrations[sample.wells.index(item[0])] = 0
            sample.remainingDNAperFraction[sample.wells.index(item[0])] = 0
            sample.area[sample.wells.index(item[0])] = 0


    #Find isotope pairs from the fractionation samples
    isotopePairs = findIsotopePairs(fractionation_samples)
    unpairedSamples = [item for item in fractionation_samples if item.id not in [item.id for sublist in isotopePairs for item in sublist]]

    fr_samples_16O = [item for item in fractionation_samples if item.experiment == '16O-7']
    fr_samples_18O = [item for item in fractionation_samples if item.experiment == '18O-7']

    if args.linear_fix:
        good_16O_densities = list()
        for sample in fr_samples_16O:
            if sample.batch == '144':
                continue
            if sample.batch == '135' and sample.plate == 'IJKL':
                continue
            good_16O_densities.append(sample.densities[:-2])

        slope, intercept, r_value, p_value, std_err = calculate_density_linear_regression(good_16O_densities, plot=True)

        bad_16O_samples = [x for x in fr_samples_16O if x.batch == '144' or x.batch == '135']
        for sample in bad_16O_samples:
            bad_slope, bad_intercept, bad_r_value, bad_p_value, bad_std_err = calculate_density_linear_regression([sample.densities[:-2]])
            print(f'{sample.id}, {sample.batch}, {sample.plate} Slope difference: {slope-bad_slope} Intercept difference: {intercept-bad_intercept}')

    #Plot fractionation data by isotope
    plotFractionationDataByIsotope(fr_samples_16O, fr_samples_18O)

    if len(unpairedSamples) > 0:
        print("Unpaired samples:")
        for sample in unpairedSamples:
            print(str(sample), sample.id)
    
    if not args.bins:
        plotAllSampleCurves(isotopePairs)
        sys.exit(0)

    # Remove fractions below the cutoff threshold
    for sample in fractionation_samples:
        sample.removeFractionsBelowDensityCutoff(args.light_cutoff)
        sample.removeFractionsAboveDensityCutoff(args.heavy_cutoff)
        #sample.removeFractionsBelowDNACutoff(args.dna_cutoff)

    #Binning first pass
    iso18Osamples = extract_heavy_samples(fractionation_samples)
    binned_18O_Samples = mergeFractions(iso18Osamples, targetDNA=args.target_dna)
    

    #Binning second pass
    improve_SIP_bins(binned_18O_Samples)

    #Merge lightest fraction bins if bin number is above 5
    for sample in binned_18O_Samples:
        if len(sample.chunked_densities) > 5:
            sample = merge_lightest_bins(sample)

    isotopePairs = findIsotopePairs(fractionation_samples)

    #Rechunk after bin improvements
    for sample in binned_18O_Samples:
        chunk_ixes = [len(x) for x in sample.chunked_densities]
        sample = chunkProperties(sample, chunk_ixes)

    #Bin the light samples based on the heavy samples
    heavy_density_thresholds = identify_density_thresholds_between_bins(binned_18O_Samples)
    print('Density thresholds:\t'+'\t'.join([f'{x:.5f}' for x in heavy_density_thresholds]))
    for sample16, sample18 in isotopePairs:
        sample16 = binLightSample(sample16, heavy_density_thresholds)
        conc_16O = ['-'.join([f'{y:.3f}' for y in x]) for x in sample16.chunked_concentrations]
        conc_18O = ['-'.join([f'{y:.3f}' for y in x]) for x in sample18.chunked_concentrations]
        dens_16O = [f'{max(x, default=0):.5f}:{min(x, default=0):.5f}' for x in sample16.chunked_densities]
        #dens_16O = ['-'.join([f'{y}' for y in x]) for x in sample16.chunked_densities]
        dens_18O = [f'{max(x, default=0):.5f}:{min(x, default=0):.5f}' for x in sample18.chunked_densities]
        #dens_18O = ['-'.join([f'{y}' for y in x]) for x in sample18.chunked_densities]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            average_conc_16O = [f'{np.average(x):.3f}' for x in sample16.chunked_concentrations]
            average_conc_18O = [f'{np.average(x):.3f}' for x in sample18.chunked_concentrations]
        bins_with_info_16O = ['_'.join(x) for x in list(zip(average_conc_16O, dens_16O, conc_16O))]
        bins_with_info_18O = ['_'.join(x) for x in list(zip(average_conc_18O, dens_18O, conc_18O))]
        #print(f'{str(sample16)}-{sample16.id}-{sample16.batch}-{sample16.plate}\t','\t'.join(bins_with_info_16O))
        #print(f'{str(sample18)}-{sample18.id}-{sample18.batch}-{sample18.plate}\t','\t'.join(bins_with_info_18O))

    print_wells_for_binning(isotopePairs)

    #Bin output
    #for sample in binned_18O_Samples:
        # print(f'{str(sample)}\t','\t'.join([' '.join(x) for x in sample.chunked_wells]))
        # print(f'{str(sample)}\t','\t'.join([' '.join([f'{y:.3f}' for y in x]) for x in sample.chunked_densities]))
        # print(f'{str(sample)}\t','\t'.join([str(sum(x)) for x in sample.chunked_fractions]))

        #print(f'{str(sample)}\t','\t'.join(['-'.join(x) for x in sample.chunked_wells]))
        #print(f'{str(sample)}\t','\t'.join([f'{np.average(x):.3f} {np.std(x):.3f}' for x in sample.chunked_densities]))
        #print(f'{str(sample)}\t','\t'.join([' '.join([f'{y:.3f}' for y in x]) for x in sample.chunked_densities]))
        #print(f'{str(sample)}\t','\t'.join([str(sum(x)) for x in sample.chunked_fractions]))

        #Print only the min max density for each chunk
        #print(f'{str(sample)}\t','\t'.join(['-'.join(x) for x in sample.chunked_wells]))
        #print(f'{str(sample)}\t','\t'.join([f'{max(x):.4f} {min(x):.4f}' for x in sample.chunked_densities]))
        #print('\t\t', '\t'.join([f'{sum(x):.4f}' for x in sample.chunked_fractions]))
        #print('\t\t', '\t'.join([f'{np.average(x):.4f}' for x in sample.chunked_concentrations]))
        #print('\t\t', '\t'.join([f'{np.average(x):.4f}' for x in sample.chunked_densities]))
        #print('\t\t', '\t'.join([f'{float(sum(x)/(len(x)*37)):.4f}' for x in sample.chunked_fractions]))
        #print('\t\t', '\t'.join(['-'.join(f'{y:.1f}' for y in x) for x in sample.chunked_fractions]))
        #print('\t\t', '\t'.join(['-'.join(f'{y:.1f}' for y in x) for x in sample.chunked_concentrations]))
        #print(f'{str(sample)}\t','\t'.join([str(sum(x)) for x in sample.chunked_fractions]))
    plotAllSampleCurves(isotopePairs, 'figures/all_samples_bins_v6_fixed.svg')

if __name__ == '__main__':
    main(sys.argv[1:])