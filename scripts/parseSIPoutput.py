#!/usr/bin/env python3


import sys, csv, argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from os import listdir
from itertools import combinations
from typing import List
from math import ceil

from scripts.waterYearSamplesInfrastructure import FractionatedSIPSample

EXCEL_LOCATION = './data/fractionation/'
SAMPLE_METADATA = '../all_samples.csv'

def parse_args(argument_list):
    parser = argparse.ArgumentParser(description='Parse SIP fractionation files, plot enrichment curves, and calculate fractionation bins.')
    parser.add_argument('-f', '--fractionation_files', help='Folder with fractionation output files in xlsx formtat.', default=EXCEL_LOCATION)
    parser.add_argument('-m', '--metadata_file', help='CSV file with sample metadata.', default=SAMPLE_METADATA)
    parser.add_argument('-o', '--output', help='Output folder for plots.', default='./figures/')
    parser.add_argument('-d', '--target_dna', help='Target DNA amount in ng for binned fractions.', default=100, type=int)
    parser.add_argument('-b', '--bins', help='Calculate fractionation bins and save to CSV file.', action='store_true')
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

def build_sample_objects(metadata_file, samples):
    fractionation_samples = list()
    with open(metadata_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            db_id = int(row['id'])
            if db_id not in samples.keys():
                continue
            if row['sample_id'][-2:] == '30':
                continue
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
                concentrations = samples[db_id][5],
                fraction_volumes = samples[db_id][6],
            )
            current_sample.batch = samples[db_id][7]
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
    return 100*((m18-m16)/(m18Max-m16))*atomPCTconst

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
    result,t = [[]],0
    for n in L:
        if t>S:
            r,v,t = (result,[n],n)
        else:
            r,v,t = (result[-1],n,t+n)
        r.append(v)
    return result

def splitFractions(sampleList, targetDNA=100):
    for sample in sampleList:
        heavy_chunk = splitsum(sample.remainingDNAperFraction, targetDNA)[0]
        light_chunk = splitsum(sample.remainingDNAperFraction[::-1], targetDNA)[0]
        mid_chunks = splitsum(sample.remainingDNAperFraction[len(heavy_chunk):-len(light_chunk)], targetDNA)
        sample.chunked_fractions = [heavy_chunk,*mid_chunks,light_chunk[::-1]]
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
            if density > lowestFractionFromSet:
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

def plotFractionationDataByIsotope(fr_samples16O, fr_samples18O=None):
    #Plot the fractionation data by isotope pair
    if fr_samples18O is not None:
        fig, axs = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(20,20))
    else:
        fig, axs = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(20,20))
        axs = [axs]
    cmap160 = cm.get_cmap('tab20', len(set([f'{x.batch}-{x.plate}' for x in fr_samples16O])))
    #make the colormap a named dictionary with keys batch and plate and values the color
    named_cmap16O = dict()
    for color, batch_plate in zip(cmap160.colors, set([f'{x.batch}-{x.plate}' for x in fr_samples16O])):
        named_cmap16O[batch_plate] = color

    if fr_samples18O is not None:
        cmap18O = cm.get_cmap('tab20', len(set([f'{x.batch}-{x.plate}' for x in fr_samples18O])))
        named_cmap18O = dict()
        for color, batch_plate in zip(cmap18O.colors, set([f'{x.batch}-{x.plate}' for x in fr_samples18O])):
            named_cmap18O[batch_plate] = color
        for i, fr_sample in enumerate(fr_samples18O):
            axs[1].plot(fr_sample.densities, fr_sample.concentrations, 
                        color=named_cmap18O[f'{fr_sample.batch}-{fr_sample.plate}'], label=f'{fr_sample.id}-{fr_sample.batch}: {fr_sample.weighted_mean_density:.3f}')
        lgd2 = axs[1].legend(loc='center left', bbox_to_anchor=(1, 0.5))
        axs[1].set_title('18O')

    for i, fr_sample in enumerate(fr_samples16O):
        axs[0].plot(fr_sample.densities, fr_sample.concentrations, 
                    color=named_cmap16O[f'{fr_sample.batch}-{fr_sample.plate}'], label=f'{fr_sample.id}-{fr_sample.dna_yield:.1f}: {fr_sample.weighted_mean_density:.3f}')
    lgd = axs[0].legend(loc='center left', bbox_to_anchor=(1, 0.5))
    axs[0].set_title('16O')
   
    for ax in axs:
        ax.set_xlim(1.65, 1.77)
        ax.set(xlabel='density (g/ml)', ylabel='DNA (ng/ul)')
    plt.tight_layout()
    plt.savefig('figures/fractionation_data_both.png', dpi=300)
    plt.cla()
    plt.clf()
    return True


def main(commandline_arguments):

    args = parse_args(commandline_arguments)

    samples = extract_fractionation_samples_from_excel_files(args.fractionation_files)

    fractionation_samples = build_sample_objects(args.metadata_file, samples)

    #Find isotope pairs from the fractionation samples
    isotopePairs = findIsotopePairs(fractionation_samples)
    unpairedSamples = [item for item in fractionation_samples if item.id not in [item.id for sublist in isotopePairs for item in sublist]]

    fr_samples_16O = [item for item in fractionation_samples if item.experiment == '16O-7']
    fr_samples_18O = [item for item in fractionation_samples if item.experiment == '18O-7']

    #Plot fractionation data by isotope
    plotFractionationDataByIsotope(fr_samples_16O, fr_samples_18O)


    #Plotting
    fig, axs = plt.subplots(7, int(ceil(len(isotopePairs)/7)), sharex=True, sharey=True, figsize=(30,20))

    print("Unpaired samples:")
    for sample in unpairedSamples:
        print(str(sample), sample.id)


    for isotope_pair, axis in zip(isotopePairs, axs.flatten()):

        iso16O = isotope_pair[0]
        iso18O = isotope_pair[1]
        isotope_shift = iso18O.weighted_mean_density-iso16O.weighted_mean_density
        atomPCTenrichment = calcAtomPCT(iso16O.weighted_mean_density, isotope_shift)
        print(f'{str(iso16O)} ape: {atomPCTenrichment:.2f}%')
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
    plt.savefig('figures/all_samples_v3.svg', dpi=300)

    if not args.bins:
        sys.exit(0)

    #Binning first pass
    samples = splitFractions(fractionation_samples, targetDNA=args.target_dna)
    iso18Osamples = list()
    for sample in samples:
        if '_18O-' in str(sample) and sample.dna_yield > 10:
            chunk_ixes = [len(x) for x in sample.chunked_fractions]
            sample = chunkProperties(sample, chunk_ixes)
            iso18Osamples.append(sample)

    #Binning second pass
    improve_SIP_bins(iso18Osamples)

    #Bin output
    for sample in iso18Osamples:
        chunk_ixes = [len(x) for x in sample.chunked_densities]
        sample = chunkProperties(sample, chunk_ixes)
        # print(f'{str(sample)}\t','\t'.join([' '.join(x) for x in sample.chunked_wells]))
        # print(f'{str(sample)}\t','\t'.join([' '.join([f'{y:.3f}' for y in x]) for x in sample.chunked_densities]))
        # print(f'{str(sample)}\t','\t'.join([str(sum(x)) for x in sample.chunked_fractions]))

        #print(f'{str(sample)}\t','\t'.join(['-'.join(x) for x in sample.chunked_wells]))
        #print(f'{str(sample)}\t','\t'.join([f'{np.average(x):.3f} {np.std(x):.3f}' for x in sample.chunked_densities]))
        #print(f'{str(sample)}\t','\t'.join([' '.join([f'{y:.3f}' for y in x]) for x in sample.chunked_densities]))
        #print(f'{str(sample)}\t','\t'.join([str(sum(x)) for x in sample.chunked_fractions]))

        #Print only the min max density for each chunk
        print(f'{str(sample)}\t','\t'.join([f'{max(x):.3f}\t{min(x):.3f}' for x in sample.chunked_densities]))

if __name__ == '__main__':
    main(sys.argv[1:])