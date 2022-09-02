#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from dataclasses import dataclass

data = pd.read_excel('data/220817_Batch_122_Water_Yr_Summary_orig.xlsx', 'Summary')

dbID_to_namedID = {2018: 'AA4_D5_16O-7',
                   2023: 'AB4_D5_16O-7',
                   2028: 'AC4_D5_16O-7',
                   2046: 'AA4_D5_16O-30',
                   2049: 'AB4_D5_16O-30',
                   2052: 'AC4_D5_16O-30',
                   2033: 'AA4_D5_18O-7',
                   2038: 'AB4_D5_18O-7',
                   2043: 'AC4_D5_18O-7',
                   2055: 'AA4_D5_18O-30',
                   2058: 'AB4_D5_18O-30',
                   2061: 'AC4_D5_18O-30'
                   }

isotopePairs = [(2018, 2033),
                (2023, 2038),
                (2028, 2043),
                (2046, 2055),
                (2049, 2058),
                (2052, 2061)]

density1 = data['Unnamed: 2']
density2 = data['Unnamed: 5']
density3 = data['Unnamed: 8']
density4 = data['Unnamed: 11']

concentration1 = data['Unnamed: 3']
concentration2 = data['Unnamed: 6']
concentration3 = data['Unnamed: 9']
concentration4 = data['Unnamed: 12']

@dataclass
class Sample:
    dbID: int
    name: str
    well_location: list
    density: list
    concentration: list

    def __post_init__(self):
        if not isinstance(self.dbID, int):
            raise ValueError("Parameter dbID must be an int!")
        if not isinstance(self.name, str):
            raise ValueError("Parameter dbID must be a str!")
        if not isinstance(self.well_location, list):
            raise ValueError("Parameter well_location must be a list!")
        if not isinstance(self.density, list):
            raise ValueError("Parameter density must be a list!")
        if not isinstance(self.concentration, list):
            raise ValueError("Parameter concentration must be a list!")
        self.area = [a*b for a,b in zip(self.density,self.concentration)]
        self.weighted_mean_density = sum(self.area)/sum(self.concentration)
        self.remainingDNAperFraction = [((abs(x)+x)/2)*37 for x in self.concentration]

samples = list()
samples.append(Sample(2018, dbID_to_namedID[2018], list(data[2018][3:22]), list(density1[3:22]), list(concentration1[3:22])))
samples.append(Sample(2023, dbID_to_namedID[2023], list(data[2023][3:22]), list(density2[3:22]), list(concentration2[3:22])))
samples.append(Sample(2028, dbID_to_namedID[2028], list(data[2028][3:22]), list(density3[3:22]), list(concentration3[3:22])))
samples.append(Sample(2046, dbID_to_namedID[2046], list(data[2046][3:22]), list(density4[3:22]), list(concentration4[3:22])))

samples.append(Sample(2049, dbID_to_namedID[2049], list(data[2018][31:50]), list(density1[31:50]), list(concentration1[31:50])))
samples.append(Sample(2052, dbID_to_namedID[2052], list(data[2023][31:50]), list(density2[31:50]), list(concentration2[31:50])))
samples.append(Sample(2033, dbID_to_namedID[2033], list(data[2028][31:48]), list(density3[31:48]), list(concentration3[31:48])))
samples.append(Sample(2038, dbID_to_namedID[2038], list(data[2046][31:50]), list(density4[31:50]), list(concentration4[31:50])))

samples.append(Sample(2043, dbID_to_namedID[2043], list(data[2018][58:77]), list(density1[58:77]), list(concentration1[58:77])))
samples.append(Sample(2055, dbID_to_namedID[2055], list(data[2023][58:77]), list(density2[58:77]), list(concentration2[58:77])))
samples.append(Sample(2058, dbID_to_namedID[2058], list(data[2028][58:77]), list(density3[58:77]), list(concentration3[58:77])))
samples.append(Sample(2061, dbID_to_namedID[2061], list(data[2046][58:77]), list(density4[58:77]), list(concentration4[58:77])))

def calcAtomPCT(wm16, isotopeShift):
    gcContent = (wm16-1.66)/0.098
    m16 = 307.691+0.496*gcContent
    m18 = ((isotopeShift/wm16)+1)*m16
    delM = m18-m16
    m18Max = 12.07747+m16
    atomPCTconst = 1-0.002000429
    return 100*((m18-m16)/(m18Max-m16))*atomPCTconst

def plotGraph(sample16, sample18, ax, title):
    density16 = sample16.density 
    density18 = sample18.density 
    conc16 = sample16.concentration 
    conc18 = sample18.concentration 
    mean16 = sample16.weighted_mean_density 
    mean18 = sample18.weighted_mean_density
    ax.set_title(title)
    ax.plot(density16, conc16, 'b-')
    ax.plot(density18, conc18, 'r-')
    ax.axvline(x = mean16, color = 'b', label = '16O')
    ax.axvline(x = mean18, color = 'r', label = '18O')
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

def splitFractions(sampleList, targetDNA=120):
    for sample in sampleList:
        heavy_chunk = splitsum(sample.remainingDNAperFraction, targetDNA)[0]
        light_chunk = splitsum(sample.remainingDNAperFraction[::-1], targetDNA)[0]
        mid_chunks = splitsum(sample.remainingDNAperFraction[len(heavy_chunk):-len(light_chunk)], targetDNA)
        sample.chunked_fractions = [heavy_chunk,*mid_chunks,light_chunk[::-1]]
    return sampleList

def chunkProperties(sample, chunk_ixes):
    chunk_ixes[0]=chunk_ixes[0]-1
    list_ixes=list(np.array(chunk_ixes).cumsum())
    chunked_fractions = lindexsplit(sample.remainingDNAperFraction, *list_ixes)[:-1]
    wells = lindexsplit(sample.well_location, *list_ixes)[:-1]
    sample.chunked_fractions = chunked_fractions
    sample.chunked_wells = wells
    return sample

def improve_SIP_bins(sampleList, passThrough=0):
    if passThrough+1 >= min([len(x.chunked_densities) for x in sampleList]):
        return sampleList
    modifiedChunks = False
    lowestFractionFromSet = min([min(x.chunked_densities[passThrough]) for x in sampleList])
    for sample in sampleList:
        remove_indices = list()
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


samplePairs = dict()
fig, axs = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(16,10))

for num,iPair in enumerate(isotopePairs):
    col = num%3
    row = round(num/6+0.1)
    iso16O = [item for item in samples if item.dbID==iPair[0]][0]
    iso18O = [item for item in samples if item.dbID==iPair[1]][0]
    isotope_shift = iso18O.weighted_mean_density-iso16O.weighted_mean_density
    atomPCTenrichment = calcAtomPCT(iso16O.weighted_mean_density, isotope_shift)
    print(f"{iso16O.dbID}:{iso16O.name}", f"{iso18O.dbID}:{iso18O.name}")
    print(atomPCTenrichment, isotope_shift, iso16O.weighted_mean_density, iso18O.weighted_mean_density)
    print(iso18O.remainingDNAperFraction)
    samplePairs[iPair] = {'isotope_shift': isotope_shift,
                          'atomPCTenrichment': atomPCTenrichment}
    incubation = 7
    if row == 1:
        incubation = 30

    plotGraph(iso16O, iso18O, axs[row,col], f'{iso16O.name[:6]} days: {incubation}  at% enrichment: {atomPCTenrichment:.2f}%')
    #print(len(iso16O.density), len(iso18O.density), len(iso16O.concentration), len(iso18O.concentration))

samples = splitFractions(samples, targetDNA=120)
iso18Osamples = list()
for sample in samples:
    if '_18O-' in sample.name:
        chunk_ixes = [len(x) for x in sample.chunked_fractions]
        chunk_ixes[0]=chunk_ixes[0]-1
        list_ixes=list(np.array(chunk_ixes).cumsum())
        density_chunks = lindexsplit(sample.density, *list_ixes)[:-1]
        wells = lindexsplit(sample.well_location, *list_ixes)[:-1]
        sample.chunked_densities = density_chunks
        sample.chunked_wells = wells
        iso18Osamples.append(sample)
        # print(sample.name)
        # print('\t'.join([' '.join(x) for x in wells]))
        # print('\t'.join([' '.join([f'{y:.3f}' for y in x]) for x in sample.chunked_densities]))
        # print('\t'.join([' '.join([f'{y:.2f}' for y in x]) for x in sample.chunked_fractions]))

improve_SIP_bins(iso18Osamples)

for sample in iso18Osamples:
    chunk_ixes = [len(x) for x in sample.chunked_densities]
    sample = chunkProperties(sample, chunk_ixes)
    print(sample.name)
    print('\t'.join([' '.join(x) for x in sample.chunked_wells]))
    print('\t'.join([' '.join([f'{y:.3f}' for y in x]) for x in sample.chunked_densities]))
    print('\t'.join([str(sum(x)) for x in sample.chunked_fractions]))

for ax in axs.flat:
    ax.set_xlim(1.6, 1.8)
    ax.set(xlabel='density (g/ml)', ylabel='DNA (ng/ul)')

for ax in axs.flat:
    ax.label_outer()

plt.plot()
#plt.savefig('figures/A4_D5_7-30_18Ocurves.svg', dpi=300)