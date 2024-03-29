#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from dataclasses import dataclass

from ..scripts.batch_125 import dbID_to_namedID, isotopePairs

excel_location = 'data/fractionation/221109 Batch 129 Water Year Summary.xlsx'
data = pd.read_excel(excel_location, 'Summary')

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
        self.remainingDNAperFraction = [((abs(x)+x)/2)*35 for x in self.concentration]

samples = list()

for dbID, name in dbID_to_namedID.items():
    samples.append(Sample(dbID, name, list(data[dbID][3:22]), list(data[f'{dbID}-density'][3:22]), list(data[f'{dbID}-conc'][3:22])))


def calcAtomPCT(wm16, isotopeShift):
    gcContent = (wm16-1.66)/0.098
    m16 = 307.691+0.496*gcContent
    m18 = ((isotopeShift/wm16)+1)*m16
    delM = m18-m16
    m18Max = 12.07747+m16
    atomPCTconst = 1-0.002000429
    return 100*((m18-m16)/(m18Max-m16))*atomPCTconst

def plotGraph(sample18, ax, title, sample16=None):

    density18 = sample18.density 
    conc18 = sample18.concentration 
    mean18 = sample18.weighted_mean_density
    ax.set_title(title)
    if sample16:
        density16 = sample16.density
        conc16 = sample16.concentration 
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

#Change to 100ng
def splitFractions(sampleList, targetDNA=120):
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
    chunked_concentrations = lindexsplit(sample.concentration, *list_ixes)[:-1]
    chunked_wells = lindexsplit(sample.well_location, *list_ixes)[:-1]
    chunked_densities = lindexsplit(sample.density, *list_ixes)[:-1]
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

#Binning first pass
samples = splitFractions(samples, targetDNA=100)
iso18Osamples = list()
for sample in samples:
    if '_18O-' in sample.name:
        chunk_ixes = [len(x) for x in sample.chunked_fractions]
        sample = chunkProperties(sample, chunk_ixes)
        iso18Osamples.append(sample)

#Binning second pass
#improve_SIP_bins(iso18Osamples)

#Bin output
for sample in iso18Osamples:
    chunk_ixes = [len(x) for x in sample.chunked_densities]
    sample = chunkProperties(sample, chunk_ixes)
    # print(f'{sample.name} \t','\t'.join([' '.join(x) for x in sample.chunked_wells]))
    # print(f'{sample.name} \t','\t'.join([' '.join([f'{y:.3f}' for y in x]) for x in sample.chunked_densities]))
    # print(f'{sample.name} \t','\t'.join([str(sum(x)) for x in sample.chunked_fractions]))

    print(f'{sample.name} \t','\t'.join(['-'.join(x) for x in sample.chunked_wells]))
    #print(f'{sample.name} \t','\t'.join([f'{np.average(x):.3f} {np.std(x):.3f}' for x in sample.chunked_densities]))
    print(f'{sample.name} \t','\t'.join([' '.join([f'{y:.3f}' for y in x]) for x in sample.chunked_densities]))
    print(f'{sample.name} \t','\t'.join([str(sum(x)) for x in sample.chunked_fractions]))

#Plotting
samplePairs = dict()
fig, axs = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(16,10))
for num,iPair in enumerate(isotopePairs):
    col = num%3
    row = round(num/6+0.1)
    iso16O = [item for item in samples if item.dbID==iPair[0]][0]
    iso18O = [item for item in samples if item.dbID==iPair[1]][0]
    isotope_shift = iso18O.weighted_mean_density-iso16O.weighted_mean_density
    atomPCTenrichment = calcAtomPCT(iso16O.weighted_mean_density, isotope_shift)
    samplePairs[iPair] = {'isotope_shift': isotope_shift,
                          'atomPCTenrichment': atomPCTenrichment}
    incubation = 7
    if row == 1:
        incubation = 30
    plotGraph(iso18O, 
              axs[row,col], 
              f'{iso16O.name[:6]} days: {incubation}  at% enrichment: {atomPCTenrichment:.2f}%',
              sample16=iso16O)

for ax in axs.flat:
    ax.set_xlim(1.6, 1.8)
    ax.set(xlabel='density (g/ml)', ylabel='DNA (ng/ul)')

for ax in axs.flat:
    ax.label_outer()

plt.plot()
plt.savefig('figures/batch129.svg', dpi=300)