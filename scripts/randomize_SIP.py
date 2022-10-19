#!/usr/bin/env python3

from tokenize import group
import numpy as np
from random import shuffle, seed
import csv, re, sys, string

def check_replicates(groupedEntries:list):
    '''Checks a list for repeated replicate samples.
    For example HB4_D2... and HA4_D2...'''
    removed_isotope_incubationLength = [x[:6] for x in groupedEntries]
    removed_replicate_identifier = [x[0]+x[2:] for x in removed_isotope_incubationLength]
    return len(set(removed_replicate_identifier)) < len(groupedEntries)

def check_weights(groupedEntries:list):
    '''Checks a list for even number of D1-2 samples.'''
    d1_or_d2, d3 = list(), list()
    for x in groupedEntries:
        if '_D1_' in x or '_D2_' in x:
            d1_or_d2.append(x)
        if '_D3_' in x:
            d3.append(x)
    return len(d1_or_d2) % 2 != 0 or len(d3) % 2 != 0

def find_and_remove_failed_shuffles(list_of_groupedEntries:list):
    '''Finds groups of entries where the repeated replicate rule fails,
    then removes these entries from the original list and returns them as 
    separate list in addition to the succeeded list'''
    removal_idx, failed_shuffles = list(), list()
    for i, x in enumerate(list_of_groupedEntries):
        if check_replicates(x) or check_weights(x):
            removal_idx.append(i)
            failed_shuffles.append(x)

    list_of_groupedEntries_withoutFails = [ele for idx, ele in enumerate(list_of_groupedEntries) if idx not in removal_idx]
    return list_of_groupedEntries_withoutFails, failed_shuffles

def reshuffle_failed_shuffle(list_of_groupedEntries:list):
    '''Flattens a list of lists which contain repeated replicate samples.
    Then reshuffles this list and groups them by n again.'''
    n = len(list_of_groupedEntries[0])
    flattened = [x for sublist in list_of_groupedEntries for x in sublist]
    shuffle(flattened)
    return [flattened[i:i + n] for i in range(0, len(flattened), n)]

def recursively_shuffle_until_noFails(sip_entries:list, n=6):
    '''Shuffles and splits the sip_entries into list of list with kength n.
    Checks for failed shuffles following the repeated replicate rule in each sublist.
    Repeats same shuffle with failed sublists until no fails occur or only one sublist keeps failing.'''

    split_sip_entries = [sip_entries[i:i + n] for i in range(0, len(sip_entries), n)]
    shuffled_sip_entries = reshuffle_failed_shuffle(split_sip_entries)
    success, fail = find_and_remove_failed_shuffles(shuffled_sip_entries)
    if len(fail) == 1:
        sys.exit("Failed to shuffle, try different seed.")
    elif len(fail) == 0:
        return success
    elif len(fail) > 1:
        success.extend(recursively_shuffle_until_noFails([x for sublist in fail for x in sublist], n))
        return success

sip_entries18, sip_entries16 = list(), list()
sip_entries_withDB18, sip_entries_withDB16 = dict(), dict()
with open('./data/all_samples_afterW9.csv') as csvFile:
    reader = csv.reader(csvFile)
    #seed(int(sys.argv[1]))
    seed(2)
    for row in reader:
        #Skip lines starting with #
        if '#' == row[0][0]:
            continue
        #Site to exclude
        #if 'H' == row[1][0]:
        #    continue
        #Weeks to exclude
        if '6_D' in row[1] or '1_D' in row[1]:
            continue
        if 'HA7_D' in row[1] or 'HB7_D' in row[1] or 'HC7_D' in row[1]:
            continue 
        #Depths to exclude
        if '_D4_' in row[1] or '_D2_' in row[1]:
            continue
        #If we are doing only D1 and D2
        # if '_D4_' in row[1]:
        #     continue
        # if '_D5_' in row[1]:
        #     continue
        #Depending on initial SIP results:
        if '-30' in row[1]:
           #seed(26)
           continue
        # if re.match(r'.{3}_D[4|5]_.{3}-7', row[1]):
        #     seed(24)
        #     continue
        
        if '18O-' in row[1]:
            sip_entries_withDB18[row[1]] = (row[0].split('/')[4], row[9])
            sip_entries18.append(row[1])
        elif '16O-' in row[1]:
            sip_entries_withDB16[row[1]] = (row[0].split('/')[4], row[9])
            sip_entries16.append(row[1])

n = 6
shuffled_16O = recursively_shuffle_until_noFails(sip_entries16, n)
shuffled_18O = recursively_shuffle_until_noFails(sip_entries18, n)

groupIDs = list(string.ascii_letters)

print("16O groups:")
for i,x in enumerate(shuffled_16O):
    print(f"Group {groupIDs[i]}")
    for oneTosix,y in enumerate(x):
        print(f"{groupIDs[i]}{oneTosix+1} {sip_entries_withDB16[y][0]}\t{y}\t{sip_entries_withDB16[y][1]}")
    #print(','.join(x))
    #print(','.join([sip_entries_withDB16[y] for y in x]))

print()
print("18O groups:")
for k,x in enumerate(shuffled_18O):
    print(f"Group {groupIDs[i+k+1]}")
    for oneTosix,y in enumerate(x):
        print(f"{groupIDs[i+k+1]}{oneTosix+1} {sip_entries_withDB18[y][0]}\t{y}\t{sip_entries_withDB18[y][1]}")
    #print(','.join(x))
    #print(','.join([sip_entries_withDB18[y] for y in x]))

print(f'Total 16O extractions: {len(sip_entries16)}. Total 16O groups: {len(shuffled_16O)}')
print(f'Total 18O extractions: {len(sip_entries18)}. Total 18O groups: {len(shuffled_18O)}')