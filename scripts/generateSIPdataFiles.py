#!/usr/bin/env python3

from pandas import read_excel
from rename_MiSeq_outputs import parseSIPoutput, createSIPSampleClasses, SIPSample


def main():
    df = read_excel("./data/fractionation/220817_Batch_122_Water_Yr_Summary_mod.xlsx", sheet_name='Summary', engine='openpyxl')

    sipData = parseSIPoutput(df)
    sipSamples = createSIPSampleClasses(SIPSample, sipData, './data/metadata/all_samples_afterW10.csv')
    print(sipSamples)

if __name__ == "__main__":
    main()