import unittest
from scripts.parseSIPoutput import extract_fractionation_samples_from_excel_files, build_sample_objects, findIsotopePairs
from scripts.parseSIPoutput import calcAtomPCT, lindexsplit, splitsum, extract_heavy_samples, mergeFractions, improve_SIP_bins
from scripts.parseSIPoutput import load_picogreen_quantification_data
from scripts.parseSIPoutput import SAMPLE_METADATA

class TestParseSIPoutput(unittest.TestCase):

    #Add a fixture to load the data
    @classmethod
    def setUpClass(cls):
        cls.fractionation_samples = extract_fractionation_samples_from_excel_files('./tests/test_data/read_files/')
        cls.paired_samples = extract_fractionation_samples_from_excel_files('./tests/test_data/paired_samples/')
        cls.picogreen_data = load_picogreen_quantification_data('./tests/test_data/pico_data/')
        cls.sample_objects = build_sample_objects(SAMPLE_METADATA ,cls.fractionation_samples, cls.picogreen_data)
        cls.paired_sample_objects = build_sample_objects(SAMPLE_METADATA ,cls.paired_samples, cls.picogreen_data)

    def testExtract_fractionation_samples_from_excel_files(self):
        # Test that the function can read the excel files
        self.assertEqual(len(self.fractionation_samples), 22)
        self.assertEqual(self.fractionation_samples[2448][0], 'ABCD')

    def testDuplicateSamples(self):
        # Test that the function raises an error when duplicate samples are found
        with self.assertRaises(ValueError) as context:
            extract_fractionation_samples_from_excel_files('./tests/test_data/duplicate_samples/')
        self.assertTrue('Duplicate samples found!' in str(context.exception))

    def testPicoGreenLoad(self):
        # Test that the function can read the pico green data
        self.assertEqual(len(self.picogreen_data), 5)
        self.assertEqual(len(self.picogreen_data['122']['IJKL']), 3)
        self.assertEqual(self.picogreen_data['122']['IJKL'][0][0][15], 'A2')

    def testBuild_sample_objects(self):
        # Test that the function can build sample objects
        self.assertEqual(len(self.sample_objects), 22)
        self.assertEqual(str(self.sample_objects[0]), 'AA3_D5_18O-7')
        self.assertEqual(self.sample_objects[0].plate, 'IJKL')
        self.assertEqual(self.sample_objects[0].id, 1429)
        self.assertEqual(self.sample_objects[0].experiment, '18O-7')
        self.assertEqual(self.sample_objects[0].core_id, 'AA3_D5')
        self.assertEqual(self.sample_objects[0].core_pH, 4.66)
        self.assertEqual(self.sample_objects[0].batch, '139')
        self.assertEqual(self.sample_objects[0].depth, ('50', '80'))
        self.assertEqual(self.sample_objects[0].depth_str, 'D5')

    def testFindIsotopePairs(self):
        # Test that the function can find isotope pairs
        sample_objects = self.paired_sample_objects
        isotope_pairs = findIsotopePairs(sample_objects)
        self.assertEqual(len(isotope_pairs), 3)
        self.assertEqual(str(isotope_pairs[0][0]), 'AA4_D5_16O-7')
        self.assertEqual(str(isotope_pairs[0][1]), 'AA4_D5_18O-7')

    def testCalcAtomPCT(self):
        # Test that the function can calculate atom percent enrichment
        wm16 = 1.7
        isotope_shift = 0.01
        atom_pct = calcAtomPCT(wm16, isotope_shift)
        self.assertEqual(atom_pct, 14.966005979112884)

    def testLindexsplit(self):
        my_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        sublists = lindexsplit(my_list, 5, 7, 9)
        self.assertEqual(sublists, [[1, 2, 3, 4, 5, 6], [7, 8], [9, 10], []])
        noList = lindexsplit(my_list)
        self.assertEqual(noList, [])

    def testSplitsum(self):
        my_list = [1, 1, 3, 4, 5, 1, 0.1, 0.1]
        split_sum = splitsum(my_list, 2)
        self.assertEqual(split_sum, [[1, 1, 3],[4],[5],[1,0.1, 0.1]])

    def testExtractHeavySamples(self):
        paired_fractionation_samples = self.paired_sample_objects
        iso18Osamples = extract_heavy_samples(paired_fractionation_samples)
        self.assertEqual(len(iso18Osamples), 3)
        self.assertEqual(str(iso18Osamples[0]), 'AA4_D5_18O-7')
        self.assertEqual(str(iso18Osamples[1]), 'AB4_D5_18O-7')
        iso16Osamples = extract_heavy_samples(paired_fractionation_samples, '_16O-')
        self.assertEqual(len(iso16Osamples), 3)
        self.assertEqual(str(iso16Osamples[0]), 'AA4_D5_16O-7')
        self.assertEqual(str(iso16Osamples[1]), 'AB4_D5_16O-7')

    def testMergeFractions(self):
        paired_fractionation_samples = self.paired_sample_objects
        iso18Osamples = extract_heavy_samples(paired_fractionation_samples)
        self.assertFalse(hasattr(iso18Osamples[0], 'chunked_concentrations'))
        binned_samples = mergeFractions(iso18Osamples, targetDNA=80)
        self.assertEqual(len(binned_samples), 3)
        self.assertTrue(hasattr(binned_samples[0], 'chunked_concentrations'))
        self.assertEqual(len(binned_samples[0].chunked_concentrations), 8)
        self.assertEqual(len(binned_samples[0].chunked_concentrations[0]), 5)
        self.assertEqual(len(binned_samples[0].chunked_wells[0]), 5)
        binned_samples = mergeFractions(iso18Osamples, targetDNA=100)
        self.assertEqual(len(binned_samples[0].chunked_concentrations), 8)
        self.assertEqual(len(binned_samples[0].chunked_concentrations[0]), 6)
        self.assertEqual(len(binned_samples[0].chunked_wells[0]), 6)

    def testimprove_SIP_bins(self):
        paired_fractionation_samples = self.sample_objects
        iso18Osamples = extract_heavy_samples(paired_fractionation_samples)
        binned_samples = mergeFractions(iso18Osamples, targetDNA=80)
        self.assertEqual(len(binned_samples[4].chunked_densities), 8)
        self.assertEqual(len(binned_samples[4].chunked_densities[0]), 6)
        improve_SIP_bins(binned_samples)
        self.assertEqual(len(binned_samples[4].chunked_densities), 5)
        self.assertEqual(len(binned_samples[4].chunked_densities[0]), 7)

if __name__ == '__main__':
    unittest.main()