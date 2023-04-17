import unittest
from scripts.parseSIPoutput import extract_fractionation_samples_from_excel_files, build_sample_objects, findIsotopePairs
from scripts.parseSIPoutput import calcAtomPCT, lindexsplit, splitsum, splitFractions, chunkProperties, improve_SIP_bins
from scripts.parseSIPoutput import SAMPLE_METADATA

class TestParseSIPoutput(unittest.TestCase):
    def testExtract_fractionation_samples_from_excel_files(self):
        # Test that the function can read the excel files
        fractionation_samples = extract_fractionation_samples_from_excel_files('./tests/test_data/read_files/')
        self.assertEqual(len(fractionation_samples), 22)
        self.assertEqual(fractionation_samples[2448][0], 'ABCD')

    def testDuplicateSamples(self):
        # Test that the function raises an error when duplicate samples are found
        with self.assertRaises(ValueError) as context:
            extract_fractionation_samples_from_excel_files('./tests/test_data/duplicate_samples/')
        self.assertTrue('Duplicate samples found!' in str(context.exception))

    def testBuild_sample_objects(self):
        # Test that the function can build sample objects
        fractionation_samples = extract_fractionation_samples_from_excel_files('./tests/test_data/read_files/')
        sample_objects = build_sample_objects(SAMPLE_METADATA ,fractionation_samples)
        self.assertEqual(len(sample_objects), 22)
        self.assertEqual(str(sample_objects[0]), 'AA3_D5_18O-7')
        self.assertEqual(sample_objects[0].plate, 'IJKL')
        self.assertEqual(sample_objects[0].id, 1429)
        self.assertEqual(sample_objects[0].experiment, '18O-7')
        self.assertEqual(sample_objects[0].core_id, 'AA3_D5')
        self.assertEqual(sample_objects[0].core_pH, 4.66)
        self.assertEqual(sample_objects[0].batch, '139')
        self.assertEqual(sample_objects[0].depth, ('50', '80'))
        self.assertEqual(sample_objects[0].depth_str, 'D5')

    def testFindIsotopePairs(self):
        # Test that the function can find isotope pairs
        fractionation_samples = extract_fractionation_samples_from_excel_files('./tests/test_data/paired_samples/')
        sample_objects = build_sample_objects(SAMPLE_METADATA ,fractionation_samples)
        isotope_pairs = findIsotopePairs(sample_objects)
        self.assertEqual(len(isotope_pairs), 3)
        self.assertEqual(str(isotope_pairs[0][0]), 'AA4_D5_16O-7')
        self.assertEqual(str(isotope_pairs[0][1]), 'AA4_D5_18O-7')



if __name__ == '__main__':
    unittest.main()