import unittest
from scripts.parseSIPoutput import extract_fractionation_samples_from_excel_files, build_sample_objects, findIsotopePairs, calcAtomPCT, lindexsplit, splitsum, splitFractions, chunkProperties, improve_SIP_bins

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


if __name__ == '__main__':
    unittest.main()