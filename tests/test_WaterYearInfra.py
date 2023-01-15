import unittest
from scripts.waterYearSamplesInfrastructure import WaterYearSample

class TestWaterYearSample(unittest.TestCase):
    def testPost(self):
        s1 = WaterYearSample('A', 'A', '1', (0, 10), 'D1', '2020-01-01', './data/metadata/all_samples.csv')
        self.assertEqual(s1.sampling_site, 'A')
        self.assertEqual(s1.replicate, 'A')
        self.assertEqual(s1.sampling_week, '1')
        self.assertEqual(s1.depth, (0, 10))
        self.assertEqual(s1.depth_str, 'D1')
        self.assertEqual(str(s1.date), '2020-01-01')
        self.assertEqual(s1.metadata_file_path, './data/metadata/all_samples.csv')
        self.assertEqual(s1.sample_id, 'AA1_D1')
        self.assertEqual(s1.oneword_id, 'AA1D1')
        self.assertEqual(s1.core_GWC, 24.7396979)
        self.assertEqual(s1.core_pH, 4.63)

        # Test error producing inputs
        with self.assertRaises(FileNotFoundError) as context:
            WaterYearSample('A', 'A', '1', (0, 10), 'D1', 
                '2020-01-01', './data/metadata/all_sampless.csv')
        self.assertTrue("Metadata file ./data/metadata/all_sampless.csv not found!" in str(context.exception))
        with self.assertRaises(TypeError) as context:
            WaterYearSample('A', 'A', '1', (0, 10), 'D1', 
                '2020-01-01', './data/metadata/all_samples.csv', core_GWC=0.1)
        self.assertTrue("__init__() got an unexpected keyword argument 'core_GWC'" in str(context.exception))
        with self.assertRaises(TypeError) as context:
            WaterYearSample('A', 'A', '1', (0, 10), 'D1', 
                '2020-01-01', './data/metadata/all_samples.csv', core_pH=0.1)
        self.assertTrue("__init__() got an unexpected keyword argument 'core_pH'" in str(context.exception))
        with self.assertRaises(ValueError) as context:
            WaterYearSample('H', 'B', '1', (20, 30), 'D3',
                '2020-01-01', './data/metadata/all_samples.csv')
        self.assertTrue('Property GWC % is empty in metadata file at row HB1_D3_GSC!' in str(context.exception))