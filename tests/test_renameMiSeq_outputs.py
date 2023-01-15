import unittest
from scripts.rename_MiSeq_outputs import WaterYearSample

class TestWaterYearSample(unittest.TestCase):
    def testPost(self):
        s1 = WaterYearSample('A', '1', '1', (0, 10), 0.1, 0.1)
        self.assertEqual(s1.sampling_site, 'A')
        self.assertEqual(s1.replicate, '1')
        self.assertEqual(s1.sampling_week, '1')
        self.assertEqual(s1.depth, (0, 10))
        self.assertEqual(s1.sample_weight, 0.1)
        self.assertEqual(s1.initial_GWC, 0.1)
