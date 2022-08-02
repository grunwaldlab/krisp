from diagvar.find_diag_var import *
from diagvar.find_diag_var import _read_group_data, _count_genotypes, _group_genotypes, _diagnostic_gts, _check_variant, GroupedVariant
import pysam
from itertools import islice
import unittest


class TestSingleVariantInternals(unittest.TestCase):
    
    def setUp(self):
        pysam.set_verbosity(0)
        self.vcf = pysam.VariantFile('test_data/unfilt_allscafs_n666.vcf.gz')
        self.variant1 = next(self.vcf)
        self.variant2 = next(self.vcf)
        self.variant3 = next(self.vcf)
        self.diag_var = next(islice(self.vcf, 1020, 1021))
        self.groups = _read_group_data('test_data/test_metadata.tsv')
    
    def test_count_genotypes(self):
        self.assertEqual(_count_genotypes(self.variant1,
                                         ["9D1", "BS2014-584", "BS96"]),
                         {'T': 3})
        self.assertEqual(_count_genotypes(self.variant1),
                         {'T': 539, '?': 119, 'C': 8})
        self.assertEqual(_count_genotypes(self.variant2),
                         {'C': 558, '?': 102, 'C/G': 5, 'G': 1})
        self.assertEqual(_count_genotypes(self.variant2, hetero = False),
                         {'C': 563, '?': 102, 'G': 6})
    
    def test_group_genotypes(self):
        self.assertEqual(_group_genotypes(self.variant1, self.groups)["NA1"],
                         {'T': 316, '?': 13})
                         
    def test_diagnostic_gts(self):
        group_subset = {k: self.groups[k] for k in ('EU1', 'NA1')}
        self.assertEqual(_diagnostic_gts(self.diag_var, group_subset),
                         {'EU1': 'C', 'NA1': 'T'})
        self.assertEqual(_diagnostic_gts(self.variant1, group_subset),
                         {'EU1': None, 'NA1': None})
    
    def test_check_variant(self):
        group_subset = {k: self.groups[k] for k in ('EU1', 'NA1')}
        for var in self.vcf:
            print(_check_variant(var, group_subset))       
    
        

if __name__ == '__main__':
    unittest.main()
