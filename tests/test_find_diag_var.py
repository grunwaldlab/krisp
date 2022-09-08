from diagvar.find_diag_var import *
from diagvar.find_diag_var import _parse_group_data, _check_variant, GroupedVariant
import pysam
from itertools import islice
import unittest


class TestGroupedVariant(unittest.TestCase):
    
    def setUp(self):
        pysam.set_verbosity(0)
        self.vcf = pysam.VariantFile('test_data/unfilt_allscafs_n666.vcf.gz')
        self.groups = _parse_group_data('test_data/test_metadata.tsv')
        self.groups = {k: v for k, v in self.groups.items() if k in ['NA1','NA2', 'EU1', 'EU2']}
        self.variant1 = GroupedVariant(next(self.vcf), self.groups)
        self.variant2 = GroupedVariant(next(self.vcf), self.groups)
        self.variant3 = GroupedVariant(next(self.vcf), self.groups)
        self.diag_var = GroupedVariant(next(islice(self.vcf, 1020, 1021)), self.groups)
        self.vcf_subset = islice(self.vcf, 10, 15)

    def test_init(self):
        self.assertEqual(self.diag_var.sample_counts,
                         {'NA1': 222, 'EU1': 202, 'EU2': 0, 'NA2': 18})
        self.assertEqual(self.diag_var.allele_counts,
                         {'NA1': {'T': 222}, 'EU1': {'C': 202}, 'EU2': {}, 'NA2': {'C': 18}})
        self.assertEqual(self.diag_var.conserved,
                         {'NA1': 'T', 'EU1': 'C', 'EU2': None, 'NA2': 'C'})
        self.assertEqual(self.diag_var.diagnostic,
                         {'NA1': None, 'EU1': None, 'EU2': None, 'NA2': None})

    def test_from_vcf(self):
        gen = GroupedVariant.from_vcf(self.vcf_subset, groups=self.groups)
        self.assertEqual(next(gen).allele_counts,
                         {'NA1': {'T': 236}, 'EU1': {'C': 207}, 'EU2': {}, 'NA2': {'C': 22}})


if __name__ == '__main__':
    unittest.main()
