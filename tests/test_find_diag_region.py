# from diagvar.find_diag_var import *
# from diagvar.find_diag_var import _read_group_data, _check_variant, GroupedVariant
# from diagvar.find_diag_region import GroupedRegion, find_diag_region
# import pysam
# from itertools import islice
# import unittest
# from Bio import SeqIO
#
# class TestGroupedRegion(unittest.TestCase):
#
#     def setUp(self):
#         pysam.set_verbosity(0)
#         self.vcf = pysam.VariantFile('test_data/unfilt_allscafs_n666.vcf.gz')
#         self.groups = _read_group_data('test_data/test_metadata.tsv')
#         self.ref = GroupedRegion._get_reference('test_data/PR-102_v3.1.fasta')
#         self.variant1 = next(self.vcf)
#         self.variant2 = next(self.vcf)
#         self.variant3 = next(self.vcf)
#         self.diag_var = next(islice(self.vcf, 1020, 1021))
#
#     def test_subset_sequence(self):
#         vars = [self.variant1, self.variant2, self.variant3]
#         self.assertEqual(GroupedRegion._subset_sequence(vars, self.groups['NA1']),
#                          ['<T316,?13>', '<C317,?12>', '<A317,?12>'])
#         self.assertEqual(GroupedRegion._subset_sequence(vars, self.groups['NA1'],
#                                                         min_geno_qual = 40,
#                                                         counts = False),
#                          ['<T>', '<C>', '<A>'])
#         self.assertEqual(GroupedRegion._subset_sequence(vars, self.groups['NA1'],
#                                                         min_geno_qual = 40,
#                                                         conserved = True),
#                          ['<T312>', '<C302>', '<A303>'])
#         self.assertEqual(GroupedRegion._subset_sequence(vars, self.groups['NA1'],
#                                                         conserved = True),
#                          [None, None, None])
#
#
# class TestFindingDiagnosticRegions(unittest.TestCase):
#
#     def setUp(self):
#         pysam.set_verbosity(0)
#         self.vcf = pysam.VariantFile('test_data/unfilt_allscafs_n666.vcf.gz')
#         self.groups = _read_group_data('test_data/test_metadata.tsv')
#         self.ref = GroupedRegion._get_reference('test_data/PR-102_v3.1.fasta')
#
#     def test_find_diag_region(self):
#         find_diag_region(self.vcf, self.groups)
