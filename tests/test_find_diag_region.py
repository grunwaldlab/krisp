from krisp_vcf.find_diag_var import *
from krisp_vcf.find_diag_var import _parse_group_data, _check_variant, GroupedVariant
from krisp_vcf.find_diag_region import GroupedRegion, find_diag_region
import pysam
from itertools import islice
import unittest
from Bio import SeqIO


class TestGroupedRegion(unittest.TestCase):
    def setUp(self):
        pysam.set_verbosity(0)
        self.vcf = pysam.VariantFile('test_data/unfilt_allscafs_n666.vcf.gz')
        self.groups = _parse_group_data('test_data/test_metadata.tsv')
        self.groups = {g: v for g, v in self.groups.items() if g in ["NA1", "NA2", "EU1", "EU2"]}
        self.ref = 'test_data/PR-102_v3.1.fasta'
        self.variant1 = next(self.vcf)
        self.variant2 = next(self.vcf)
        self.variant3 = next(self.vcf)
        self.diag_var = next(islice(self.vcf, 1020, 1021))
        self.vcf_subset = islice(self.vcf, 10, 15)

    def test_init(self):
        x = GroupedRegion(GroupedVariant.from_vcf(self.vcf_subset, groups=self.groups), group="NA1", reference=self.ref)
        print(x.region_length())

    # def test_sliding_window(self):
    #     windower = GroupedRegion.sliding_window(GroupedVariant.from_vcf(self.vcf, groups=self.groups),
    #                                             groups=self.groups.keys(),
    #                                             reference=self.ref,
    #                                             span=30)
    #     for x in windower:
    #         print(x.group, x.region_length())
    #         print(len(x.variants))

    def test_find_diag_region(self):
        find_diag_region(self.vcf,
                         groups=self.groups,
                         reference=self.ref)

