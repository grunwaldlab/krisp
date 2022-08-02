from diagvar.find_diag_var import *
import pysam

vcf = pysam.VariantFile('test_data/unfilt_allscafs_n666.vcf.gz')

# count_genotypes
print(count_genotypes(next(vcf), ["9D1", "BS2014-584", "BS96"]))
