def my_new_function(my_arg1, my_arg2, my_opt = TRUE):
    # All the code
    return None



ref =            ["A", "C", "GTC", "T",   "T", "G",  "G",         "C", "C"]
seqs = {"seq 1": ["A", "G", "GTC", "T/G", "T", "GC", "G",         "A", "C"],
        "seq 2": ["A", "G", "G",   "T",   "",  "G",  "<123G12?>", "A", "C"]}

my_new_function(ref, seqs, mask_same = TRUE)

# Should print:
#
# Reference : ACGTC T TG-    G    CC
#     seq 1 : .G...T/GTGC    G    A.
#     seq 2 : .G.-- T -G-<1234G2?>A. 
              
