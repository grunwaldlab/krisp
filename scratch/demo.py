



	
def _mask_same(seqs, reference):
	"""
	Parameters
	----------
	seqs : dict of list of str
		sequecnes to align, named by group
	reference : list of str
		The refenrece used to call variants	
	
	Returns
	-------
	list str
		Modified version of seqs
	"""
	for group in seqs.keys():
		for index in range(len(seqs[group])):
			if seqs[group][i] == reference[i]:
				seqs[group][i] = "."
	return(seqs)

def _pad_columns(seqs, reference):
	"""Adding space to account for formatting (non indels)."""
	return seqs, reference
	
def _pad_indel(seqs, reference):
	"""Add dashes to pad indels."""
	return seqs, reference

def _print_alignment(seqs, reference, ref_name = "Reference"):
	ref_name = ref_name.rjust(max_len, " ")
	print(f'{ref_name} : ' + ''.join(reference))
	for group_name, seq in seqs.items():
		group_name = group_name.rjust(max_len, " ")
		print(group_name + ' : ' + ''.join(seq))


def render_variant(seqs, reference):
	"""my short summary
	
	my llong summaryyyyyyyyyyyyyyyyyyyyyyyyyyyyy
	yyyyyyyyyyyyyyyyyyyyyyyy
	yyyyyyyyyyyy
	
	Parameters
	----------
	seqs : dict of list of str
		sequecnes to align, named by group
	reference : list of str
		The refenrece used to call variants
	"""
	
	seqs = _mask_same(seqs, reference)
	seqs, reference = _pad_columns(seqs, reference)
	seqs = _pad_indel(seqs, reference)
	
	_print_alignment(seqs, reference)
	
	


if __name__ == "__main__":
	
	render_variant()
