#!/usr/bin/python

import sys
import os
import pandas as pd
import numpy as np
import argparse

def get_options():

	parser = argparse.ArgumentParser()

	parser.add_argument('--state', action='store', dest='state',
			required=True, help='chromatin state model (currently 15 or 18)')
	parser.add_argument('--meta', action='store', dest='meta_file',
			required=True, help='metadata file linking eid with cell type')
	parser.add_argument('--rename', action='store', dest='rename',
			required=True, help='label to rename output file (EID, Mnemonic, EDACC')

	args =  parser.parse_args()
	return args

def main():
	# define input/output
	options = get_options()

	# read metadata file
	metadata = pd.read_table(options.meta_file, sep=",")
	metadata.rename(columns={'Epigenome Mnemonic':'Mnemonic'}, inplace=True)
	metadata.rename(columns={'Epigenome name (from EDACC Release 9 directory)':'EDACC'}, inplace=True)

	# process each file
	for i,o in zip(metadata.EID, metadata[options.rename]):

		# define files
		raw_file = "source_data/roadmap_r9/%s_state_model/raw/%s_%s_coreMarks_dense.bed.gz" % (options.state, i, options.state)
		bed_file = "source_data/roadmap_r9/%s_state_model/bed/%s.bed.gz" % (options.state, o)

		# sort by chromosome and bp
		cmd = "gunzip -c %s | sed '1 d' | sort -k1,1 -k2,2n -k3,3n | bgzip > %s" % (raw_file, bed_file)
		os.system(cmd)

		#index
		cmd = "tabix -p bed %s" % (bed_file)
		os.system(cmd)

if __name__=='__main__':
        main()
