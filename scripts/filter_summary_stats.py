#!/usr/bin/python

import sys
import os
import pandas as pd
import numpy as np
import argparse

def get_options():

        parser = argparse.ArgumentParser()

        parser.add_argument('--regions', action='store', dest='regions', required=True,
                    help='File containing region information')
	parser.add_argument('--stats', action='store', dest='stats', required=True,
                    help='File containing summary statistics')
        parser.add_argument('--out', action='store', dest='out_file', required=True,
                    help='Output file')

        args =  parser.parse_args()
        return args

def main():

        # define input/output
        options = get_options()

        # read summary stats file
        stats_file = options.stats
        stats = pd.read_table(stats_file, sep="\t")

	# read region file
	regions_file = options.regions
	regions = pd.read_table(regions_file, header=None, names=['index_snp','locus'])

	regions['chr'] = regions['locus'].str.extract('chr(.*):').astype(int)
	regions['start'] = regions['locus'].str.extract(':(.*)-').astype(int)
	regions['end'] = regions['locus'].str.extract('-(.*)').astype(int)

	# iterate through regions
	subset_list = []

	for index, row in regions.iterrows():
		
		subset = stats[np.logical_and(stats['CHROM'] == row['chr'],
			np.logical_and(stats['POS'] >= row['start'], stats['POS'] <= row['end']))]
		subset.is_copy = False
		subset['index_snp'] = row['index_snp']
		subset_list.append(subset)

	# concatenate list of dataframes
	stats_subset = pd.concat(subset_list)

        # write exclusion file
        out_file = options.out_file
        stats_subset.to_csv(out_file, sep=' ', index=False)

if __name__=='__main__':
        main()
