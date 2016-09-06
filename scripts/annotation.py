#!/usr/bin/python

import sys
import os
import re
import pandas as pd
import numpy as np
import argparse
import vcf as pyvcf

# define data paths

bed_path = "source_data/roadmap_r9/15_state_model/bed/"
meta_file = "source_data/roadmap_r9/meta_data/roadmap_consolidated_epigenome_ids.txt"

def get_options():

        parser = argparse.ArgumentParser()

        parser.add_argument('--input', action='store', dest='input_file', required=True,
                    help='Inputcredible snp file')
        parser.add_argument('--output', action='store', dest='output_file', required=True,
                    help='Output file')
        parser.add_argument('--epi_names', action='store', dest='epi_names', required=False, default='ENCODE2012',
                    help='Group name or list of epigenome names (seperated by space or comma')
        parser.add_argument('--epi_type', action='store', dest='epi_type', required=False, default='group',
                    help='How are epigenome names specified: list or group')

        args =  parser.parse_args()
        return args

def get_group_ids(group_list, meta_file):
    """
    Return a list of epigenome IDs based on a group label.
    """

    group_list = group_list.split(",")

    meta = pd.read_table(meta_file, sep="\t")
    meta = meta[meta['Group'].isin(group_list)]

    meta_list = meta['Epigenome name (from EDACC Release 9 directory)'].tolist()

    return meta_list

def write_config_file(eid_list):
    """
    Write a VEP configuration file.
    """

    config_arguments = ("biotype            1"      "\n"
                        "cache              1"      "\n"
                        "canonical          1"      "\n"
                        "assembly           GRCh37" "\n"
                        "cache_version      84"     "\n"
                        "core_type          core"   "\n"
                        "dir                source_data/ensembl/cache" "\n"
                        "dir_cache          source_data/ensebl/cache" "\n"
                        "host               130.88.97.228"  "\n"
                        "port               3337"   "\n"
                        "plugin             CADD,source_data/ensembl/plugins/CADD/1000G.tsv.gz" "\n"
                        "force_overwrite    1"      "\n"
                        "numbers            1"      "\n"
                        "polyphen           p"      "\n"
                        "sift               p"      "\n"
                        "symbol             1"      "\n"
                        "chr                1-22"   "\n"
                        "vcf                1"      "\n")

    custom_arguments = "custom\t" 

    for id in eid_list:
        bed_file = bed_path + id + ".bed.gz"
        custom_arguments = custom_arguments + "\t" + bed_file + "," + id + "," + "bed,overlap,0"

    config_arguments = config_arguments + custom_arguments

    with open("temp.config", "w") as config_file:
        config_file.write(config_arguments)


def run_vep(input, output):
    """
    Run VEP

    Currently use incline module.
    """

    vep_vcf = output + ".vcf"   

    vep_cmd = "variant_effect_predictor.pl -i %s -o %s --config temp.config -v" % (input, vep_vcf)
    os.system(vep_cmd)
    os.system("rm temp.config")

def tabulate_vcf(vcf_file, eid_list):
    """
    Tabulate vcf file in csv file.
    """
    vcf = pyvcf.Reader(open(vcf_file, 'r'))

    # define results table df
    info_fields = re.search('Format: (.+)', vcf.infos['CSQ'].desc).group(1).split("|")
    header = ['ID', 'CHROM', 'POS', 'REF', 'ALT'] + info_fields + eid_list
    df = pd.DataFrame(columns=header)
    
    for record in vcf:
        snp_info = [record.ID, record.CHROM, record.POS, record.REF, record.ALT]
        eid_info = []
        for eid in eid_list:
            eid_info.append(record.INFO[eid])
        for csq in record.INFO['CSQ']:
            csq_info = csq.split("|")
            test = snp_info + csq_info + eid_info
            result = pd.Series(test, index=header)
            df = df.append(result, ignore_index=True)

    # remove unwanted fields
    fields_to_remove = ['FLAGS','SYMBOL_SOURCE','HGNC_ID','CANONICAL','HGVSc','HGVSp','cDNA_position','CDS_position','Protein_position','Amino_acids','Codons','Existing_variation']
    df.drop(fields_to_remove, axis=1, inplace=True)

    df.replace("", ".", inplace=True)

    return df

def main():

    # define input/output
    options = get_options()

    # process input file
    credible_snps = pd.read_table(options.input_file, sep = " ")
    credible_snps = credible_snps[['SNPID', 'PVAL', 'pp', 'cpp']]

    # create a SNP list for VEP input
    snp_list_file = options.output_file + "snp_list.txt"
    credible_snps.to_csv(snp_list_file, columns = ['SNPID'], index = False, header = False)

    # get list of epigenomes
    if options.epi_type == "list":
        epigenomes = options.epi_names.split(",")
    else:
        epigenomes = get_group_ids(options.epi_names, meta_file)

    # write epigenomes to file
    epigenome_file = options.output_file + '.epigenomes'
    epi_string = ','.join(epigenomes)
    with open(epigenome_file, "w") as file:
        file.write(epi_string)

    # run VEP annotation
    write_config_file(epigenomes)
    run_vep(snp_list_file, options.output_file)

    # tabulate vcf
    table = tabulate_vcf(options.output_file + '.vcf', epigenomes)

    # merge posterior probabilities
    table = pd.merge(table, credible_snps, how = 'left', left_on = 'ID', right_on = 'SNPID')

    # write table to file
    table_file = options.output_file + ".csv"
    table.to_csv(table_file, index=False)

if __name__=='__main__':
    main()