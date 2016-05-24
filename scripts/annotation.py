#!/usr/bin/python

import sys
import os
import pandas as pd
import numpy as np
import argparse
import vcf as pyvcf

# define data paths

bed_path = "source_data/roadmap_r9/15_state_model/bed/"

def get_options():

        parser = argparse.ArgumentParser()

        parser.add_argument('--input', action='store', dest='input_file', required=True,
                    help='Input vcf file')
        parser.add_argument('--output', action='store', dest='output_file', required=True,
                    help='Output file')
        parser.add_argument('--eid', nargs="+", action='store', dest='eid_name', required=False, default=['E116'],
                    help='Epigenome name for annotation (default = GM12878_Lymphoblastoid')

        args =  parser.parse_args()
        return args

def write_config_file(eid_list):
    """
    Write a VEP configuration file.
    """

    config_arguments = ("biotype            1"      "\n"
                        "cache              1"      "\n"
                        "canonical          1"      "\n"
                        "assembly           GRCh38" "\n"
                        "core_type          core"   "\n"
                        "dir                source_data/ensembl/cache" "\n"
                        "dir_cache          source_data/ensebl/cache" "\n"
                        "host               130.88.97.228"  "\n"
                        "port               3306"   "\n"          
                        "force_overwrite    1"      "\n"
                        "numbers            1"      "\n"
                        "polyphen           p"      "\n"
                        "sift               p"      "\n"
                        "symbol             1"      "\n"
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

    vep_cmd = "variant_effect_predictor.pl -i %s -o %s --config temp.config" % (input, vep_vcf)
    os.system(vep_cmd)
    os.system("rm temp.config")

def main():

    # define input/output
    options = get_options()

    write_config_file(options.eid_name)
    run_vep(options.input_file, options.output_file)


    # format vcf into table

    header = ['ID', 'CHROM', 'POS', 'REF', 'ALT']

    df = pd.DataFrame(columns=header)

    vcf = pyvcf.Reader(open('output/annotation/annotation.vcf', 'r'))

    csq_subset = [1,2,3,5,7,8,9,23,24]

    for record in vcf:
        result = pd.Series([record.ID, record.CHROM, record.POS, record.REF, record.ALT], index=header)
        df = df.append(result, ignore_index=True)
        for csq in record.INFO['CSQ']:
            #print record.ID, csq
            details = csq.split("|")
            #print [details[i] for i in csq_subset]

    # write table
    table_file = options.output_file + ".csv"
    df.to_csv(table_file, index=False)


if __name__=='__main__':
    main()