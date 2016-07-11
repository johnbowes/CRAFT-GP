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

    # config_arguments = ("biotype            1"      "\n"
    #                     "cache              1"      "\n"
    #                     "canonical          1"      "\n"
    #                     "assembly           GRCh38" "\n"
    #                     "core_type          core"   "\n"
    #                     "dir                source_data/ensembl/cache" "\n"
    #                     "dir_cache          source_data/ensebl/cache" "\n"
    #                     "host               130.88.97.228"  "\n"
    #                     "port               3306"   "\n"          
    #                     "force_overwrite    1"      "\n"
    #                     "numbers            1"      "\n"
    #                     "polyphen           p"      "\n"
    #                     "sift               p"      "\n"
    #                     "symbol             1"      "\n"
    #                     "vcf                1"      "\n")

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
    

    df.replace("", ".", inplace=True)

    return df

def main():

    # define input/output
    options = get_options()

    # run VEP annotation
    write_config_file(options.eid_name)
    run_vep(options.input_file, options.output_file)

    # tabulate vcf and write to file
    table = tabulate_vcf(options.output_file + '.vcf', options.eid_name)
    table_file = options.output_file + ".csv"
    table.to_csv(table_file, index=False)


if __name__=='__main__':
    main()