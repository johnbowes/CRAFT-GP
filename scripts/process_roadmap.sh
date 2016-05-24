#!/bin/bash

# Simple script to sort and index roadmap data for use in CRAFT
# annotation. Input is model state and epigenome id (eid).
#
# Requires tabix: module load apps/gcc/tabix/0.2.6
#
# Example usage for 15 state model and CD4_Memory_Primary_Cells:
#
# sh scripts process_roadmap.sh 15 E037 

model=$1
eid=$2

input="source_data/roadmap_r9/${model}_state_model/raw/${eid}_${model}_coreMarks_dense.bed.gz"
output="source_data/roadmap_r9/${model}_state_model/bed/${eid}.bed.gz"

gunzip -c $input | sed '1 d' | sort -k1,1 -k2,2n -k3,3n | bgzip > $output
tabix -p bed $output