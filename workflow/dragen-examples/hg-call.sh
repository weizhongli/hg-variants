#!/bin/bash

my_host=`hostname`
my_pid=$$
my_core=18
my_queue=sh1
my_time_start=`date +%s`

cd PATH/hg-workflow/NA12878_XXX
mkdir hg-call
if ! [ -f PATH/hg-workflow/NA12878_XXX/hg-call/WF.start.date ]; then date +%s > PATH/hg-workflow/NA12878_XXX/hg-call/WF.start.date;  fi


dragen -f -r /ephemeral/data/refs/GRch38cnv \
    --fastq-list pull-fastq/fastq.csv \
    --fastq-list-sample-id NA12878_XXX \
    --lic-server=XXXYY_KEY \
    --enable-map-align true \
    --enable-map-align-output true \
    --enable-variant-caller true \
    --vc-emit-ref-confidence GVCF \
    --enable-sv true \
    --enable-cnv true \
    --cnv-enable-self-normalization true \
    --output-directory hg-call \
    --output-file-prefix HG > hg-call/WF.log

if [[ -s hg-call/HG.cnv.vcf.gz && -s hg-call/HG.sv.vcf.gz && -s hg-call/HG.hard-filtered.gvcf.gz ]]; then 
  #### remove fastq files
  rm -rf /ephemeral/fastq/NA12878_XXX
  rm -f  hg-call/HG.bam*
  rm -rf hg-call/sv/workspace/genomeDepth
  rm -rf hg-call/sort_spill
  rm -rf hg-call/vc
fi



if ! [ -s hg-call/HG.cnv.vcf.gz ]; then echo "zero size hg-call/HG.cnv.vcf.gz"; exit; fi
if ! [ -s hg-call/HG.sv.vcf.gz ]; then echo "zero size hg-call/HG.sv.vcf.gz"; exit; fi
if ! [ -s hg-call/HG.hard-filtered.gvcf.gz ]; then echo "zero size hg-call/HG.hard-filtered.gvcf.gz"; exit; fi

date +%s > PATH/hg-workflow/NA12878_XXX/hg-call/WF.complete.date
my_time_end=`date +%s`;
my_time_spent=$((my_time_end-my_time_start))
echo "sample=NA12878_XXX job=hg-call host=$my_host pid=$my_pid queue=$my_queue cores=$my_core time_start=$my_time_start time_end=$my_time_end time_spent=$my_time_spent" >> PATH/hg-workflow/NA12878_XXX/hg-call/WF.cpu

