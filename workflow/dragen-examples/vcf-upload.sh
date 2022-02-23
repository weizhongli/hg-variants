#!/bin/bash

my_host=`hostname`
my_pid=$$
my_core=4
my_queue=sh1
my_time_start=`date +%s`

cd PATH/hg-workflow/NA12878_XXX
mkdir vcf-upload
if ! [ -f PATH/hg-workflow/NA12878_XXX/vcf-upload/WF.start.date ]; then date +%s > PATH/hg-workflow/NA12878_XXX/vcf-upload/WF.start.date;  fi


for i in HG.hard-filtered.gvcf.gz HG.cnv.vcf.gz HG.sv.vcf.gz; do
  echo "s3 cp hg-call/$i aws_s3_path_for_vcf_and_results_upload_e.g.s3://bucket/folder_level1/folder_level2/NA12878_XXX/$i"
  echo "s3 cp hg-call/$i.tbi aws_s3_path_for_vcf_and_results_upload_e.g.s3://bucket/folder_level1/folder_level2/NA12878_XXX/$i.tbi"
  aws s3 cp hg-call/$i aws_s3_path_for_vcf_and_results_upload_e.g.s3://bucket/folder_level1/folder_level2/NA12878_XXX/$i
  aws s3 cp hg-call/$i.tbi aws_s3_path_for_vcf_and_results_upload_e.g.s3://bucket/folder_level1/folder_level2/NA12878_XXX/$i.tbi
done 

for i in HG-replay.json WF.log; do
  echo "s3 cp hg-call/$i aws_s3_path_for_vcf_and_results_upload_e.g.s3://bucket/folder_level1/folder_level2/NA12878_XXX/$i"
  aws s3 cp hg-call/$i aws_s3_path_for_vcf_and_results_upload_e.g.s3://bucket/folder_level1/folder_level2/NA12878_XXX/$i
done 

tar cvf work-files.tar.gz --exclude=hg-call/*vcf.gz* --gzip pull-fastq hg-call
aws s3 cp work-files.tar.gz aws_s3_path_for_vcf_and_results_upload_e.g.s3://bucket/folder_level1/folder_level2/NA12878_XXX/work-files.tar.gz

aws s3 ls aws_s3_path_for_vcf_and_results_upload_e.g.s3://bucket/folder_level1/folder_level2/NA12878_XXX/ > vcf-upload/upload.list

rm -f hg-call/HG.hard-filtered.gvcf.gz


if ! [ -s vcf-upload/upload.list ]; then echo "zero size vcf-upload/upload.list"; exit; fi

date +%s > PATH/hg-workflow/NA12878_XXX/vcf-upload/WF.complete.date
my_time_end=`date +%s`;
my_time_spent=$((my_time_end-my_time_start))
echo "sample=NA12878_XXX job=vcf-upload host=$my_host pid=$my_pid queue=$my_queue cores=$my_core time_start=$my_time_start time_end=$my_time_end time_spent=$my_time_spent" >> PATH/hg-workflow/NA12878_XXX/vcf-upload/WF.cpu

