#!/bin/bash

my_host=`hostname`
my_pid=$$
my_core=4
my_queue=sh0
my_time_start=`date +%s`

cd PATH/hg-workflow/NA12878_XXX
mkdir pull-fastq
if ! [ -f PATH/hg-workflow/NA12878_XXX/pull-fastq/WF.start.date ]; then date +%s > PATH/hg-workflow/NA12878_XXX/pull-fastq/WF.start.date;  fi


#### sleep up to 10 s
sleep $(( $RANDOM % 10  ))
#### wait unless there is less than 2 folderis in /ephemeral/fastq/
#### so that fastq files for all samples to be downloaded be <=3 
#### fastq folder will be deleted in the next job
if [ ! -e /ephemeral/fastq ]
then
  mkdir -p /ephemeral/fastq
fi

num=$(ls -1 /ephemeral/fastq | grep -c .)
while [ $num -gt 2 ]; do
  echo wait 
  sleep 60
  num=$(ls -1 /ephemeral/fastq | grep -c .)
done
mkdir -p /ephemeral/fastq/NA12878_XXX

aws s3 cp aws_s3_path_for_flowcell_info_tsv_file_e.g.s3://bucket/folder_level1/flowcell_all.tsv pull-fastq/sample-all.tsv
grep NA12878_XXX pull-fastq/sample-all.tsv | grep $DATA.0 > pull-fastq/sample.tsv
rm -f pull-fastq/sample-all.tsv

if [ ! -s pull-fastq/sample.tsv ]
then
  exit
fi

#### embeded python to select lanes up to 40GB
#### to be run as: env python job1.py NA12878_XXX  40
cat << ENDPY > job1.py
import re
import sys

input_sample_name = sys.argv[1]
target_depth = float(sys.argv[2])
f = open('pull-fastq/sample.tsv')
f_csv = open('pull-fastq/fastq.csv', 'w')
fout = open('pull-fastq/fastq.list', 'w')
f_csv.write('RGID,RGSM,RGLB,Lane,Read1File,Read2File\n')

depth_total = 0
for line in f:
  ll = re.split('\t', line.rstrip())
  depth_total = depth_total + float(ll[7])
  rgid = 'NA12878_XXX' + 'L' + str(ll[3])
  #### write fastq.csv
  f_csv.write('{0},{1},{2},{3},/ephemeral/fastq/NA12878_XXX/{4},/ephemeral/fastq/NA12878_XXX/{5}\n'.format(rgid, 'NA12878_XXX', 'LIB1', ll[3], ll[10], ll[12] ))

  fout.write('{0}\t{1}\n'.format(ll[10], ll[11]))
  fout.write('{0}\t{1}\n'.format(ll[12], ll[13]))
  if depth_total > target_depth: break

f.close()
f_csv.close()
fout.close()
ENDPY
env python job1.py NA12878_XXX 40 
sort pull-fastq/fastq.list > pull-fastq/fastq.list.1
mv -f pull-fastq/fastq.list.1 pull-fastq/fastq.list

#### download fastqs
for i in $(cut -f 1 pull-fastq/fastq.list | sort) 
  do
  echo "aws s3 cp aws_s3_path_for_input_fastq_e.g.s3://bucket/folder_level1/fastq/$DATA.0/$i /ephemeral/fastq/NA12878_XXX/$i --dryrun"
  aws s3 cp aws_s3_path_for_input_fastq_e.g.s3://bucket/folder_level1/fastq/$DATA.0/$i /ephemeral/fastq/NA12878_XXX/$i
done
ls -l /ephemeral/fastq/NA12878_XXX | grep fastq.gz | \
  perl -e '@_=<>; print map{chop(); @t=split(/\s+/,$_); "$t[8]\t$t[4]\n" } @_' | sort > pull-fastq/fastq-download.list

if diff -b pull-fastq/fastq.list pull-fastq/fastq-download.list | grep -q . ; then echo "Error, fastq size don't match"; else echo "OK" > pull-fastq/fastq.ok; fi


if ! [ -s pull-fastq/fastq.csv ]; then echo "zero size pull-fastq/fastq.csv"; exit; fi
if ! [ -s pull-fastq/fastq.ok ]; then echo "zero size pull-fastq/fastq.ok"; exit; fi

date +%s > PATH/hg-workflow/NA12878_XXX/pull-fastq/WF.complete.date
my_time_end=`date +%s`;
my_time_spent=$((my_time_end-my_time_start))
echo "sample=NA12878_XXX job=pull-fastq host=$my_host pid=$my_pid queue=$my_queue cores=$my_core time_start=$my_time_start time_end=$my_time_end time_spent=$my_time_spent" >> PATH/hg-workflow/NA12878_XXX/pull-fastq/WF.cpu

