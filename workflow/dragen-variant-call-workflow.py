#!/usr/bin/python

################################################################################
# NGS workflow by Weizhong Li, http://weizhongli-lab.org
################################################################################

queue_system = 'SGE'

########## local variables etc. Please edit
ENV={
  'NGS_root' : '/ephemeral/data',
}

########## computation resources for execution of jobs
NGS_executions = {}

NGS_executions['sh0'] = {
  'type'                : 'sh',
  'cores_per_node'      : 12,
  'number_nodes'        : 1,
  'template'            : '''#!/bin/bash

'''
}


NGS_executions['sh1'] = {
  'type'                : 'sh',
  'cores_per_node'      : 32,
  'number_nodes'        : 1,
  'template'            : '''#!/bin/bash

'''
}


NGS_batch_jobs = {}

NGS_batch_jobs['pull-fastq'] = {
  ## CMD_opts
  ## 1st option: aws s3 folder of input fastq files
  ## 2nd option: path of local aws instance to save fastq files
  ## 3rd option: aws s3 folder of flowcell data
  ## 4th option: depth of reads cutoff (GB), if this sample has fastq files from multiple flowcell lanes and total reads 
  ##             is > than the cutoff, only a subset of lanes will be used for mapping to save cost
  'CMD_opts'         : ['aws_s3_path_for_input_fastq_e.g.s3://bucket/folder_level1/fastq','/ephemeral/fastq','aws_s3_path_for_flowcell_info_tsv_file_e.g.s3://bucket/folder_level1/flowcell_all.tsv','40'],
  'non_zero_files'   : ['fastq.csv','fastq.ok'],
  'execution'        : 'sh0',               # where to execute
  'cores_per_cmd'    : 4,                    # number of threads used by command below
  'no_parallel'      : 1,                    # number of total jobs to run using command below
  'command'          : '''

#### sleep up to 10 s
sleep $(( $RANDOM % 10  ))
#### wait unless there is less than 2 folderis in $CMDOPTS.1/
#### so that fastq files for all samples to be downloaded be <=3 
#### fastq folder will be deleted in the next job
if [ ! -e $CMDOPTS.1 ]
then
  mkdir -p $CMDOPTS.1
fi

num=$(ls -1 $CMDOPTS.1 | grep -c .)
while [ $num -gt 2 ]; do
  echo wait 
  sleep 60
  num=$(ls -1 $CMDOPTS.1 | grep -c .)
done
mkdir -p $CMDOPTS.1/$SAMPLE

aws s3 cp $CMDOPTS.2 $SELF/sample-all.tsv
grep $SAMPLE $SELF/sample-all.tsv | grep $DATA.0 > $SELF/sample.tsv
rm -f $SELF/sample-all.tsv

if [ ! -s $SELF/sample.tsv ]
then
  exit
fi

#### embeded python to select lanes up to 40GB
#### to be run as: env python job1.py $SAMPLE  $CMDOPTS.3
cat << ENDPY > job1.py
import re
import sys

input_sample_name = sys.argv[1]
target_depth = float(sys.argv[2])
f = open('$SELF/sample.tsv')
f_csv = open('$SELF/fastq.csv', 'w')
fout = open('$SELF/fastq.list', 'w')
f_csv.write('RGID,RGSM,RGLB,Lane,Read1File,Read2File\\n')

depth_total = 0
for line in f:
  ll = re.split('\\t', line.rstrip())
  depth_total = depth_total + float(ll[7])
  rgid = '$SAMPLE' + 'L' + str(ll[3])
  #### write fastq.csv
  f_csv.write('{0},{1},{2},{3},$CMDOPTS.1/$SAMPLE/{4},$CMDOPTS.1/$SAMPLE/{5}\\n'.format(rgid, '$SAMPLE', 'LIB1', ll[3], ll[10], ll[12] ))

  fout.write('{0}\\t{1}\\n'.format(ll[10], ll[11]))
  fout.write('{0}\\t{1}\\n'.format(ll[12], ll[13]))
  if depth_total > target_depth: break

f.close()
f_csv.close()
fout.close()
ENDPY
env python job1.py $SAMPLE $CMDOPTS.3 
sort $SELF/fastq.list > $SELF/fastq.list.1
mv -f $SELF/fastq.list.1 $SELF/fastq.list

#### download fastqs
for i in $(cut -f 1 $SELF/fastq.list | sort) 
  do
  echo "aws s3 cp $CMDOPTS.0/$DATA.0/$i $CMDOPTS.1/$SAMPLE/$i --dryrun"
  aws s3 cp $CMDOPTS.0/$DATA.0/$i $CMDOPTS.1/$SAMPLE/$i
done
ls -l $CMDOPTS.1/$SAMPLE | grep fastq.gz | \\
  perl -e '@_=<>; print map{chop(); @t=split(/\\s+/,$_); "$t[8]\\t$t[4]\\n" } @_' | sort > $SELF/fastq-download.list

if diff -b $SELF/fastq.list $SELF/fastq-download.list | grep -q . ; then echo "Error, fastq size don't match"; else echo "OK" > $SELF/fastq.ok; fi

'''
}

NGS_batch_jobs['hg-call'] = {
  ## CMD_opts
  ## 1st option: path of local aws instance to save fastq files
  'CMD_opts'         : ['/ephemeral/fastq'],
  'injobs'           : ['pull-fastq'],
  'non_zero_files'   : ['HG.cnv.vcf.gz','HG.sv.vcf.gz','HG.hard-filtered.gvcf.gz'],
  'execution'        : 'sh1',               # where to execute
  'cores_per_cmd'    : 18,                    # number of threads used by command below
  'no_parallel'      : 1,                    # number of total jobs to run using command below
  'command'          : '''

#### ver 1, bam + gvcf + sv + cnv - with gvcf, there is no un-hard-filtered gvcf output
dragen -f -r $ENV.NGS_root/refs/GRch38cnv \\
    --fastq-list $INJOBS.0/fastq.csv \\
    --fastq-list-sample-id $SAMPLE \\
    --lic-server=0xU1SGJKuGc=:NUH56PNZxxEIgPDmEChSGcGJYEgbY2mY@license.edicogenome.com \\
    --enable-map-align true \\
    --enable-map-align-output true \\
    --enable-variant-caller true \\
    --vc-emit-ref-confidence GVCF \\
    --enable-sv true \\
    --enable-cnv true \\
    --cnv-enable-self-normalization true \\
    --output-directory $SELF \\
    --output-file-prefix HG > $SELF/WF.log

if [[ -s $SELF/HG.cnv.vcf.gz && -s $SELF/HG.sv.vcf.gz && -s $SELF/HG.hard-filtered.gvcf.gz ]]; then 
  #### remove fastq files
  rm -rf $CMDOPTS.0/$SAMPLE
  rm -f  $SELF/HG.bam*
  rm -rf $SELF/sv/workspace/genomeDepth
  rm -rf $SELF/sort_spill
  rm -rf $SELF/vc
fi


'''
}


NGS_batch_jobs['vcf-upload'] = {
  'CMD_opts'         : ['aws_s3_path_for_vcf_and_results_upload_e.g.s3://bucket/folder_level1/folder_level2'],
  'injobs'           : ['pull-fastq','hg-call'],
  'non_zero_files'   : ['upload.list'],
  'execution'        : 'sh1',               # where to execute
  'cores_per_cmd'    : 4,                    # number of threads used by command below
  'no_parallel'      : 1,                    # number of total jobs to run using command below
  'command'          : '''

for i in HG.hard-filtered.gvcf.gz HG.cnv.vcf.gz HG.sv.vcf.gz; do
  echo "s3 cp $INJOBS.1/$i $CMDOPTS.0/$SAMPLE/$i"
  echo "s3 cp $INJOBS.1/$i.tbi $CMDOPTS.0/$SAMPLE/$i.tbi"
  aws s3 cp $INJOBS.1/$i $CMDOPTS.0/$SAMPLE/$i
  aws s3 cp $INJOBS.1/$i.tbi $CMDOPTS.0/$SAMPLE/$i.tbi
done 

for i in HG-replay.json WF.log; do
  echo "s3 cp $INJOBS.1/$i $CMDOPTS.0/$SAMPLE/$i"
  aws s3 cp $INJOBS.1/$i $CMDOPTS.0/$SAMPLE/$i
done 

tar cvf work-files.tar.gz --exclude=$INJOBS.1/*vcf.gz* --gzip $INJOBS.0 $INJOBS.1
aws s3 cp work-files.tar.gz $CMDOPTS.0/$SAMPLE/work-files.tar.gz

aws s3 ls $CMDOPTS.0/$SAMPLE/ > $SELF/upload.list

rm -f $INJOBS.1/HG.hard-filtered.gvcf.gz

'''
}


#### skip vcf-ann, very slow
#### will do vcf-ann after gvcf merge and filtering
NGS_batch_jobs['vcf-ann'] = {
  'CMD_opts'         : [],
  'injobs'           : ['hg-call'],
  'non_zero_files'   : ['HG.Nirvana.json.gz'],
  'execution'        : 'sh1',               # where to execute
  'cores_per_cmd'    : 6,                    # number of threads used by command below
  'no_parallel'      : 1,                    # number of total jobs to run using command below
  'command'          : '''

# Nirvana
/opt/edico/share/nirvana/Nirvana -c $ENV.NGS_root/refs/nirvana/Cache/GRCh38/Both \\
  -r $ENV.NGS_root/refs/nirvana/References/Homo_sapiens.GRCh38.Nirvana.dat \\
  --sd $ENV.NGS_root/refs/nirvana/SupplementaryAnnotation/GRCh38 -i $INJOBS.0/HG.vcf.gz -o $SELF/HG.Nirvana

'''
}






