#!/usr/bin/env python3

import re
import os
import sys
import argparse
from argparse import RawTextHelpFormatter

import numpy as np
import pandas as pd
import json
import logging
import boto3
from botocore.exceptions import ClientError

__author__ = 'Weizhong Li'

#### global variables
json_s3_source = False
rec = []

def process_FC_stat_json_list(input_list):
  with open(input_list, 'r') as f:
    for line in f:
      if line[0] == '#':
        continue
      if not re.match('^\w', line):
        continue
      ll = re.split('\s+', line.rstrip());
      process_FC_stat_json(ll[0])
    f.close()

def process_FC_stat_json(input_file):
  json_file = input_file
  global json_s3_source

  s3_fq_info = []
  ## aws s3 object s3://bucket_name/folder1/folder1_1/RUN_ID/Stats/Stats.json
  if input_file[0:5] == "s3://":
    json_s3_source = True
    json_file = 'tmp' + str(os.getpid()) + '.json'
    bucket1, obj1 = re.split('\/', input_file[5:], maxsplit=1)
    s3_prefix = re.sub('\/Stats\/Stats.json$','', obj1)  ## to be folder1/folder1_1/RUN_ID
    if s3_prefix == obj1:
      s3_prefix = '/'

    s3 = boto3.client('s3')
    with open(json_file, 'wb') as f:
      s3.download_fileobj(bucket1, obj1, f)
    f.close()

    res = s3.list_objects_v2(Bucket=bucket1, Prefix=s3_prefix, MaxKeys=10000)
    for i in range(0, len(res['Contents'])):
      obj_i = res['Contents'][i]
      path1 = obj_i['Key']
      if not re.search('\.fastq\.gz$', path1): continue

      ll = re.split('\/', path1)
      file1 = ll.pop()
      ll = re.split('_', file1)
      sample = ll[0]
      if re.search('_R1_001.fastq.gz$', file1): R = "R1" 
      else: R = "R2"
      lane = int(re.findall('_S\d_L00\d_',file1)[0][7:8])
      size1 = obj_i['Size']
      s3_fq_info.append([sample, lane, R, file1, size1])

  with open(json_file) as f:
    data = json.load(f)
  f.close()

  fc_id = data['Flowcell']
  run_id = data['RunId']

  for i in range(0,len(data['ConversionResults'])):
    lane = int(data['ConversionResults'][i]['LaneNumber'])
    num_reads = int(data['ConversionResults'][i]['TotalClustersPF'])

    for j in range(0, len( data['ConversionResults'][i]['DemuxResults'])):
      dataj = data['ConversionResults'][i]['DemuxResults'][j]
      sample0 = dataj['SampleId']
      sample = dataj['SampleId']
      if sample == 'NA12878':
        sample = sample + '_' + fc_id
      sample_reads = int(dataj['NumberReads'])
      if sample_reads <= 0: continue
      reads_idx_mis_0 = int(dataj['IndexMetrics'][0]['MismatchCounts']['0'])
      reads_idx_mis_1 = int(dataj['IndexMetrics'][0]['MismatchCounts']['1'])
      perfect_idx = 1.0 * reads_idx_mis_0 / (0.0 + reads_idx_mis_0 + reads_idx_mis_1)
      yield_1       = int(dataj['ReadMetrics'][0]['Yield'])
      yield_1_Q30   = int(dataj['ReadMetrics'][0]['YieldQ30'])
      trimmed_1     = int(dataj['ReadMetrics'][0]['TrimmedBases'])
      yield_2       = int(dataj['ReadMetrics'][1]['Yield'])
      yield_2_Q30   = int(dataj['ReadMetrics'][1]['YieldQ30'])
      trimmed_2     = int(dataj['ReadMetrics'][1]['TrimmedBases'])
      num_bases     = yield_1+yield_2
      num_bases_post_trimming = num_bases - trimmed_1 - trimmed_2
      depth         = 1.0 * num_bases_post_trimming / 3000000000.
      Q30 = (0.0 + yield_1_Q30 + yield_2_Q30) / (0.0 + yield_1 + yield_2)

      if json_s3_source: ## with fq_info 
        R1_f = ''
        R2_f = ''
        R1_size = ''
        R2_size = ''
        for k in range(0, len(s3_fq_info)):
          ki = s3_fq_info[k]
          if not ( ki[0] == sample0 and ki[1] == lane): continue
          if ki[2] == "R1":
            R1_f = ki[3]
            R1_size = ki[4]
          elif ki[2] == "R2":
            R2_f = ki[3]
            R2_size = ki[4]
        if R1_f and R2_f: ## skip records without R1 & R2
          rec.append([sample, fc_id, run_id, lane, sample_reads, num_bases, num_bases_post_trimming, depth, perfect_idx, Q30, R1_f, R1_size, R2_f, R2_size])
        else:
          logging.warning('no fastq file\t{0}\t{1}\t{2}'.format(sample, lane, run_id))
      else: ## without fq_info
        rec.append([sample, fc_id, run_id, lane, sample_reads, num_bases, num_bases_post_trimming, depth, perfect_idx, Q30, R1_f, R1_size, R2_f, R2_size])

def output_rec(output_file1, output_file2):
  rec.sort(key=lambda x: (x[2], x[0], x[3]), reverse=False) ## by run_id, sample, lane
  fout = open(output_file1, 'w')

  if json_s3_source: ## with fq_info
    fout.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\n'.format('sample', \
      'fc_id', 'run_id', 'lane', 'num_reads', 'num_bases', 'num_bases_post_trimming', 'depth', 'frac_perfect_index','frac_Q30','R1','R1_file_size','R2','R2_file_size'))
    for i in range(0, len(rec)):
      reci = rec[i]
      fout.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\n'.format(reci[0], \
        reci[1], reci[2], reci[3], reci[4], reci[5], reci[6], reci[7], reci[8], reci[9], reci[10], reci[11], reci[12], reci[13] ))
  else:
    fout.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n'.format('sample', \
      'fc_id', 'run_id', 'lane', 'num_reads', 'num_bases', 'num_bases_post_trimming', 'depth', 'frac_perfect_index','frac_Q30'))
    for i in range(0, len(rec)):
      reci = rec[i]
      fout.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n'.format(reci[0], \
        reci[1], reci[2], reci[3], reci[4], reci[5], reci[6], reci[7], reci[8], reci[9] ))
  fout.close()

  #### merge by sample
  if json_s3_source: ## with fq_info
    df = pd.DataFrame(rec, columns = ['sample', 'fc_id', 'run_id', 'lane', 'num_reads', 'num_bases', 'num_bases_post_trimming', 'depth', 'frac_perfect_idx', 'frac_Q30', 'z1','z2','z3','z4'])
  else:
    df = pd.DataFrame(rec, columns = ['sample', 'fc_id', 'run_id', 'lane', 'num_reads', 'num_bases', 'num_bases_post_trimming', 'depth', 'frac_perfect_idx', 'frac_Q30' ])
  df['num_reads_perfect'] = (df['num_reads'] * df['frac_perfect_idx'])
  df['num_reads_Q30'] = df['num_reads'] * df['frac_Q30']
  df_sample = df[['fc_id', 'run_id', 'sample', 'num_reads', 'num_bases', 'num_bases_post_trimming', 'depth', 'num_reads_perfect','num_reads_Q30']].groupby(['sample','fc_id', 'run_id'], as_index=False).sum()
  df_sample['frac_perfect_idx'] = 1.0 * df_sample['num_reads_perfect'] / df_sample['num_reads']
  df_sample['frac_Q30'] = 1.0 * df_sample['num_reads_Q30'] / df_sample['num_reads']
  df_sample = df_sample[['sample', 'fc_id', 'run_id', 'num_reads', 'num_bases', 'num_bases_post_trimming', 'depth', 'frac_perfect_idx', 'frac_Q30']]
  df_sample.to_csv(output_file2, index=False, sep='\t')


if __name__ == "__main__":
  parser = argparse.ArgumentParser(formatter_class = RawTextHelpFormatter,
                                   description     = 'Process flowcell demultiplexing stat json file')

  parser.add_argument('-i', '--input',       required=False, help='''Input AWS s3 URL or Path of a FC/Stats/Stats.json file.
  ''')

  parser.add_argument('-f', '--list_files',       required=False, help='''Filename of a list of input json files.
  ''')

  parser.add_argument('-o', '--output',       required=True, help='''Output filename prefix, required. 
output file 1: prefix.tsv, per sample / lane statistics
output file 2: perfix-sample.tsv, per sample statistics
output file 3: perfix.log, log file
  ''')
  args = parser.parse_args()

  out_file1 = args.output + '.tsv'
  out_file2 = args.output + '-sample.tsv'
  out_log1  = args.output + '.log'

  logging.basicConfig(format='%(levelname)s:%(message)s', filename=out_log1)

  if args.input:
    process_FC_stat_json(args.input)
    output_rec(out_file1, out_file2)
  elif args.list_files:
    process_FC_stat_json_list(args.list_files)
    output_rec(out_file1, out_file2)
  else:
    print("Error: no input")
    exit(1)



