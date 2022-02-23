#!/usr/bin/env python3

import sys
import re
import argparse
import pandas as pd

qc_stat = []
time_stat = []
var_stat = []

def process_sample(sample):
  ## time
  with open(sample + '/hg-call/HG.time_metrics.csv', 'r') as f:
    t_sv, t_align, t_output, t_total = 0, 0, 0, 0
    for line in f:
      ll = re.split(',', line.rstrip())
      if ll[2] == 'Time aligning reads':
        t_align = ll[4]
      elif ll[2] == 'Time saving map/align output':
        t_output =ll[4]
      elif ll[2] == 'Time structural variant calling':
        t_sv = ll[4]
      elif ll[2] == 'Total runtime':
        t_total = ll[4]
    f.close()
    time_stat.append([sample, t_sv, t_align, t_output, t_total])

  ## variant
  n_variant = 0 
  muti_allelic = 0
  snp = 0
  ins = 0
  dele = 0
  indel = 0
  with open(sample + '/hg-call/HG.vc_metrics.csv', 'r') as f:
    for line in f:
      ll = re.split(',', line.rstrip())
      if   ll[2] == 'Total':        n_variant = ll[3]
      elif ll[2] == 'Multiallelic': muti_allelic = ll[4]
      elif ll[2] == 'SNPs':         snp = ll[3]
      elif ll[2] == 'Insertions (Hom)': ins = ins + int(ll[3])
      elif ll[2] == 'Insertions (Het)': ins = ins + int(ll[3])
      elif ll[2] == 'Deletions (Hom)': dele = dele + int(ll[3])
      elif ll[2] == 'Deletions (Het)': dele = dele + int(ll[3])
      elif ll[2] == 'Indels (Het)': indel = ll[3]
    f.close()

  cnv = 0
  cnv_del = 0
  cnv_ins = 0
  cnv_del_pass = 0
  cnv_ins_pass = 0
  with open(sample + '/hg-call/HG.cnv_metrics.csv', 'r') as f:
    for line in f:
      ll = re.split(',', line.rstrip())
      if   ll[2] == 'Number of amplifications': cnv_ins = ll[3]
      elif ll[2] == 'Number of deletions':      cnv_del = ll[3]
      elif ll[2] == 'Number of passing amplifications': cnv_ins_pass = int(ll[3])
      elif ll[2] == 'Number of passing deletions': cnv_del_pass = int(ll[3])
    f.close()
    cnv = cnv_del_pass + cnv_ins_pass

  sv = 0
  sv_del = 0
  sv_ins = 0
  sv_dup = 0
  sv_bnd = 0
  with open(sample + '/hg-call/HG.sv_metrics.csv', 'r') as f:
    for line in f:
      ll = re.split(',', line.rstrip())
      if   ll[2] == 'Number of deletions (PASS)':  sv_del = int(ll[3])
      elif ll[2] == 'Number of insertions (PASS)': sv_ins = int(ll[3])
      elif ll[2] == 'Number of duplications (PASS)': sv_dup = int(ll[3])
      elif ll[2] == 'Number of breakend pairs (PASS)': sv_bnd = int(ll[3])
    f.close()
    sv = sv_ins + sv_del + sv_dup + sv_bnd

  var_stat.append([sample, n_variant, muti_allelic, snp, ins, dele, indel, cnv, cnv_ins, cnv_del, cnv_ins_pass, cnv_del_pass, sv, sv_ins, sv_del, sv_dup, sv_bnd])

  ## QC
  ploidy = ''
  with open(sample + '/hg-call/HG.ploidy_estimation_metrics.csv', 'r') as f:
    for line in f:
      ll = re.split(',', line.rstrip())
      if   ll[2] == 'Ploidy estimation':  ploidy = ll[3]
    f.close()

  depth = 0
  with open(sample + '/hg-call/HG.wgs_overall_mean_cov.csv', 'r') as f:
    for line in f:
      ll = re.split(',', line.rstrip())
      if   ll[0] == 'Average alignment coverage over wgs':  depth = ll[1]
    f.close()

  num_reads = 0
  num_bases = 0
  read_len = 0
  with open(sample + '/hg-call/HG.trimmer_metrics.csv', 'r') as f:
    for line in f:
      ll = re.split(',', line.rstrip())
      if   ll[2] == 'Total input reads':  num_reads = ll[3]
      elif ll[2] == 'Total input bases':  num_bases = ll[3]
      elif ll[2] == 'Average input read length':  read_len = ll[3]
    f.close()

  mapped_reads = 0
  mapped_R1 = 0
  mapped_R2 = 0
  Q30 = 0
  Q30_R1 = 0
  Q30_R2 = 0
  ins_size_mean = 0
  ins_size_median = 0
  with open(sample + '/hg-call/HG.mapping_metrics.csv', 'r') as f:
    for line in f:
      ll = re.split(',', line.rstrip())
      if   ll[0] == 'MAPPING/ALIGNING PER RG': continue
      if   ll[2] == 'Mapped reads':  mapped_reads = ll[4]
      elif ll[2] == 'Mapped reads R1':  mapped_R1 = ll[4]
      elif ll[2] == 'Mapped reads R2':  mapped_R2 = ll[4]
      elif ll[2] == 'Q30 bases':  Q30 = ll[4]
      elif ll[2] == 'Q30 bases R1':  Q30_R1  = ll[4]
      elif ll[2] == 'Q30 bases R2':  Q30_R2  = ll[4]
      elif ll[2] == 'Insert length: mean':  ins_size_mean = ll[3]
      elif ll[2] == 'Insert length: median':  ins_size_median = ll[3]
    f.close()

  qc_stat.append([sample, ploidy, depth, num_reads, num_bases, read_len, mapped_reads, mapped_R1, mapped_R2, Q30, Q30_R1, Q30_R2, ins_size_mean, ins_size_median])



def process_all(sample_file, output_prefix):
  with open(sample_file, 'r') as f:
    for line in f:
      ll = re.split('\s+', line.rstrip())
      process_sample(ll[0])
    f.close()

  ## time
  df = pd.DataFrame(time_stat, columns=['sample','time_sv','time_map','time_output','time_total'])
  df.to_csv(output_prefix + '.time.tsv', index=False, sep='\t')

  ## variants
  df = pd.DataFrame(var_stat, columns=['sample', 'n_variant','muti_allelic','snp','ins','del','indel', \
                                       'cnv', 'cnv_ins', 'cnv_del', 'cnv_ins_pass', 'cnv_del_pass', \
                                       'sv', 'sv_ins', 'sv_del', 'sv_dup', 'sv_bnd'])
  df.to_csv(output_prefix + '.variants.tsv', index=False, sep='\t')

  ## QC
  df = pd.DataFrame(qc_stat, columns=['sample', 'ploidy', 'depth', 'num_reads', 'num_bases', 'read_len', \
                                      'mapped_reads', 'mapped_R1', 'mapped_R2', 'Q30', 'Q30_R1', 'Q30_R2', 'ins_size_mean', 'ins_size_median'])
  df.to_csv(output_prefix + '.qc.tsv', index=False, sep='\t')


def main():
  parser = argparse.ArgumentParser(description     = 'Post Dragen summary')

  parser.add_argument('-i', '--input',       required=False, help='''Input NGS-samples file''')
  parser.add_argument('-o', '--output',      required=True, help='''Output filename prefix, required. 
output file 1: prefix.time.tsv, cpu time results
output file 2: perfix.qc.tsv, QC results
output file 3: perfix.var.tsv, variant results
  ''')
  args = parser.parse_args()

  if args.input and args.output:
    process_all(args.input, args.output)
  else:
    print("Error: no input")
    exit(1)


if __name__ == "__main__":
  main()

