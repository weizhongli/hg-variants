#!/usr/bin/env python3

import sys
import re
import argparse
import gzip

## global variables
vcf_in = None
vcf_out = None
vcf_headers = [] ## #CHROM  POS ID  REF ALT QUAL  FILTER  INFO  FORMAT
vcf_samples = []


def snp_indel_test(ref1, alt1):
  '''given ref allele and alt allele, test if this is a SNP or INDEL
  after vcf merging (e.g. by GLnexus), SNP may look at indel
  chr1_12345_T_A	TTTTATTTA	TTTTATTAA  SNP        => T to A
  '''
  if len(ref1) != len(alt1): return ['indel', 0]
  pos, num_diff = 0, 0
  for i in range(len(ref1)):
    if ref1[i] != alt1[i]:
      pos = i
      num_diff = num_diff + 1
  return (['snp', pos] if num_diff == 1 else ['indel', 0])
 

class vcf_record:
  ## class variables, shared by all instances
  allowed_gt_fields = ['GT', 'DP', 'AD', 'GQ', 'RNC']
  num_sample = 0

  def __init__(self, line):
    #    0    1  2    3   4    5       6     7       8        9       10
    #CHROM  POS ID  REF ALT QUAL  FILTER  INFO  FORMAT  sample1  sample2
    lls = re.split('\t', line.rstrip())
    self.chrom, self.pos, self.id, self.ref, self.alt, self.qual, self.filter, self.info, self.format = lls[0:9]
    self.samples    = lls[9:]
    self.num_sample = len(self.samples)     ## class variables
    self.alt_lst    = re.split(',', self.alt)
    self.format_lst = re.split(':', self.format) ## e.g. GT:DP:AD:GQ:PL:RNC
    self.alt_type = []
    self.alt_pos  = []
    self.num_alt_type = 0


  def allele_typing(self):
    for alt in self.alt_lst:
      t_type, t_pos = snp_indel_test(self.ref, alt)
      self.alt_type.append(t_type)
      self.alt_pos.append(t_pos)
    self.num_alt_type = len( list(set(self.alt_type)))


  def reduce_gt_fields(self):
    '''this script can only handle some genotype fields, remove others'''
    allowed_gt_idx = [ x for x in range(len(self.format_lst)) if self.format_lst[x] in self.allowed_gt_fields]
    if len(allowed_gt_idx) == len(self.format_lst): return

    for i in range(self.num_sample):
      fields = re.split(":", self.samples[i])
      self.samples[i] = str.join(':', [ fields[x] for x in allowed_gt_idx ] )
    self.format_lst = [ self.format_lst[x] for x in allowed_gt_idx ]
    self.format = str.join(':', self.format_lst)

    
  def merge_indel(self):
    '''merge different indel alleles into a single indel allele, re-assign GT'''
    if len(self.alt_type) == 0: self.allele_typing()
    num_alt = len(self.alt_lst)

    self.reduce_gt_fields()
    if num_alt == 1: return # only 1 alt allele
    if not ('indel' in self.alt_type): return # no indel at this site
    if len([x for x in self.alt_type if x == 'indel']) == 1: return # only 1 indel at this site

    first_indel_idx = [x for x in range(num_alt) if self.alt_type[x] == 'indel'][0]
    new_alt_idx = [x for x in range(num_alt) if self.alt_type[x] != 'indel']
    new_alt_idx.append(first_indel_idx)
    new_alt_idx.sort()

    new_idx = []
    j = 0
    for i in range(num_alt):
      if self.alt_type[i] == 'indel':
        new_idx.append(first_indel_idx)
        if i == first_indel_idx: j = j+1
      else:
        new_idx.append(j)
        j = j+1
    # now new_idx[i] is new idx of allele i
    num_alt_new = j

    # update self.alt and self.alt_lst
    self.alt_lst = [self.alt_lst[x] for x in new_alt_idx]
    self.alt = str.join(',', self.alt_lst)

    # update info - if having AF= and AQ=, such as by GLnexus
    info_lst = re.split(';', self.info)
    for i in range(len(info_lst)):
      f1 = info_lst[i]
      if re.search('AF=', f1):
        old_AFs = re.split(',', f1[3:])
        old_AFs = [float(x) for x in old_AFs]
        new_AFs = [0.0 for x in range(num_alt_new) ] ## list of 0.0
        for k in range(num_alt):
          m = new_idx[k]
          new_AFs[m] = new_AFs[m] + old_AFs[k]
        f1 = 'AF=' + str.join(',', [ str(x) for x in new_AFs])
        info_lst[i] = f1
      elif re.search('AQ=', f1):
        old_AQs = re.split(',', f1[3:])
        old_AQs = [int(x) for x in old_AQs]
        new_AQs = [0 for x in range(num_alt_new) ] ## list of 0
        for k in range(num_alt):
          m = new_idx[k]
          new_AQs[m] = new_AQs[m] + old_AQs[k]
        f1 = 'AQ=' + str.join(',', [ str(x) for x in new_AQs])
        info_lst[i] = f1
    self.info = str.join(';', info_lst)

    # update id - if id is the join of all alts, such as by GLnexus, then use the join of new alts
    ids = re.split(';', self.id)
    if len(ids) == num_alt:
      ids = [ ids[x] for x in new_alt_idx ]
      self.id = str.join(';', ids) 

    # now update sample genotypes data
    for i in range(self.num_sample):
      ffs = re.split(':', self.samples[i])
      new_v = []
      # GT:DP:AD:GQ:PL:RNC  0/0:23:12,0,0:0:0,0,114,0,114,114:..      see VCFv4.2.pdf page 6/28
      #                      GT:DP:    AD:GQ:              PL:RNC
      for j in range(len(ffs)):
        i_type = self.format_lst[j]

        ## GT : genotype, encoded as allele values separated by /
        if i_type == "GT":
          old_gt = re.split('\/|\|', ffs[j])    # such as ['.', '2'] or ['0','1']
          #new_gt = ['.' if x== '.' else ( 0 if x == 0 else new_idx[int(x)-1]+1 ) for x in old_gt]
          new_gt = [ x if x in ['.','0'] else new_idx[int(x)-1] + 1 for x in old_gt]
          new_gt = [str(x) for x in new_gt] # convert to str
          new_v.append(str.join('/', new_gt))
        
        ## DP : read depth at this position for this sample (Integer)
        elif i_type == "DP":
          new_v.append(ffs[j])
        
        ## AD : Allelic depths for the ref and alt alleles in the order listed
        ## can be a '.' or a list of int (or .)
        ## if list of mixture of int and .
        ## . will be converted into 0 
        elif i_type == "AD":
          if ffs[j] == '.':
            new_v.append(ffs[j])
            continue

          if re.search('\.', ffs[j]): ffs[j] = re.sub('\.', '0', ffs[j])
          old_ad = re.split(',', ffs[j])
          old_ad = [ int(x) for x in old_ad]
          new_ad = [0 for x in range(num_alt_new+1)] # list of 0
          new_ad[0] = old_ad[0] # first element is for ref allele, copy over

          for k in range(num_alt):
            m = new_idx[k - 1] + 1
            new_ad[m] = new_ad[m] + old_ad[k]
          new_v.append(str.join(',', [str(x) for x in new_ad] ))
        
        ## GQ : conditional genotype quality
        elif i_type == "GQ":
          new_v.append(ffs[j])

        ## skip PL, it is too complicated to process
        elif i_type == "PL":
          continue

        ## RNC: Reason for No Call in GT (not standard in VCF)
        elif i_type == "RNC":
          new_v.append(ffs[j])

        ## skip everythig else  
        else:
          continue
      self.samples[i] = str.join(':', new_v)
    
    self.alt_type = []
    self.alt_pos  = []
    self.num_alt_type = 0
    self.allele_typing()

  
  def write(self, fout):
    fout.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t'.format(self.chrom, self.pos, \
      self.id, self.ref, self.alt, self.qual, self.filter, self.info, self.format))
    fout.write(str.join('\t', self.samples) + '\n')

  def split_INFO(self, n_alt):
    '''return a list of INFO, currently only handle AF, AQ'''
    ## AF=0.187123,0.01509;AQ=149,197;XX=yyy -> [AF=0.187123;AQ=149;XX=yyy, AF=0.01509;AQ=197;XX=yyy]
    if n_alt == 1:
      return [self.info]
  
    split_info = ['' for x in range(n_alt)]
    n_alt_match_flag = True
    ff = re.split(';',self.info) ## [AF=0.187123,0.01509 AQ=149,197 XX=yyy]
    for i in ff:
      if re.search('=', i):
        key1, value1 = re.split('=', i)
        if key1 in ['AF','AQ']:
          vv = re.split(',', value1)
          if len(vv) != n_alt:
            n_alt_match_flag = False
            break
          split_info = [ key1+'='+vv[x]  if len(split_info[x])==0 else split_info[x]+';'+key1+'='+vv[x]  for x in range(n_alt) ]
        else:
          split_info = [ i               if len(split_info[x])==0 else split_info[x]+';'+i               for x in range(n_alt) ]
      else:
        split_info = [ i               if len(split_info[x])==0 else split_info[x]+';'+i               for x in range(n_alt) ]
  
    if n_alt_match_flag:
      return split_info
    else:
      return [ self.info for x in range(n_alt)]


  def write_multi(self, fout, args):
    '''if a site has multi alleles, output them into seperate lines'''
    num_alt = len(self.alt_lst)
    if num_alt == 1: 
      self.write(fout)
      return


    sorted_alt_idx = list(range(num_alt))
    sorted_alt_idx.sort(key=lambda x: self.alt_pos[x])
    list_id  = re.split(';', self.id)
    info_split = self.split_INFO(num_alt)
    for idx in sorted_alt_idx:
      idx_1 = str(idx+1) ## +1 since 0 is ref
      t_pos = self.alt_pos[idx]+int(self.pos) if len(self.alt_pos) == num_alt else self.alt_pos
      t_id  = list_id[idx]      if len(list_id)      == num_alt else self.id
      t_ref = self.ref
      t_alt = self.alt_lst[idx]
      if (self.alt_type[idx] == 'snp') and (len(t_ref)>1):
        t_ref = t_ref[ self.alt_pos[idx] ]
        t_alt = t_alt[ self.alt_pos[idx] ]
      if args.GLnexus_ID: ##chr1_12345_T_A
        id_split = re.split('_', t_id)
        if len(id_split) == 4:
          if (id_split[0] == self.chrom) and (re.fullmatch('\d+', id_split[1])) and \
            (re.fullmatch('[ATCGatcg]+', id_split[2])) and (re.fullmatch('[ATCGatcg]+', id_split[3])):
            t_pos = int(id_split[1])
            t_ref = id_split[2]
            t_alt = id_split[3]

      fout.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}'.format(self.chrom, t_pos, \
        t_id, t_ref, t_alt, self.qual, self.filter, info_split[idx], self.format))

      for i in range(self.num_sample):
        ffs = re.split(':', self.samples[i])
        new_v = []
        # GT:DP:AD:GQ:PL:RNC  0/0:23:12,0,0:0:0,0,114,0,114,114:..      see VCFv4.2.pdf page 6/28
        #                      GT:DP:    AD:GQ:              PL:RNC
        for j in range(len(ffs)):
          i_type = self.format_lst[j]

          ## GT : genotype, encoded as allele values separated by /
          if i_type == "GT":
            old_gt = re.split('\/|\|', ffs[j])    # such as ['.', '2'] or ['0','1']
            new_gt = [ x if x in ['0','.']  else ( '1' if x==idx_1 else '.' )  for x in old_gt]
            new_v.append(str.join('/', new_gt))
        
          ## DP : read depth at this position for this sample (Integer)
          elif i_type == "DP":
            new_v.append(ffs[j])
        
          ## AD : Allelic depths for the ref and alt alleles in the order listed
          ## can be a '.' or a list of int (or .)
          ## if list of mixture of int and .
          ## . will be converted into 0 
          elif i_type == "AD":
            if ffs[j] == '.':
              new_v.append(ffs[j])
              continue
            old_ad = re.split(',', ffs[j])
            if len(old_ad) < num_alt+1:
              new_v.append('.')
              continue
            new_ad = [ old_ad[x] for x in [0, idx+1]]
            new_v.append(str.join(',', new_ad))
       
          ## GQ : conditional genotype quality
          elif i_type == "GQ":
            new_v.append(ffs[j])

          ## skip PL, it is too complicated to process
          elif i_type == "PL":
            continue

          ## RNC: Reason for No Call in GT (not standard in VCF)
          elif i_type == "RNC":
            new_v.append(ffs[j])

          ## skip everythig else  
          else:
            continue

        fout.write('\t' + str.join(':', new_v))
      fout.write('\n')



def vcf_view(args):
  global vcf_in, vcf_out
  vcf_out.write(str.join('\t', vcf_headers) + '\t')
  vcf_out.write(str.join('\t', vcf_samples) + '\n')

  for line in vcf_in:
    site_i = vcf_record(line)

    if args.merge_INDEL: 
      site_i.merge_indel()
      
    if args.seperate_allele: 
      site_i.write_multi(vcf_out, args)
    else:
      site_i.write(vcf_out)
    


def vcf_header_pass_through():
  global vcf_in, vcf_out
  global vcf_headers, vcf_samples
  for line in vcf_in:
    if re.search('^##', line):
      vcf_out.write(line)
    elif re.search('^#\w', line):
      vcf_out.write('##vcf_addon_command=' + str.join(' ', sys.argv) + '\n')
      lls = re.split('\t', line.rstrip())
      vcf_headers = lls[0:9]
      vcf_samples = lls[9:]
      return
    else:
      return


def vcf_input_open(vcf_in_file):
  global vcf_in
  if vcf_in_file is None:               vcf_in = sys.stdin
  elif re.search('\.gz$', vcf_in_file): vcf_in = gzip.open(vcf_in_file, 'rt', encoding='utf-8')
  else:                                 vcf_in =      open(vcf_in_file, 'r')


def vcf_output_open(vcf_out_file):
  global vcf_out
  if   vcf_out_file is None:             vcf_out = sys.stdout
  elif re.search('\.gz$', vcf_out_file): vcf_out = gzip.open(vcf_out_file, 'wt', encoding='utf-8')
  else:                                  vcf_out =      open(vcf_out_file, 'w')


def main():
  parser = argparse.ArgumentParser(description     = '''Additional tools for working with VCF files.
    Note, this is only tested with vcf files merged by GLnexus''')

  parser.add_argument('command',  help='''Sub commands: view''')
  parser.add_argument('input',    nargs='?',  help='''Input VCF file or gzipped VCF file, not required if input is from STDIN''')
  parser.add_argument('region',   nargs='?',  help='''Genomic region''')
  parser.add_argument('-o', '--output',    required=False,      help='''Output VCF filename, not required if output is STDOUT''')
  parser.add_argument('--seperate_allele', action='store_true', help='''If a site has both SNP and INDEL alleles, 
    seperate this site into a SNP site and an INDEL site''')
  parser.add_argument('--merge_INDEL',     action='store_true', help='''If a site has multiple INDEL alleles, merge them,
    and update GT,AF and other fields''')
  parser.add_argument('--merge_SNP',       action='store_true', help='''If a site has multiple SNP alleles, merge them,
    and update GT,AF and other fields''')
  parser.add_argument('--GLnexus_ID',      action='store_true', help='''Trust IDs in vcf files by GLNexus, given ID like chr1_12345_T_A 
    trust POS to be 12345, REF to be T, ALT to be A, print POS,REF,ALT coded in GLNexus ID with --seperate_allele function''')

  args = parser.parse_args()
  
  if args.command is None:
    print('No command')
    exit(1)

  vcf_input_open(args.input)
  vcf_output_open(args.output)

  if args.command == 'view': 
    vcf_header_pass_through()
    vcf_view(args)
 
  global vcf_in, vcf_out
  if not (args.input is None):  vcf_in.close()
  if not (args.output is None): vcf_out.close()

if __name__ == "__main__":
  main()

