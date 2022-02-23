#!/bin/bash

#### follow genotype.sh

WORKPWD=$(pwd)
CHRLST="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX.PARs"

cd $WORKPWD/MAF.05.vcf
#### merge multiple indels at the same loci to the most abundant indel - by vcf_addons.py --merge_INDEL
#### seperate ALT alleles at the same loci to seperate records - by vcf_addons.py --seperate_allele
#### fix POS, REF, ALT data due to changes made by GLnexus during merging - by vcf_addons.py --GLnexus_ID
#### filter out variants with call rate < 0.95 - by bcftools view -i 'F_MISSING < 0.05'
#### sort vcf files
for CHR in $CHRLST; do
  if [ ! -e "MAF.05.$CHR.long.vcf.gz" ] ; then
    /vol1/data/apps/hg-variants/workflow/vcf_addons.py view MAF.05.$CHR.vcf.gz --merge_INDEL --seperate_allele --GLnexus_ID \
      | bcftools view -i 'F_MISSING < 0.10' -O z -o MAF.05.$CHR.long-unsorted.vcf.gz
    bcftools sort -m 8G -o MAF.05.$CHR.long.vcf.gz -O z -T . MAF.05.$CHR.long-unsorted.vcf.gz
    rm MAF.05.$CHR.long-unsorted.vcf.gz
  fi
done

#### concat VCF files from each chr to a single one
if [ -e auto.list ]; then
  rm auto.list
fi
for CHR in $CHRLST; do
  echo "MAF.05.$CHR.long.vcf.gz" >> auto.list
done
bcftools concat -o MAF.05.auto.vcf.gz -O z -f auto.list
bcftools index MAF.05.auto.vcf.gz


#### ANNOVAR annotation
#### first generate a vcf file with only 1 sample to save space and time
bcftools view -H MAF.05.auto.vcf.gz | cut -f 1-10 > MAF.05.auto.1s.vcf
table_annovar.pl MAF.05.auto.1s.vcf /vol1/data/apps/annovar/humandb38/ \
  -buildver hg38 -out MAF.05.auto-ann -remove \
  -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation g,r,f,f,f -nastring . -vcfinput -polish
rm MAF.05.auto.1s.vcf

#### Plink2 analysis
cd $WORKPWD
if [ ! -e "MAF.05.plink2" ] ; then
  mkdir MAF.05.plink2
fi

cd $WORKPWD/MAF.05.plink2
#### need to prepare MAF.05.plink2/samples.psam, with columns: #IID	SEX
#### need to prepare MAF.05.plink2/plink2.pheno.tsv, with columns: #IID	meta_data1	metadata2	metadata3(e.g. race)	group
if [[ -e "samples.psam" && -e "plink2.pheno.tsv" ]] ; then
  plink2 --vcf ../MAF.05.vcf/MAF.05.auto.vcf.gz --psam samples.psam --vcf-half-call h  --out MAF.05.auto.merged.long
  plink2 --pfile MAF.05.auto.merged.long --psam samples.psam --pheno iid-only plink2.pheno.tsv --pheno-name group --glm --out MAF.05.auto.merged.long-glm --adjust
else
  echo "file samples.psam or plink2.pheno.tsv  missing"
fi

#### MAF calc by Plink2
# 1=control/2=case
grep "2$" plink2.pheno.tsv | cut -f 1 > plink2.pheno.group2.IIDs
grep "1$" plink2.pheno.tsv | cut -f 1 > plink2.pheno.group1.IIDs
plink2 --pfile MAF.05.auto.merged.long --psam samples.psam --keep plink2.pheno.group1.IIDs --freq --out MAF.05.auto.merged.long.group1.freq
plink2 --pfile MAF.05.auto.merged.long --psam samples.psam --keep plink2.pheno.group2.IIDs --freq --out MAF.05.auto.merged.long.group2.freq

rlist="black white"
for RACE in $rlist; do
  grep $RACE plink2.pheno.tsv | cut -f 1 > plink2.pheno.$RACE.IIDs
  plink2 --pfile MAF.05.auto.merged.long --psam samples.psam --keep plink2.pheno.$RACE.IIDs --freq --out MAF.05.auto.merged.long.$RACE.freq
done


