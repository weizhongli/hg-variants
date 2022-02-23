#!/bin/bash

WORKPWD=$(pwd)
AWSS3=s3://bucket/folder/GLnexus
CHRLST="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM"

################################################################################
#### gVCF merging and joint variant calling

mkdir GLnexus.input
cd $WORKPWD/GLnexus.input

#### step 0, rename gvcf files into $sample.g.vcf.gz, required by GLnexus 
#### the first column of file samples contain the sample name
for sample in $(cut -f 1 samples); do
    echo $sample 
    if [ ! -e $sample.g.vcf.gz ] ; then 
        echo "process $sample.g.vcf.gz"
        ln -s ../$sample/HG.hard-filtered.gvcf.gz     $sample.g.vcf.gz
        ln -s ../$sample/HG.hard-filtered.gvcf.gz.tbi $sample.g.vcf.gz.tbi
    fi 
done

#### step 1, split a whole gvcf file into 25 gvcf files for each chromsome, done in parallel
batchsize=15
for CHR in $CHRLST; do
    ind=0
    for sample in $(cut -f 1 samples); do
        if [ ! -e $sample.$CHR.g.vcf.gz ] ; then
            echo "process $sample.$CHR.g.vcf.gz"
            bcftools view -o $sample.$CHR.g.vcf.gz -O z $sample.g.vcf.gz $CHR &
            ind=$(expr $ind + 1)
            if [ $(expr $ind % $batchsize) == 0 ] ; then
                echo waiting threads    
                wait
            fi
        fi
    done
    wait
    ind=0
    for sample in $(cut -f 1 samples); do
        if [ ! -e $sample.$CHR.g.vcf.gz.tbi ] ; then
            echo "process tabix $sample.$CHR.g.vcf.gz"
            tabix $sample.$CHR.g.vcf.gz &
            ind=$(expr $ind + 1)
            if [ $(expr $ind % $batchsize) == 0 ] ; then
                echo waiting threads    
                wait
            fi
        fi
    done
    wait
done
#### now $sample.g.vcf.gz $sample.g.vcf.gz.tbi and source HG.hard-filtered.gvcf.gz HG.hard-filtered.gvcf.gz.tbi
#### can be deleted


cd $WORKPWD/GLnexus.input
#### step 2, gVCF merging and joint variant calling using GLnexus
CHRLST="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22"
for CHR in $CHRLST; do
    if [ ! -e GLnexus.$CHR.bcf]; then
        #### GLnexus.input contains all the $sample.$CHR.g.vcf.gz 
        ls GLnexus.input/*$CHR.g.vcf.gz > GLnexus.$CHR.input.list

        #### with --mem-gbytes 64 --threads 8 will actually need 256GB RAM
        /vol1/data/apps/glnexus/glnexus_cli --config DeepVariant --mem-gbytes 64 --threads 8 \
            --dir GLnexus.$CHR.DB --list GLnexus.$CHR.input.list > GLnexus.$CHR.bcf 2>GLnexus.$CHR.bcf.log 

        rm GLnexus.$CHR.input.list
        rm GLnexus.$CHR.DB/*sst
        aws s3 cp GLnexus.$CHR.bcf     $AWSS3/GLnexus.$CHR.bcf
        aws s3 cp GLnexus.$CHR.bcf.log $AWSS3/GLnexus.$CHR.bcf.log
    fi 
done

#### step 3, gVCF merging and joint variant calling using GLnexus for chrX, chrY and chrM
echo -e 'chrX\t1\t2781479'           >  chrX.PARs.bed
echo -e 'chrX\t155701383\t156030895' >> chrX.PARs.bed
echo -e 'chrX\t2781480\t155701382'   >  chrX.nonPARs.bed
echo -e 'chrY\t1\t2781479'           >  chrY.PARs.bed
echo -e 'chrY\t56887903\t57217415'   >> chrY.PARs.bed
echo -e 'chrY\t2781480\t56887902'    >  chrY.nonPARs.bed

CHR=chrX
ls GLnexus.input/*$CHR.g.vcf.gz > GLnexus.$CHR.input.list
/vol1/data/apps/glnexus/glnexus_cli --config DeepVariant --mem-gbytes 16 --threads 4 --bed $CHR.PARs.bed \
  --dir GLnexus.$CHR.PARs.DB --list GLnexus.$CHR.input.list > GLnexus.$CHR.PARs.bcf 2>GLnexus.$CHR.PARs.bcf.log 
rm GLnexus.$CHR.input.list
rm GLnexus.$CHR.PARs.DB/*sst
aws s3 cp GLnexus.$CHR.PARs.bcf     $AWSS3/GLnexus.$CHR.PARs.bcf
aws s3 cp GLnexus.$CHR.PARs.bcf.log $AWSS3/GLnexus.$CHR.PARs.bcf.log
#### no need to run chrY.PARs

#### GLnexus.input/samples.XX is a subset of samples with XX type
CHR=chrX
SEX=F
rm GLnexus.$CHR.$SEX.input.list
for sample in $(cat GLnexus.input/samples.XX); do echo "GLnexus.input/$sample.$CHR.g.vcf.gz" >> GLnexus.$CHR.$SEX.input.list; done
/vol1/data/apps/glnexus/glnexus_cli --config DeepVariant --mem-gbytes 24 --threads 6 --bed $CHR.nonPARs.bed \
  --dir GLnexus.$CHR.$SEX.nonPARs.DB --list GLnexus.$CHR.$SEX.input.list > GLnexus.$CHR.$SEX.nonPARs.bcf 2>GLnexus.$CHR.$SEX.nonPARs.bcf.log
rm GLnexus.$CHR.$SEX.input.list
rm GLnexus.$CHR.$SEX.nonPARs.DB/*sst
aws s3 cp GLnexus.$CHR.$SEX.nonPARs.bcf     $AWSS3/GLnexus.$CHR.$SEX.nonPARs.bcf
aws s3 cp GLnexus.$CHR.$SEX.nonPARs.bcf.log $AWSS3/GLnexus.$CHR.$SEX.nonPARs.bcf.log



################################################################################
#### filter rare variants

CHRLST="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX.PARs chrX.F.nonPARs"
## need chrX.M.nonPARs chrY.M.nonPARs chrM

#### MAF > 0.05
mkdir MAF.05.vcf
for CHR in $CHRLST; do
   if [ ! -e MAF.05.vcf/MAF.05.$CHR.vcf.gz ]; then
       bcftools view -i 'AF>0.05' -o MAF.05.vcf/MAF.05.$CHR.vcf.gz -O z --threads 2 GLnexus.$CHR.bcf
   fi
done
aws s3 sync MAF.05.vcf $AWSS3/MAF.05.vcf


################################################################################
#### convert vcf to tsv
#### these are not really useful, there are other tools for analysis
mkdir MAF.05.vcf.tsv
for CHR in $CHRLST; do
    ## GT
    if [[ (! -e MAF.05.vcf.tsv/MAF.05.$CHR.GT.tsv) && (! -e MAF.05.vcf.tsv/MAF.05.$CHR.GT.tsv.gz) ]]; then
        echo -en 'CHR\tPOS\tID\tEND\tREF\tALT\tTYPE\tAF\tAQ' > MAF.05.vcf.tsv/MAF.05.$CHR.GT.tsv
        bcftools query -l MAF.05.vcf/MAF.05.$CHR.vcf.gz | perl -e 'while(<>){chop; print "\t$_";} print "\n";' >> MAF.05.vcf.tsv/MAF.05.$CHR.GT.tsv
        bcftools query -f'%CHROM\t%POS\t%ID\t%END\t%REF\t%ALT\t%TYPE\t%INFO/AF\t%INFO/AQ\t[%GT\t]\n' MAF.05.vcf/MAF.05.$CHR.vcf.gz >> MAF.05.vcf.tsv/MAF.05.$CHR.GT.tsv
        gzip MAF.05.vcf.tsv/MAF.05.$CHR.GT.tsv
    fi

    ## depth
    if [[ (! -e MAF.05.vcf.tsv/MAF.05.$CHR.DP.tsv) && (! -e MAF.05.vcf.tsv/MAF.05.$CHR.DP.tsv.gz) ]]; then
        echo -en 'ID' > MAF.05.vcf.tsv/MAF.05.$CHR.DP.tsv
        bcftools query -l MAF.05.vcf/MAF.05.$CHR.vcf.gz | perl -e 'while(<>){chop; print "\t$_";} print "\n";' >> MAF.05.vcf.tsv/MAF.05.$CHR.DP.tsv
        bcftools query -f'%ID\t[%DP\t]\n' MAF.05.vcf/MAF.05.$CHR.vcf.gz >> MAF.05.vcf.tsv/MAF.05.$CHR.DP.tsv
        gzip MAF.05.vcf.tsv/MAF.05.$CHR.DP.tsv
    fi
done
aws s3 sync MAF.05.vcf.tsv $AWSS3/MAF.05.vcf.tsv






