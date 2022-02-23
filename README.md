# hg-variants

This is a GWAS pipeline for whole genome analysis starting from FASTQ files to find significant variants. The pipeline was tested on AWS cloud. It has the following major steps: 
  1. Prepare fastq files
  2. Run variant calls using the Illumina Dragen Bio-IT platform
  3. Join variants (gvcf) using GLnexus
  4. GWAS analysis using PLINK2
  5. Variant annotation using ANNOVAR

## Software requirement
  1. Python 3
  2. Illumina Dragen Bio-IT platform - https://www.illumina.com/products/by-type/informatics-products/dragen-bio-it-platform.html
  3. bcftools - https://samtools.github.io/bcftools/bcftools.html
  4. GLnexus - https://github.com/dnanexus-rnd/GLnexus
  5. PLINK2 - https://www.cog-genomics.org/plink/2.0/
  6. ANNOVAR - https://annovar.openbioinformatics.org/en/latest/
  7. ngomicswf - https://github.com/weizhongli/ngomicswf, only the NG-Omics-WF.py3 is needed

## Prepare fastq files
The pipeline starts with demultiplexed fastq files. These files should be stored in a AWS s3 folder. Besides the fastq files, The Stats.json produced during the demultiplexing process is also needed. The Stats.json includes information about sample names, number of reads per lane, etc. The FC_stat.py can be used to parse these json files and generate sample stats in tsv format (example/aws_folder/flowcell_all.tsv)

<pre>
FC_stat.py -f list_of_Stats.json -o output_prefix

Here, list_of_Stats.json is a text file with the s3 locations of these Stats.json, one per line
s3://bucket_name/folder1/folder1_1/RUN_ID_XYZ1/Stats/Stats.json
s3://bucket_name/folder1/folder1_1/RUN_ID_XYZ2/Stats/Stats.json
s3://bucket_name/folder1/folder1_1/RUN_ID_XYZ3/Stats/Stats.json
s3://bucket_name/folder1/folder1_1/RUN_ID_XYZ4/Stats/Stats.json
...
</pre>

## Run variant calls using the Illumina Dragen Bio-IT platform
This is done using the Dragen Bio-IT platform. An AWS FPGA instance (e.g. f1.4xlarge) must be used. The human genome hg38 (e.g. GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna) can be used as the reference genome. 

First the reference need to be indexed:
<pre>
cd 
mkdir GRch38cnv
dragen --build-hash-table true --ht-reference GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
    --enable-cnv true --output-directory GRch38cnv --ht-num-threads 32
</pre>

Second the reference need to be loaded:
<pre>
dragen -l -r GRch38cnv
</pre>

The next step is run the Dragen pipeline for the samples. With a large number of samples (e.g. 1000 samples), a workflow engine is needed. Here, the the NG-Omics-WF.py from https://github.com/weizhongli/ngomicswf is used. We can run these samples by batch (e.g. 100 samples/batch). To run:

<pre>
nohup /ephemeral/data/apps/ngomicswf/NG-Omics-WF.py3 -s NGS-samples -i dragen-variant-call-workflow.py -j vcf-upload -f &

Here, NGS-samples is a text file list the sample names, one per line. the sample name is the first column in example/aws_folder/flowcell_all.tsv
</pre>

Here, the file dragen-variant-call-workflow.py describes the steps and command lines to pull fastq files from AWS s3, run dragen and save results back to AWS s3. dragen-variant-call-workflow.py need to be edited to specify the parameters, e.g. file path, s3 folder, dragen license key etc. NG-Omics-WF.py3 will generate needed SH scritps. Some scripts (workflow/dragen-examples/*sh) are provided for demo purpose. 

After this step, the gVCF for SNVs and small indels, VCF for CNVs and VCF for SVs, other files and LOG files will be saved to AWS s3 folder.

