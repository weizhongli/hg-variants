#!/bin/bash

for sample in $(cut -f 1 NGS-samples); do
  if [ ! -e "$sample/hg-call/HG.sv.pass.tsv" ]; then
    bcftools query -i 'FILTER="PASS"' -f '%CHROM\t%POS\t%ID\t%END\t%REF\t%ALT\t%TYPE\t%INFO/END\t%INFO/SVTYPE\t%INFO/MATEID\t%INFO/SVLEN\t[%GT\t]\n' $sample/hg-call/HG.sv.vcf.gz > $sample/hg-call/HG.sv.pass.tsv
  fi
done

