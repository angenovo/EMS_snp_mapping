#!/bin/bash
# genetic mapping for mutants from EMS screen
# first the mutant vcf file is filtered with background vcf file
# then filtered with genotype and quality 
# if the mutant is sterial the heterozygous mutant could be sequence
# 
# usage: snpFA -b=/paht/to/background.vcf
#              -m=/path/to/mutant.vcf
#              -o=/path/to/outdir
# -bg|--background
# -m|--mutant
# -o|--outdir
#
# output files are filtered vcf file end with _rmbg_homo300.vcf
# SNP annotation file end with _ann_filter300.vcf
# out.log from snpEff
# snpEff.gene.txt from snpEff
# snpEff_summary.html from snpEff

#
# example:
# /home/ubuntu/test/code/snpFA.sh -bg=/home/ubuntu/test/data/vcffiles/1182F3_rmdupl.raw.vcf -m=/home/ubuntu/test/data/vcffiles/1182F1_rmdupl.raw.vcf -o=/home/ubuntu/test/data/vcffiles
#

## the snpEff v4 is not suppourt ce10 database, I installed the snpEff v3_6 
snpEff="/home/ubuntu/Software/snpEff_v3_6/snpEff/snpEff.jar"
snpSift="/home/ubuntu/Software/snpEff_v3_6/snpEff/SnpSift.jar"

for i in "$@"
do
case $i in 
  -bg=*|--background=*)
  BACKGROUND="${i#*=}"
  ;;
  -m=*|--mutant=*)
  MUTANT="${i#*=}"
  ;;
  -o=*|--outdir=*)
  OUTDIR="${i#*=}"
  ;;
  *)
  ;;
esac
done

## filter vcf files with the background line 1274F3

cd ${OUTDIR}
vcftools --vcf ${MUTANT} \
         --exclude-positions ${BACKGROUND} \
         --stdout --recode | java -jar $snpSift \
         filter "(GEN[*].GT = '1/1') & (QUAL > 300)" \
         > ${MUTANT}_rmbg_homo300.vcf
java -Xmx600M -jar $snpEff WS220.66 \
         -c /home/ubuntu/Software/snpEff_v3_6/snpEff/snpEff.config \
         ${MUTANT}_rmbg_homo300.vcf > ${MUTANT}_ann_filter300.vcf
cd





