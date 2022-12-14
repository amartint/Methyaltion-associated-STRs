#!/bin/bash


while getopts d:i: flag
do
    case "${flag}" in
        d) directory=${OPTARG};;
        i) input=${OPTARG};;
    esac
done
echo "Directory: $directory"; echo "Input: $input";


echo "Extractintg genotypes from $input.filtered.vcf.gz"
cd $directory
bcftools query -f '[%CHROM\t%START\t%END\t%POS\t%ID\t%PERIOD\t%AN\t%SAMPLE\t%GT\t%REF\t%ALT\n]' $input.filtered.vcf.gz | awk '$9!="./."'| tr "|" "\t" |\
	awk 'BEGIN{ FS = "\t"; OFS = "\t"}{
        str = $11","$12;
        split(str,arr,",");
        $9=length(arr[$9+1])/$6;
        $10=length(arr[$10+1])/$6;
        print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' > $input.filtered.txt
        

echo "The end"

