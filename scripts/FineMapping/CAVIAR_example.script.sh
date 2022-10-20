

#!/bin/bash


while getopts d:i: flag
do
    case "${flag}" in
        d) directory=${OPTARG};;
        i) input=${OPTARG};;
    esac
done
echo "Directory: $directory"; echo "Input: $input";

#myinput.pearson.txt. Matrix containing Pearson's correlation obtained from association between genotypes of mSTRs (average allelic size) and SNVs (dosage alternate allele) 
#myinput.zscore.txt. File containing two columns: geneticvariantID and zscore obtained from association between genotypes and DNA methlaytion values of target CpG

CAVIAR -o $input.caviar_c3.txt -l $input.pearson.txt -z $input.zscores.txt -r 0.95 -g 0.01 -c 3
