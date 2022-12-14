
##################################################################################################################
#	usage  STRgenotype.hipstr.sh -i $mycramlist -b $mystrs -d $mydir -p $myprefix -r $myreadlength 
##################################################################################################################

#!/bin/bash

while getopts i:b:d:p:r: flag
do
    case "${flag}" in
        i) mycramlist=${OPTARG};;
        b) mystrs=${OPTARG};;
        d) mydir=${OPTARG};;
		p) myprefix=${OPTARG};;
		r) myreadlength=${OPTARG};;
    esac
done


echo "Genotyping TRs with HipSTR"
echo "WD: $mydir"													# provide full path for working directory
echo "Input: $mycramlist"									# provide file with full path for cram files (if multiple samples will be genotyped, please include individual sample per row)
echo "STRs: $mystrs"        							# provide bed file with regions of interest (columns as follows chr,start,end,motif length, copies in the reference genome,STRId,motif)
echo "Prefix: $myprefix"									# prefix (eg. batch01)
echo "Read Length: $myreadlength bp"			# provide max length of STRs to be genotyped (smaller than sequencing read length).

ref=<reference genome> 										# provide full path for reference genome in fasta format
echo "ReferenceGenome: $ref"

mkdir -p ${mydir}/hipstr/
mkdir -p ${mydir}/hipstr/raw/

echo "... Runnning HipSTR for $mycramlist ...\n"

cd ${mydir}
HipSTR \
--bam-files ${mycramlist} \
--fasta ${ref} \
--regions ${mystrs} \
--str-vcf ${mydir}/hipstr/raw/${myprefix}.${mycramlist}.HipSTR.vcf.gz \
--log ${mydir}/hipstr/raw/${myprefix}.${mycramlist}.HipSTR.out \
--max-str-len ${myreadlength}


echo "The end"
