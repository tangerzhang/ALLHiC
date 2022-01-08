#!/bin/bash

usage()
{
	echo "    Usage: `basename $0` -r reference -1 R1.fq -2 R2.fq -k group_count [-e enzyme] [-t threads] [-b bin_size]"
	echo "          -r: reference genome"
	echo "          -1: Lib_R1.fq.gz"
	echo "          -2: Lib_R2.fq.gz"
	echo "          -k: group_count"
	echo "          -e: enzyme_sites (HindIII: AAGCTT; MboI: GATC), default: HindIII"
	echo "          -t: threads, default: 10"
	echo "          -b: bin_size for hic heatmap, can be divided with comma, default: 500k"
	exit 0
}

### get options
while getopts ':r:1:2:k:e:t:b:' OPT; do
	case $OPT in
		r)
			ref="$OPTARG";;
		1)
			R1="$OPTARG";;
		2)
			R2="$OPTARG";;
		e)
			enzyme="$OPTARG";;
		k)
			group_count="$OPTARG";;
		t)
			threads="$OPTARG";;
		b)
			bin_size="$OPTARG";;
		?)
			usage;;
	esac
done
bwa="bwa"

### check required variants
if [ -z $ref ] || [ -z $R1 ] || [ -z $R2 ] || [ -z $group_count ]; then
	usage
fi

### set default values while optional variants were not set
if [ -z $threads ]; then
	threads=10
fi

if [ -z $bin_size ]; then
	bin_size=500k
fi

if [ -z $enzyme ]; then
	enzyme=AAGCTT
fi

enzyme=`echo $enzyme | tr '[a-z]' '[A-Z]'`

if [ $enzyme = HINDIII ]; then
	enzyme=AAGCTT
fi

if [ $enzyme = MBOI ]; then
	enzyme=GATC
fi

### link required files
ln -s ${ref} ./seq.fasta
ln -s ${R1} ./Lib_R1.fastq.gz
ln -s ${R2} ./Lib_R2.fastq.gz

### index reference genome
bwa index seq.fasta
samtools faidx seq.fasta


### 1st round of mapping
bwa mem -SP5M -t $threads seq.fasta Lib_R1.fastq.gz Lib_R2.fastq.gz \
     | samtools view -hF 256 - \
     | samtools sort -@ $threads -o sorted.bam -T tmp.ali
samtools index sorted.bam

### correct contig
ALLHiC_corrector -m sorted.bam -r seq.fasta -o seq.HiCcorrected.fasta -t $threads

### 2nd round of mapping
bwa index seq.HiCcorrected.fasta
samtools faidx seq.HiCcorrected.fasta
bwa mem -SP5M -t $threads seq.HiCcorrected.fasta Lib_R1.fastq.gz Lib_R2.fastq.gz \
     | samtools view -hF 256 - \
     | samtools sort -@ $threads -o sample.bwa_mem.bam -T tmp.ali


### filter bam
samtools view -bq 40 sample.bwa_mem.bam  |samtools view -bt seq.HiCcorrected.fasta.fai > sample.unique.bam
PreprocessSAMs.pl sample.unique.bam seq.HiCcorrected.fasta $enzyme

### partition
ALLHiC_partition -r seq.HiCcorrected.fasta -e $enzyme -k $group_count -b sample.unique.REduced.paired_only.bam

### optimize
rm cmd.list
for((K=1;K<=$group_count;K++));do echo "allhic optimize sample.unique.REduced.paired_only.counts_${enzyme}.${group_count}g${K}.txt sample.unique.REduced.paired_only.clm" >> cmd.list;done
ParaFly -c cmd.list -CPU $threads

### build
ALLHiC_build seq.HiCcorrected.fasta

### plot
samtools faidx groups.asm.fasta
cut -f1,2 groups.asm.fasta.fai|grep sample > chrn.list
ALLHiC_plot sample.bwa_mem.bam groups.agp chrn.list $bin_size pdf


