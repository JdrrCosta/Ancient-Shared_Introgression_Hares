#Mapping
cat ./0_clean/Ltim_paired | xargs -n 1 -P 2 sh -c 'bwa-mem2 mem -M -t 20 ~/Z_refs/GCA_033115175.1_mLepEur2.pri_genomic.fna \
./0_clean/$0_R1.fq.gz ./0_clean/$0_R2.fq.gz | samtools view -b --threads 20 | samtools sort --threads 20 -T $0 > ./2_bwa/$0.bam'

#Flagstats
ls ./2_bwa/*.bam | sed 's|./2_bwa/||g' | sed 's/.bam//g' | xargs -n 1 -P 9 sh -c 'samtools flagstat ./2_bwa/$0.bam >./3_flagstats/$0.stats'

#Remove duplicates and compile libraries for WGS (ie. merge all files from the same individual)
#diffrent codes for diffrent amount of bam files
cat Ltim_2_bam.txt | xargs -n 3 -P 1 sh -c 'picard -Xmx10g MarkDuplicates REMOVE_DUPLICATES=true \
ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 TMP_DIR=./4_remdup/$2.tmp \
INPUT=./2_bwa/$0.bam INPUT=./2_bwa/$1.bam OUTPUT=./4_remdup/$2.rmd.bam \
METRICS_FILE=./4_remdup/$2.stats.rmd.bam.metrics'

cat Ltim_1_bam.txt | xargs -n 2 -P 1 sh -c 'picard -Xmx10g MarkDuplicates REMOVE_DUPLICATES=true \
ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 TMP_DIR=./4_remdup/$1.tmp \
INPUT=./2_bwa/$0.bam OUTPUT=./4_remdup/$1.rmd.bam \
METRICS_FILE=./4_remdup/$1.stats.rmd.bam.metrics'

cat Ltim_8bam | xargs -n 9 -P 1 sh -c 'picard -Xmx10g MarkDuplicates REMOVE_DUPLICATES=true \
ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 TMP_DIR=./4_remdup/$8.tmp \
INPUT=./2_bwa/$0.bam INPUT=./2_bwa/$1.bam INPUT=./2_bwa/$2.bam INPUT=./2_bwa/$3.bam INPUT=./2_bwa/$4.bam \
INPUT=./2_bwa/$5.bam INPUT=./2_bwa/$6.bam INPUT=./2_bwa/$7.bam OUTPUT=./4_remdup/$8.rmd.bam \
METRICS_FILE=./4_remdup/$8.stats.rmd.bam.metrics'

#Run FlagStats on new bams
ls ./2_bwa/*.bam | sed 's|./2_bwa/||g' | sed 's/.bam//g' | xargs -n 1 -P 1 sh -c 'samtools flagstat ./2_bwa/$0.bam >./3_flagstats/$0.stats'

#Index all the new bam files
ls ./4_remdup/*.bam | xargs -n 1 -P 1 sh -c 'samtools index $0'

# In order to reduce the number of miscalls of INDELs in your data it is helpful to realign your raw gapped alignment with the Broad’s GATK Realigner.
#Adding read groups
ls ./4_remdup/*.rmd.bam | sed 's|./4_remdup/||' | sed 's/.rmd.bam//' | xargs -n 1 -P 3 sh -c 'picard AddOrReplaceReadGroups RGSM=$0 \
RGPL=ILLUMINA RGPU=$0 RGLB=$0 INPUT=./4_remdup/$0.rmd.bam OUTPUT=./4_remdup/$0.rmd.RG.bam'

#Index the bam files again
ls ./4_remdup/*.rmd.RG.bam | xargs -n 1 -P 6 sh -c 'samtools index $0'

#RealignerTargetCreator (Define intervals to target for local realignment)(using GATK3 since it is no longer available as of GATK4)
#Notes: The -Xmx4g flag is a java flag to specify the max allocated memory in this case 4Gb
ls ./4_remdup/*.rmd.RG.bam | sed 's|./4_remdup/||' | sed 's/.rmd.RG.bam//' | xargs -n 1 -P 1 sh -c 'gatk3 -Xmx4g -T RealignerTargetCreator \
-R ~/Z_refs/GCA_033115175.1_mLepEur2.pri_genomic.fna -I ./4_remdup/$0.rmd.RG.bam -o ./5_realign/$0.intervals'

## IndelRealigner (Perform local realignment of reads around indels)
ls ./4_remdup/*.rmd.RG.bam | sed 's|./4_remdup/||' | sed 's/.rmd.RG.bam//' | xargs -n 1 -P 6 sh -c 'gatk3 -Xmx4g -T IndelRealigner \
-R ~/Z_refs/GCA_033115175.1_mLepEur2.pri_genomic.fna -I ./4_remdup/$0.rmd.RG.bam -targetIntervals ./5_realign/$0.intervals -o ./5_realign/$0.realigned.bam'

#Index again
ls ./5_realign/*.realigned.bam | xargs -n 1 -P 6 sh -c 'samtools index $0'

#rename chr in bam to intuitive names
for i in {1..6}; do
    samtools view -H Ltim${i}.realigned.bam | \
    sed -e 's/SN:CM065551.1/SN:chr1/' \
        -e 's/SN:CM065552.1/SN:chr2/' \
        -e 's/SN:CM065553.1/SN:chr3/' \
        -e 's/SN:CM065554.1/SN:chr4/' \
        -e 's/SN:CM065555.1/SN:chr5/' \
        -e 's/SN:CM065556.1/SN:chr6/' \
        -e 's/SN:CM065557.1/SN:chr7/' \
        -e 's/SN:CM065558.1/SN:chr8/' \
        -e 's/SN:CM065559.1/SN:chr9/' \
        -e 's/SN:CM065560.1/SN:chr10/' \
        -e 's/SN:CM065561.1/SN:chr11/' \
        -e 's/SN:CM065562.1/SN:chr12/' \
        -e 's/SN:CM065563.1/SN:chr13/' \
        -e 's/SN:CM065564.1/SN:chr14/' \
        -e 's/SN:CM065565.1/SN:chr15/' \
        -e 's/SN:CM065566.1/SN:chr16/' \
        -e 's/SN:CM065567.1/SN:chr17/' \
        -e 's/SN:CM065568.1/SN:chr18/' \
        -e 's/SN:CM065569.1/SN:chr19/' \
        -e 's/SN:CM065570.1/SN:chr20/' \
        -e 's/SN:CM065571.1/SN:chr21/' \
        -e 's/SN:CM065572.1/SN:chr22/' \
        -e 's/SN:CM065573.1/SN:chr23/' \
        -e 's/SN:CM065574.1/SN:chrX/' \
        -e 's/SN:CM065575.1/SN:chrY/' \
        -e 's/SN:AJ421471.1/SN:chrM/' | \
    samtools reheader - Ltim${i}.realigned.bam > Ltim${i}.realigned_chr.bam
done
