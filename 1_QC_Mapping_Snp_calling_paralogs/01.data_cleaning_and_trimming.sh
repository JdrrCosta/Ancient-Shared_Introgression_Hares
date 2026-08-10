#Data Cleaning and Trimming
#Raw data files list
ls * > 0_Ltim_samples

#Run fastqc
ls ./0_Raw/*/*.gz | xargs -n 1 -P 10 sh -c 'fastqc $0 -o ./1_fastqc/'

#Run trimmomatic
cat 0_Ltimidus_samples | xargs -n 4 -P 10 sh -c 'trimmomatic PE $0 $1 /data/evochange/jmarques/B_Ltim/0_clean/$2.fq.gz /data/evochange/jmarques/B_Ltim/0_clean/$2_unpaired.fq.gz /data/evochange/jmarques/B_Ltim/0_clean/$3.fq.gz /data/evochange/jmarques/B_Ltim/0_clean/$3_unpaired.fq.gz ILLUMINACLIP:/data/evochange/jmarques/Z_REFs/TruSeq3-PE-2.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36'

#Run Fastqc on the trimmed reads
ls ./0_clean/*.gz | xargs -n 1 -P 10 sh -c 'fastqc $0 -o ./1_fastqc/'

