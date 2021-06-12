Bulk ATAC-seq analysis of human T cells (naive and activated, obtained from blood)
The command line analysis is depicted exemplary.


#QC, cutadapt
fastqc *fastq.gz

cutadapt -m 13 --trim-n -g AATGATACGGCGACCACCGAGATCTACAC -a GTGTAGATCTCGGTGGTCGCCGTATCATT -G CAAGCAGAAGACGGCATACGAGAT -A ATCTCGTATGCCGTCTTCTGCTTG -o s301_TN601_R1_cutadapt.fastq.gz -p s301_TN601_R2_cutadapt.fastq.gz samp301_S301_LOUT_R2_001.fastq.gz samp301_S301_LOUT_R2_001.fastq.gz

fastqc *cutadapt.fastq.gz

#alignment with bowtie 2
wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gzip -dk Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
bowtie2-build Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38

bowtie2 -p 12 --very-sensitive -k 10 -x /corgi/lisabeth/ATAC050521/ref/GRCh38 -1 ../s301_TN601_R1_cutadapt.fastq.gz -2 ../s301_TN601_R2_cutadapt.fastq.gz -S s301_TN601_cutadapt_bt2.sam

samtools view -S -b s301_TN601_cutadapt_bt2.sam > s301_TN601_cutadapt_bt2.bam
samtools sort s301_TN601_cutadapt_bt2.bam -o s301_TN601_cutadapt_bt2_sorted.bam

#PCR dedup
picard MarkDuplicates -I s301_TN601_cutadapt_bt2_sorted.bam -O s301_TN601_cutadapt_bt2_sorted_dedup.bam -M s301_TN601_cutadapt_bt2_sorted_dedupMetrics.txt -REMOVE_DUPLICATES TRUE

#QC
conda deactivate
java -jar /usr/local/bin/picard.jar BuildBamIndex -INPUT s301_TN601_cutadapt_bt2_sorted_dedup.bam
java -jar /usr/local/bin/picard.jar CollectInsertSizeMetrics -I s301_TN601_cutadapt_bt2_sorted_dedup.bam -O s301_TN601_cutadapt_bt2_sorted_dedup_SizeMetrics.txt -H s301_TN601_cutadapt_bt2_sorted_dedup_SizeMetrics.pdf
conda activate

#pileup
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
cut -f1,2 Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai > sizes.genome

bedtools genomecov -bg -ibam s301_TN601_cutadapt_bt2_sorted_dedup.bam > s301_TN601_cutadapt_bt2_sorted_dedup.bedgraph
bedGraphToBigWig s301_TN601_cutadapt_bt2_sorted_dedup.bedgraph /corgi/lisabeth/ATAC050521/ref/sizes.genome pileup301_TN601.bigwig

cp *bigwig /data/public_http/public/henriksson/temp/lis/final

#technical repl
samtools merge TN30_cutadapt_bt2_sorted_dedup.bam s302_TN301_cutadapt_bt2_sorted_dedup.bam s306_TN302_cutadapt_bt2_sorted_dedup.bam

#peak calling with homer
bedtools bamtobed -i TN30_cutadapt_bt2_sorted_dedup.bam > TN30_cutadapt_bt2_sorted_dedup.bed
makeTagDirectory TN30/ TN30_cutadapt_bt2_sorted_dedup.bed -format bed
findPeaks TN30/ -style factor -o auto

#python prep
bedtools bamtobed -i s301_TN601_cutadapt_bt2_sorted.bam > s301_TN601_cutadapt_bt2_sorted.bed

wget http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.chr.gtf.gz
gzip -dk Homo_sapiens.GRCh38.104.chr.gtf.gz

#differential ana
cut -f2,3,4 peaksTN60.txt > peaksTN30.bed
tail -n +36 peaksTN30.bed > peaksTN30x.bed
bedtools sort -i peaksTN30x.bed > peaksTN30xs.bed
bedtools merge -i peaksTA30xs.bed -i peaksTN30xs.bed > peaksT30.bed
bedtools coverage -a peaksT30.bed -b TN30_cutadapt_bt2_sorted_dedup.bed > peaksT30_TN30.bed

#MT contamination
samtools view -F 0x40 s301_TN601_cutadapt_bt2.bam | cut -f 1 | uniq | wc -l
samtools view -F 0x40 s301_TN601_cutadapt_bt2.bam | grep MT | cut -f 1 | uniq | wc -l
