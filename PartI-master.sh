#!/usr/bin/bash
#SBATCH --mem 32768
#SBATCH -p short
#SBATCH --output=/storage/goodell/home/u244406/slurm_out/221213_ma_Mut_ckit_H2AZ_%j.out
#SBATCH -e /storage/goodell/home/u244406/slurm_out/221213_ma_Mut_ckit_H2AZ_%j.err # Standard output and error log
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=u244406@bcm.edu # Email to which notifications will be sent
pwd; hostname; date


BCLDIR="/storage/goodell/bcl/220901_Chunwei_TEMseq/Files"
GENOMDIR="/storage/goodell/home/chunweic/mm10/STARgenome"
FASTQDIR="/storage/goodell/fastq/220901_Chunwei_TEMseq"
PROJECTDIR="/storage/goodell/projects/chunweic/220901_Chunwei_TEMseq"
SAMPLE="Mut_ckit_H2AZ"
SAMPLE1="Mut_ckit_H2AZ2_S22"
SAMPLE2="Mut_ckit_H2AZ1_S21"
CONTROL="Mut_ckit_mIgG_S20" #igG


for FILE in $PROJECTDIR/*R1.fastq ;
do
  bsmap -a $PROJECTDIR/${FILE%_R1.fastq}_R1.fastq -b $PROJECTDIR/${FILE%_R1.fastq}_R2.fastq \
  -d $HOME/mm10/mm10_reference_genome.fa -o $FASTQDIR/$FILE.bam -s 16 -n 0 -r 2 -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 2 -z 64 -v 7
done

samtools sort $TEMOUTPUT/${SAMPLE1}_bsmap.bam -o $TEMOUTPUT/${SAMPLE1}_bsmap_sort.bam #sort aligment
samtools index $TEMOUTPUT/${SAMPLE1}_bsmap_sort.bam # index aligment
samtools rmdup -S $TEMOUTPUT/${SAMPLE1}_bsmap_sort.bam $TEMOUTPUT/${SAMPLE1}_bsmap_drm.bam #remove PCR duplicate
samtools sort $TEMOUTPUT/${SAMPLE1}_bsmap_drm.bam -o $TEMOUTPUT/${SAMPLE1}_bsmap_drm_sort.bam #sort aligment
samtools index $TEMOUTPUT/${SAMPLE1}_bsmap_drm_sort.bam #index aligment

bedtools bamtobed -i $TEMOUTPUT/${SAMPLE1}_bsmap_drm_sort.bam > $TEMOUTPUT/${SAMPLE1}_bsmap.bed # .bam to .bed convertion
bedtools intersect -v -a $TEMOUTPUT/${SAMPLE1}_bsmap.bed -b $HOME/mm10-blacklist.v2.bed > $TEMOUTPUT/${SAMPLE1}_bsmap_bl.bed #Find overlapping intervals in various ways

bamCoverage -b $TEMOUTPUT/${SAMPLE1}_bsmap_drm_sort.bam -bl $HOME/mm10-blacklist.v2.bed -o $TEMOUTPUT/${SAMPLE1}_bsmap_cov.bw --normalizeUsing RPKM #computes read coverage per bins or regions

macs3 callpeak \
-g mm -f BAMPE \
-t $TEMOUTPUT/${SAMPLE1}_bsmap_drm_sort.bam $TEMOUTPUT/${SAMPLE2}_bsmap_drm_sort.bam \
-c $TEMOUTPUT/${CONTROL}_bsmap_drm_sort.bam \
-q 0.05 --outdir $MACSOUTPUT \
-n ${SAMPLE} --nomodel --keep-dup=all --call-summits

annotatePeaks.pl $MACSOUTPUT/${SAMPLE}_summits.bed mm10 -gtf /storage/goodell/home/u244406/macs3/mm10/gencode.vM10.chr_patch_hapl_scaff.annotation.gtf > $MACSOUTPUT/${SAMPLE}.txt
cut -f 1-3 $MACSOUTPUT/${SAMPLE}_peaks.narrowPeak > $MACSOUTPUT/${SAMPLE}_peaks.bed

