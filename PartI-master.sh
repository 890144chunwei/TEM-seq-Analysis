#!/usr/bin/bash
#SBATCH --mem 100768
#SBATCH -p mhgcp
#SBATCH --output=/storage/goodell/projects/chunweic/slurm_out/220923_bsmap_%j.out
#SBATCH -e /storage/goodell/projects/chunweic/slurm_out/220923_bsmap_%j.err # Standard output and error log
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=chunweic@bcm.edu # Email to which notifications will be sent

pwd; hostname; date

FASTQDIR="/storage/goodell/projects/chunweic/220914_Chunwei_TEMseq_2/220914"
OUTDIR="/storage/goodell/projects/chunweic/220914_Chunwei_TEMseq_2/Align"
HOMEDIR="/storage/goodell/home/chunweic"


for FILE in $PROJECTDIR/*R1.fastq ;
do
  bsmap -a $PROJECTDIR/${FILE%_R1.fastq}_R1.fastq -b $PROJECTDIR/${FILE%_R1.fastq}_R2.fastq \
  -d $HOME/mm10/mm10_reference_genome.fa -o $FASTQDIR/$FILE.bam -s 16 -n 0 -r 2 -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 2 -z 64 -v 7

#bsmap -a $FASTQDIR/293T_WT_doxo_flag1_S10_L001_R1_001.fastq -b $FASTQDIR/293T_WT_doxo_flag1_S10_L001_R2_001.fastq \
#-d $HOMEDIR/hg19/GRCh37.primary_assembly.genome.fa -o $OUTDIR/WT_doxo_Flag1.bam -s 16 -n 1 -r 2 -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 2

done

samtools sort $OUTDIR/WT_doxo_Flag1.bam $OUTDIR/WT_doxo_Flag1_sort
samtools index $OUTDIR/WT_doxo_Flag1_sort.bam
samtools rmdup -S $OUTDIR/WT_doxo_Flag1_sort.bam $OUTDIR/WT_doxo_Flag1_drm.bam
samtools sort $OUTDIR/WT_doxo_Flag1_drm.bam $OUTDIR/WT_doxo_Flag1_drm_sort
samtools index $OUTDIR/WT_doxo_Flag1_drm_sort.bam
bedtools bamtobed -i $OUTDIR/WT_doxo_Flag1_drm_sort.bam > $OUTDIR/WT_doxo_Flag1.bed
bedtools intersect -v -a $OUTDIR/WT_doxo_Flag1.bed -b $HOMEDIR/hg19-blacklist.v2.bed > $OUTDIR/WT_doxo_Flag1_bl.bed
bamCoverage -b $OUTDIR/WT_doxo_Flag1_drm_sort.bam -bl $HOMEDIR/hg19-blacklist.v2.bed -o $OUTDIR/WT_doxo_Flag1_cov.bw --normalizeUsing RPKM

bamCoverage -b $TEMOUTPUT/${SAMPLE1}_bsmap_drm_sort.bam -bl $HOME/mm10-blacklist.v2.bed -o $TEMOUTPUT/${SAMPLE1}_bsmap_cov.bw --normalizeUsing RPKM #computes read coverage per bins or regions

macs3 callpeak -g mm -f BAMPE -t $OUTDIR/WT_NT_Flag1_drm_sort.bam $OUTDIR/WT_NT_Flag2_drm_sort.bam $OUTDIR/WT_NT_Flag3_drm_sort.bam \ 
-c $OUTDIR/WT_NT_mIgG_drm_sort.bam -q 0.05 -n $OUTDIR/WT_NT_Flag --nomodel --keep-dup=all --call-summits

cat $OUTDIR/WT_NT_Flag_peaks.bed $OUTDIR/WT_doxo_Flag_peaks.bed > $OUTDIR/WT_Flag_combined_peaks.bed
bedtools sort -i $OUTDIR/WT_Flag_combined_peaks.bed > $OUTDIR/WT_Flag_combined_peaks_sort.bed
bedtools merge -i $OUTDIR/WT_Flag_combined_peaks_sort.bed > $OUTDIR/WT_Flag_combined_peaks_merge.bed  

bedtools intersect -v -a $OUTDIR/WT_NT_Flag_peaks.bed -b $OUTDIR/WT_doxo_Flag_peaks.bed -wa > $OUTDIR/WT_Flag_a-b.bed
bedtools intersect -v -a $OUTDIR/WT_doxo_Flag_peaks.bed -b $OUTDIR/WT_NT_Flag_peaks.bed -wa > $OUTDIR/WT_Flag_b-a.bed

awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}' $OUTDIR/WT_Flag_combined_peaks_merge.bed > $OUTDIR/WT_Flag_combined_merge.saf

featureCounts -T 7 -p -F SAF -a $OUTDIR/WT_Flag_combined_merge.saf -o $OUTDIR/WT_NT_Flag1_fc.txt $OUTDIR/WT_NT_Flag1_drm_sort.bam
featureCounts -T 7 -p -F SAF -a $OUTDIR/WT_Flag_combined_merge.saf -o $OUTDIR/WT_NT_Flag2_fc.txt $OUTDIR/WT_NT_Flag2_drm_sort.bam

annotatePeaks.pl $OUTDIR/WT_NT_Flag_summits.bed mm10 -gtf $HOMEDIR/hg19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf > $OUTDIR/WT_NT_Flag.txt
cut -f 1-3 $OUTDIR/WT_NT_Flag_peaks.narrowPeak > $OUTDIR/WT_NT_Flag_peaks.bed

annotatePeaks.pl $MACSOUTPUT/${SAMPLE}_summits.bed mm10 -gtf /storage/goodell/home/u244406/macs3/mm10/gencode.vM10.chr_patch_hapl_scaff.annotation.gtf > $MACSOUTPUT/${SAMPLE}.txt
cut -f 1-3 $MACSOUTPUT/${SAMPLE}_peaks.narrowPeak > $MACSOUTPUT/${SAMPLE}_peaks.bed

