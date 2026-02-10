# Genome_Analysis_Pipeline_bash
Highlight of genomic analysis procedures performed in bash


**1.Demultiplex Illumina reads (Cutadapt)**

#Trim leading bases and demultiplex paired-end reads using Golay barcodes.
cutadapt -u 4 -U 4 \
  -o illumina_multiplexed_trimmed4_R1.fastq \
  -p illumina_multiplexed_trimmed4_R2.fastq \
  illumina_multiplexed_R1.fastq illumina_multiplexed_R2.fastq

  #Demultiplex by barcode pair:
  cutadapt -e 0 --discard-untrimmed --cores=4 --no-indels \
  -g file:golay_barcodes_F.fasta \
  -G file:golay_barcodes_R.fasta \
  -o "{name1}-{name2}.1.fastq" \
  -p "{name1}-{name2}.2.fastq" \
  illumina_multiplexed_trimmed4_R1.fastq illumina_multiplexed_trimmed4_R2.fastq

  **2. Rename demultiplexed samples using mapping file**

  #Mapping file:
#Golay_L1-Golay_R5	sampleAB
#Golay_L2-Golay_R7	sampleBJ
#Golay_L3-Golay_R3	sampleJK
#Golay_L4-Golay_R8	sampleTT
#Golay_L5-Golay_R9	sampleXZ

while IFS=$'\t' read -r barcode_combo sample_id
do
  mv "${barcode_combo}.1.fastq" "${sample_id}_R1.fastq"
  mv "${barcode_combo}.2.fastq" "${sample_id}_R2.fastq"
done < barcode_sample_identifiers.tsv

rm -f Golay*.fastq



**3. Generate FASTQ statistics (SeqKit)**

#Creating combined stats table
echo -e "filename\tnumber_of_seqs\tnumber_of_bases\tavg_length\tN50\tGC_content" > seq_stats.tsv

for i in sample*_R*.fastq
do
  seqkit stats -a -T "$i" | sed '1d' | awk '{print $1"\t"$4"\t"$5"\t"$7"\t"$13"\t"$18}' >> seq_stats.tsv
done



**Download Fructobacillus genomes programmatically**

#Accession list:
GCF_963580165.1
GCF_963579665.1
GCF_963580705.1
GCF_963580715.1
GCF_963580745.1
GCF_963580755.1
GCF_963580805.1



#Download loop:
for i in $(cat fructo_accessions.txt)
do
  datasets download genome accession $i --filename ${i}.zip
done



**5. Calculate ANI similarity (FastANI)**

ls fructo_genomes_only/*.fna | fastANI -t 4 \
  -r fructo_genomes_only/GCF_963580165.1_LMG_32999_genomic.fna \
  --ql /dev/stdin \
  -o fructo_fastani_results.tsv

  #Fornat output table:
  echo -e "genome_query\tgenome_reference\tANI" > fructo_fastani_summary.tsv

awk '{print $1"\t"$2"\t"$3}' fructo_fastani_results.tsv \
  | sed 's/\.fna//g' \
  | sed 's/_genomic//g' \
  | sed 's/fructo_genomes_only\///g' \
  >> fructo_fastani_summary.tsv


**6. ANI visualization**
ANIclustermap \
  -i fructo_genomes_only \
  -o fructobacillus_heatmap \
  -fig_width 12 \
  -cmap_gamma 0.5

**7. Restriction enzyme digestion simulation (VSEARCH)**

#Download genome:
datasets download genome accession GCF_001077535.1 --filename clostr_genome.zip
unzip clostr_genome.zip

#Digestion script:
bash RE_digest.sh genome.fna RE_cuts



**8. Restriction enzyme summary statistics**
bash RE_digest_summary.sh genome.fna RE_cuts


#Example output:
enzyme  GCF_001077535

EcoRI   78

BamHI   15

PstI    258

AluI    17691


