
# now after assembly we can map reads and get a coverage estimate. 

samtools faidx phage/contigs.fasta NODE_1_length_141366_cov_30.585153 > phage_contig.fasta
bowtie2-build phage_contig.fasta phage_bowtie_db
bowtie2 -1 SRR6764339_1.fastq -2 SRR6764339_2.fastq -x phage_bowtie_db -S mapping_output.SAM
samtools view -bS -F 4 mapping_output.SAM > mapping_output.bam
samtools sort mapping_output.bam mapping_output.sort
samtools index mapping_output.sort.bam
samtools idxstats mapping_output.sort.bam
bedtools genomecov -d -g phage_contig.fasta -ibam mapping_output.sort.bam > phage_coverage.txt



# download the genome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/328/885/GCA_002328885.1_ASM232888v1/GCA_002328885.1_ASM232888v1_genomic.fna.gz
gunzip GCA_002328885.1_ASM232888v1_genomic.fna.gz

# now we need to use bowtie2. If we need to install first we can run:
sudo apt install bowtie2

# now we need to index the genome so we can map reads
bowtie2-build GCA_002328885.1_ASM232888v1_genomic.fna UBA2153


# get some sample reads
fastq-dump -X 1000000 --gzip --split-3 SRR5322088


# now we can map reads and create a SAM file
bowtie2 -1 SRR5322088_1.fastq -2 SRR5322088_2.fastq -x UBA2153 -S mapping_output.SAM

# running Samtools Version: 0.1.19-96b5f2294a
samtools view -bS -F 4 mapping_output.SAM > mapping_output.bam
samtools sort mapping_output.bam mapping_output.sort
samtools index mapping_output.sort.bam
samtools idxstats mapping_output.sort.bam


