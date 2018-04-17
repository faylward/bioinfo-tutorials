
# today we will be working with the raw reads from the Staphylococcus phage 812 genome sequencing project.

# First we need to install a package manager called Miniconda. Go to the website here and download the 64-bit Linux distribution:
https://conda.io/miniconda.html
# then, once you have located this file on your computer, type 
bash Miniconda2-latest-Linux-x86_64.sh
# and follow the installation instructions. When asked to append the PATH info say "yes".
# After this you will need to close your terminal and re-open it before beginning again. Then you should be able to install tools using conda. For the module today we will need to install two tools: the sra-toolkit and an assembler called spades. 
conda install -c bioconda sra-tools
conda install -c bioconda spades

# To get the raw reads we need to use a program called the sra-toolkit. This is a tool maintained by NCBI to allow users to download data easily from the command line. 
# The toolkit allows users to specify the unique accession number for a project and then download the associated reads. 
# Here we will also use a command -X, which specifies how many reads we want to download. Since these datasets can be quite large we want to start with a small number. Here we will use 10,000
# we also want to use the --split-3 flag, which for Illumina data makes sure the forward and reverse reads are split into separate files. 
fastq-dump -X 10000 --split-3 SRR6764339

# Now once we have the raw reads we can begin assembling them using a program calls SPAdes. 
# Now actually running the assembly is fairly easy- we just specify the read files, an output folder name, the number of threads we want to use, and the k-mer length. 
spades.py -1 SRR6764339_1.fastq -2 SRR6764339_2.fastq -o phage -t 4 -k 21 &> log.txt



# now after assembly we can map reads and get a coverage estimate. 
# first we can get the longest contig and index it so we can map reads. 
samtools faidx phage/contigs.fasta NODE_1_length_141366_cov_30.585153 > phage_contig.fasta
bowtie2-build phage_contig.fasta phage_bowtie_db

# now we can map reads against the reference using bowtie2
bowtie2 -1 SRR6764339_1.fastq -2 SRR6764339_2.fastq -x phage_bowtie_db -S mapping_output.SAM

# after that we need to process the .sam file by converting it into a .bam file, and then sorting and indexing it before we can get read mapping statistics. 
samtools view -bS -F 4 mapping_output.SAM > mapping_output.bam
samtools sort mapping_output.bam mapping_output.sort
samtools index mapping_output.sort.bam
samtools idxstats mapping_output.sort.bam

# Now that we have our indexed .bam file, we can use bedtools to get a coverage estimate. I will pipe the output and cut it to retain only the third column, which will contain the coverage estimates for each bp of the genome. 
bedtools genomecov -d -g phage_contig.fasta -ibam mapping_output.sort.bam | cut -f 3 > phage_coverage.txt

