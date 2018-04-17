
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


