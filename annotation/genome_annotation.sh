
####################################################################
###### Part I : Gene Prediction using Prodigal and Barrnap #########
####################################################################

# Today we will be predicting both protein-coding and rRNA genes from prokaryotic genomes
# first step, install Prodigal, a useful tool for gene prediction
sudo apt install prodigal

# Prodigal will predict genes from chromosomes (or contigs), translate those genes into amino acids, and produce annotation summary files (such as "gene feature format", or gff, files), depending on what options you use. 
prodigal -i med4_genome.fna -a med4.proteins.faa -d med4.genes.fna -f gff -o med4.prodigal.gff

# or use GenBank output file for a summary
prodigal -i med4_genome.fna -a med4.proteins.faa -d med4.genes.fna -f gbk -o med4.prodigal.gbk

# Usually I am mainly interested in the amino acid sequences of the predicted genes, so in this case I would only really use the med4.proteins.faa file. 

# Prodigal is only useful for predicting protein coding genes. To predict rRNA genes such as 16S rRNA we need to use barrnap
# First, to install:
sudo apt install barrnap

# And then run:
barrnap med4_genome.fna > med4.rRNA.gff
 
# unfortunately barrnap only provides the summary files (in this case gff). So we need to do a bit more work to get the actual sequences
# In this case I will use a tool called bedtools to extract the nucleic acid sequences from particular regions of the chromosome, using the gff file from barrnap as a guide. 
# If you need to install bedtools try "sudo apt install bedtools". 
bedtools getfasta -fi med4_genome.fna -bed med4.rRNA.gff -fo med4.rRNA.fasta

# 16S genes are extremely useful for classification. If you ever have a genome and you don't know what it is, a good first step is to identify any 16S ribosomal genes in the chromosome and use them for classification. 

####################################################################
###### Part II : Gene Annotation using Hidden Markov Models ########
####################################################################

# Today we will be using a tool called HMMER3 to annotate genes in a genome. 
# We will be discussing Hidden Markov Models (HMMs). The HMMs that we go over today will be present in the file cog_hmms.hmm
# The descriptions of these models can be found in cog_descriptions.txt

# HMMs are models that allow us to compare a protein to a profile of other proteins rather than just individual proteins. 
# This is useful if we want to know if a protein belongs to a certain protein family. For example, we may wish to know if a protein is an RNA polymerase. We don't necessarily care which RNA polymerse it is most similar to, we just want to know if the protein belongs to that general family. This is useful for functional annotations. 


# The first step is to predict genes from the genome. We did this last week. In this case we will once again start with the Prochlorococcus MED4 genome, and predict genes using Prodigal. 
prodigal -i med4_genome.fna -a med4.proteins.faa -f gff -o med4.genes.gff 

# How many genes were found? Use seqtk to find out. 
seqtk comp med4.proteins.faa | wc

# You should get 1921 

# Now that we have the amino acid sequences of the predicted genes (in the .faa file) we can compare these genes to a set of Hidden Markove Models using HMMER3. 
# First step, install HMMER:
sudo apt install hmmer

# Now run the HMMER command. Note that the last two arguments are "positional arguments" since they do not have flags in front of them. The .hmm file and the query protein file must always be provided at the end, in that order. 
hmmsearch --tblout med4.hmmout -o med4.output cog_hmms.hmm med4.proteins.faa

# The main tabulated output we want is in med4.hmmout. Unfortunately the authors of HMMER made the output space-delimited, so it's a bit hard to look at or put in an Excel spreadsheet. 
# I made a small Python script that will parse through this output, pull out the best hit for each query protein, and put it in a tab-delimited output.  
python parse_hmmout.py med4.hmmout > med4.hmmout.parsed

# Now we should have the best hit for each protein, followed by the HMM that it hit to and the score. 

# Now you may notice some hits that have very low bit scores. This is because we did not use any quality cutoffs when we ran HMMER. Just like with BLAST, there is an e-value cutoff option that we can use. 
# For that we can use the following command:
hmmsearch -E 1e-10 --tblout med4.hmmout -o med4.output cog_hmms.hmm med4.proteins.faa

# Now let's practice again with another genome and see what we get. A new genome is in the file N_maritimus.fna. This is an Archaea called Nitrosopumilis maritimus, an abundant ammonia-oxidizing microbe in the ocean. 
# Here is the overall workflow:
prodigal -i N_maritimus.fna -a N_maritimus.faa -f gff -o N_maritimus.gff
hmmsearch -E 1e-10 --tblout N_maritimus.hmmout -o N_maritimus.output cog_hmms.hmm N_maritimus.faa
python parse_hmmout.py N_maritimus.hmmout > N_maritimus.hmmout.parsed

# what functional genes are present here that are not present in Prochlorococcus? Which genes are present in both?
