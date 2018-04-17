
# Today we will be predicting genes from prokaryotic genomes

# first step, install Prodigal, a useful tool for gene prediction
sudo apt install prodigal

# Prodigal will predict genes from chromosomes (or contigs), translate those genes into amino acids, and produce annotation summary files such as gff, depending on what options you use. 
prodigal -i med4_genome.fna -a med4.proteins.faa -d med4.genes.fna -f gff -o med4.prodigal.gff

# or use GenBank output file for a summary
prodigal -i med4_genome.fna -a med4.proteins.faa -d med4.genes.fna -f gbk -o med4.prodigal.gbk

# Prodigal is only useful for predicting protein coding genes. What other kind of genes are there in genomes?

# Barrnap is useful for predicting rRNA genes
barrnap med4_genome.fna > med4.rRNA.gff
 
# unfortunately barrnap only provides the summary files (in this case gff). So we need to do a bit more legwork to get the actual sequences
bedtools getfasta -fi med4_genome.fna -bed med4.rRNA.gff -fo med4.rRNA.fasta

# 16S genes are extremely useful for classification. If you ever have a genome and you don't know what it is, a good first step is to identify any 16S ribosomal genes in the chromosome and use them for classification. 
