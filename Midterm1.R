#1. import and align your dna sequences
#first i have to set the working directory
setwd ("Users/kylecameron/Documents/GitHub/kyle-bioinformatics")
#load necessary packages
library(msa)
library(seqinr)
library(phangorn)
library(Biostrings)
library(UniprotR)
library(protti)

#download the new fasta file
mySequences <- Biostrings::readDNAStringSet("midtermfasta")
#i have to turn the file into an alignment
midtermAlignment <- msa::msa(mySequences)
#show alignment
print(midtermAlignment, show= "complete")

#2.Check to see how different...
#i'm going to add color to better detect the differences
#create a DNA string set from the alignment
seqhs<-DNAMultipleAlignment(midtermAlignment)
print(seqhs)

#now there is color to visualize the differences, 6 and 10 both look different
#they have deletion and substitutions.

#3. You suspect that an individual...
#i uploaded the fasta file onto genbank using the blast feature 
#The top match was "Homo sapiens hbb gene for beta globin, partial cds," 
#the accession number was GenBank LC121775.

#4. Find the individual that is the most different...
#the most different one is homo sapiens 6
#translate into a protein.
#define DNA sequence
dna_string <- DNAString("AATCTACTCCCAGGAGCAGGGAGGGCAGGAGCCAGGGCTGGGCATGAAAGTCAGGGCAGAGCCATCTATTGCTTACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACCATGGTGCACCTGACTCCTGTGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGACAGGTTTAAGGAGACCAATAGAAACTGGGCATGTGGAGACAGAGAAGACTCTTGGGTTTCTGATAGGCACTGACTCTCTCTGCCTATTGGTCTATTTTCCCACCCTTAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGGTGAGTCTATGGGACCCTTGATGTTTTCTTTCCCCTTCTTTTCTATGGTTAAGTTCATGTCATAGGAAGGGG")

#translate in AA sequence
homo6AA <- Biostrings::translate(dna_string)
print(homo6AA)

#define sequence header
header <- ">Homo_sapiens_6"

#define AA sequence
sequence <- as.character(homo6AA)

#write sequence header and sequence to a FASTA file
writeLines(c(header, sequence), "homo_sapiens_6.fasta")

#5 use a database to figure out...
#i used the BLAST tool to identify the match for the protein. It is "Hemoglobin subunit beta, HBB, Homo sapiens (human)."
#accession number is A0A0J9YWK4.

#6 using R or database what disease...
#HBB (hemoglobin subunit beta) is associated with Beta thalassemia (HBB/LCRB), Beta-zero thalassemia, Beta-plus thalassemia, Dominant-beta thalassemia, sickle cell disease, Methemoglobinemia, beta-globin type, Hemoglobin C disease
#clearly this is a blood related disease, my guess would be one tof the variants of thalassemia, most likely Beta (HBB/LCRB)

#7 what is the 3-dimensional structure of this protein?




