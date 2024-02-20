#Load Packages that i might need
library(msa)
library(seqinr)
library(Biostrings)
library(dplyr)
library(tidyverse)
library(tidyr)
library(UniprotR)
library(protti)

#set the working directorty
setwd("/Users/kylecameron/Documents/GitHub/kylebioinformatics")

#load the files into r and align them and print them to view it all.
seqsmt1 <- readDNAStringSet("midtermk.fasta")
seqsmt1aligned <- msa(seqsmt1)
print(seqsmt1aligned, show="complete")

#now i'm going to see if there are any differences between them
#since I'm in biostrings I'm going to create a seqinr file
allhumans <- readDNAStringSet("midtermk.fasta")
allhumansaligned <- msa(allhumans)
allhumansaligned
homosapiensal <- msaConvert(allhumansaligned, type = "seqinr::alignment")
homosapiensal
DistanceMatrix <- dist.alignment(homosapiensal, "identity")
DistanceMatrix

#from here i see that homo sapiens 6 is the odd one out not only is it the least similar it also has deletions and missmatches.
#I added the file to blast and found Homo sapiens hbb gene for beta globin to be the closest match.
#accession number GenBank LC121775

# homo_sapiens_6 is the most different
#now I need to turn it into a protein seq
string6 <- DNAString("AATCTACTCCCAGGAGCAGGGAGGGCAGGAGCCAGGGCTGGGCATGAAAGTCAGGGCAGAGCCATCTATTGCTTACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACCATGGTGCACCTGACTCCTGTGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGACAGGTTTAAGGAGACCAATAGAAACTGGGCATGTGGAGACAGAGAAGACTCTTGGGTTTCTGATAGGCACTGACTCTCTCTGCCTATTGGTCTATTTTCCCACCCTTAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGGTGAGTCTATGGGACCCTTGATGTTTTCTTTCCCCTTCTTTTCTATGGTTAAGTTCATGTCATAGGAAGGGG")
hs6AA <- Biostrings::translate(string6)
print(hs6AA)

# now i need to turn it into a fasta file, define header and and write it all to a fasta
header <- ">Homo_sapiens_6"
sequence <- as.character(hs6AA)
writeLines(c(header, sequence), "homo_sapiens_6.fasta")

#using blast my protein matches Hemoglobin subunit beta, HBB, Homo sapiens (A0A0J9YWK4)
#HBB is associated with some form of Beta thalassemia
#I do believe they have some form of Beta thalassemia


