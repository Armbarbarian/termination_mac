# Replication protein alignments for Ecoli 
# from tblastn
# needs to be in tabular output of csv, it will be too difficult to see the sequences in raw pairwise

#libraries
library(seqinr)
library(Biostrings)
library(msa)
library(DECIPHER)

# set the AAStringSet with names based on the order they appear in the OG csv file, 
## may have to reorder the OG data frames with match and order.
ecoli_strains <- c( "042","11368","APECO78","BW2952","CE10","E2348/69 ","EDL933",   
                    "IAI1","IAI39","MG1655","REL606","S88","Sakai DNA","SMS-3-5",  
                    "TW14359","UMN026","UTI89" 
                   
)

ecoli_strains

# read in protein csv from tblastn program
dnaA  <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\BLAST_results\\Ecoli\\Replication_Proteins\\DnaA_Ecoli_tblastn.csv', 
                      header=FALSE)

dnaB  <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\BLAST_results\\Ecoli\\Replication_Proteins\\DnaB_Ecoli_tblastn.csv', 
                  header=FALSE)

recA  <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\BLAST_results\\Ecoli\\Replication_Proteins\\RecA_Ecoli_tblastn.csv', 
                  header=FALSE)

recB  <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\BLAST_results\\Ecoli\\Replication_Proteins\\RecB_Ecoli_tblastn.csv', 
                  header=FALSE)

recC  <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\BLAST_results\\Ecoli\\Replication_Proteins\\RecC_Ecoli_tblastn.csv', 
                  header=FALSE)

recD  <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\BLAST_results\\Ecoli\\Replication_Proteins\\RecD_Ecoli_tblastn.csv', 
                  header=FALSE)

recG  <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\BLAST_results\\Ecoli\\Replication_Proteins\\RecG_Ecoli_tblastn.csv', 
                  header=FALSE)

rnase  <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\BLAST_results\\Ecoli\\Replication_Proteins\\RNaseHI_Ecoli_tblastn.csv', 
                  header=FALSE)

ruvA  <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\BLAST_results\\Ecoli\\Replication_Proteins\\RuvA_Ecoli_tblastn.csv', 
                  header=FALSE)

ruvB  <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\BLAST_results\\Ecoli\\Replication_Proteins\\RuvB_Ecoli_tblastn.csv', 
                  header=FALSE)

ruvC  <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\BLAST_results\\Ecoli\\Replication_Proteins\\RuvC_Ecoli_tblastn.csv', 
                  header=FALSE)

xerC  <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\BLAST_results\\Ecoli\\Replication_Proteins\\XerC_Ecoli_tblastn.csv', 
                  header=FALSE)

xerD  <- read.csv('C:\\Users\\Danie\\Documents\\R\\termination\\BLAST_results\\Ecoli\\Replication_Proteins\\XerD_Ecoli_tblastn.csv', 
                  header=FALSE)

# MATCH V1 AGAINST ECOLI_STRAINS SO THEY ARE IN THE SAME ORDER
dnaA <- dnaA[order(dnaA$V1),]
dnaB <- dnaB[order(dnaB$V1),]
recA <- recA[order(recA$V1),]
recB <- recB[order(recB$V1),]
recC <- recC[order(recC$V1),]
recD <- recD[order(recD$V1),]
recG <- recG[order(recG$V1),]
rnase <- rnase[order(rnase$V1),]
ruvA <- ruvA[order(ruvA$V1),]
ruvB <- ruvB[order(ruvB$V1),]
ruvC <- ruvC[order(ruvC$V1),]
xerC <- xerC[order(xerC$V1),]
xerD <- xerD[order(xerD$V1),]

dnaA$V1


# turn the protein sequences into DNAStrings
dnaA <- AAStringSet(dnaA$V6)
dnaB <- AAStringSet(dnaB$V6)
recA <- AAStringSet(recA$V6)
recB <- AAStringSet(recB$V6)
recC <- AAStringSet(recC$V6)
recD <- AAStringSet(recD$V6)
recG <- AAStringSet(recG$V6)
rnase <- AAStringSet(rnase$V6)
ruvA <- AAStringSet(ruvA$V6)
ruvB <- AAStringSet(ruvB$V6)
ruvC <- AAStringSet(ruvC$V6)
xerC <- AAStringSet(xerC$V6)
xerD <- AAStringSet(xerD$V6)

# give names
names(dnaA) <- ecoli_strains
names(dnaB) <- ecoli_strains
names(recA) <- ecoli_strains
names(recB) <- ecoli_strains
names(recC) <- ecoli_strains
names(recD) <- ecoli_strains
names(recG) <- ecoli_strains
names(rnase) <- ecoli_strains
names(ruvA) <- ecoli_strains
names(ruvB) <- ecoli_strains
names(ruvC) <- ecoli_strains
names(xerC) <- ecoli_strains
names(xerD) <- ecoli_strains


#sanity
dnaA
dnaB
recA
recB
recC
recD
recG
rnase
ruvA
ruvB
ruvC
xerC
xerD


# align dnaA
dnaA_align <- AlignSeqs(dnaA)
dnaA_align
dnaA_msa <-  msa(dnaA_align, method = 'Muscle', order = 'input')

msaPrettyPrint(x=dnaA_msa, file="DnaA_Alignment.pdf",
               shadingMode="identical",
               shadingColors="grays",
               logoColors="chemical",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE)


# align dnaB
dnaB_align <- AlignSeqs(dnaB)
dnaB_align
dnaB_msa <-  msa(dnaB_align, method = 'Muscle', order = 'input')

msaPrettyPrint(x=dnaB_msa, file="dnaB_Alignment.pdf",
               shadingMode="identical",
               shadingColors="grays",
               logoColors="chemical",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE)

# align RecA
recA_align <- AlignSeqs(recA)
recA_align
recA_msa <-  msa(recA_align, method = 'Muscle', order = 'input')

msaPrettyPrint(x=recA_msa, file="recA_Alignment.pdf",
               shadingMode="identical",
               shadingColors="grays",
               logoColors="chemical",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE)

# align RecB
recB_align <- AlignSeqs(recB)
recB_align
recB_msa <-  msa(recB_align, method = 'Muscle', order = 'input')

msaPrettyPrint(x=recB_msa, file="recB_Alignment.pdf",
               shadingMode="identical",
               shadingColors="grays",
               logoColors="chemical",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE)

# align RecC
recC_align <- AlignSeqs(recC)
recC_align
recC_msa <-  msa(recC_align, method = 'Muscle', order = 'input')

msaPrettyPrint(x=recC_msa, file="recC_Alignment.pdf",
               shadingMode="identical",
               shadingColors="grays",
               logoColors="chemical",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE)

# align RecD
recD_align <- AlignSeqs(recD)
recD_align
recD_msa <-  msa(recD_align, method = 'Muscle', order = 'input')

msaPrettyPrint(x=recD_msa, file="recD_Alignment.pdf",
               shadingMode="identical",
               shadingColors="grays",
               logoColors="chemical",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE)

# align RecG
recG_align <- AlignSeqs(recG)
recG_align
recG_msa <-  msa(recG_align, method = 'Muscle', order = 'input')

msaPrettyPrint(x=recG_msa, file="recG_Alignment.pdf",
               shadingMode="identical",
               shadingColors="grays",
               logoColors="chemical",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE)


# align RNaseHI
rnase_align <- AlignSeqs(rnase)
rnase_align
rnase_msa <-  msa(rnase_align, method = 'Muscle', order = 'input')

msaPrettyPrint(x=rnase_msa, file="rnase_Alignment.pdf",
               shadingMode="identical",
               shadingColors="grays",
               logoColors="chemical",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE)

# align RuvA
ruvA_align <- AlignSeqs(ruvA)
ruvA_align
ruvA_msa <-  msa(ruvA_align, method = 'Muscle', order = 'input')

msaPrettyPrint(x=ruvA_msa, file="ruvA_Alignment.pdf",
               shadingMode="identical",
               shadingColors="grays",
               logoColors="chemical",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE)

# align RuvB
ruvB_align <- AlignSeqs(ruvB)
ruvB_align
ruvB_msa <-  msa(ruvB_align, method = 'Muscle', order = 'input')

msaPrettyPrint(x=ruvB_msa, file="ruvB_Alignment.pdf",
               shadingMode="identical",
               shadingColors="grays",
               logoColors="chemical",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE)

# align RuvC
ruvC_align <- AlignSeqs(ruvC)
ruvC_align
ruvC_msa <-  msa(ruvC_align, method = 'Muscle', order = 'input')

msaPrettyPrint(x=ruvC_msa, file="ruvC_Alignment.pdf",
               shadingMode="identical",
               shadingColors="grays",
               logoColors="chemical",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE)

# align XerC
xerC_align <- AlignSeqs(xerC)
xerC_align
xerC_msa <-  msa(xerC_align, method = 'Muscle', order = 'input')

msaPrettyPrint(x=xerC_msa, file="xerC_Alignment.pdf",
               shadingMode="identical",
               shadingColors="grays",
               logoColors="chemical",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE)


# align XerD
xerD_align <- AlignSeqs(xerD)
xerD_align
xerD_msa <-  msa(xerD_align, method = 'Muscle', order = 'input')

msaPrettyPrint(x=xerD_msa, file="xerD_Alignment.pdf",
               shadingMode="identical",
               shadingColors="grays",
               logoColors="chemical",
               showLogoScale="right",
               showConsensus=c("bottom"),
               consensusColors=c("ColdHot"),
               showLegend=FALSE,
               askForOverwrite=FALSE, verbose = TRUE)




