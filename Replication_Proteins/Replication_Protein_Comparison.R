# Replication_Protein Genome Comparison

# The goal here is to view the location of the proteins genes

library(seqinr)
library(Biostrings)
library(genoPlotR)



# set the AAStringSet with names based on the order they appear in the OG csv file, 
## may have to reorder the OG data frames with match and order.
ecoli_strains <- c( "042","11368","APECO78","BW2952","CE10","E2348/69 ","EDL933",   
                    "IAI1","IAI39","MG1655","REL606","S88","Sakai DNA","SMS-3-5",  
                    "TW14359","UMN026","UTI89" 
                    
)

ecoli_strains

# set wd to the folder with the csv files in
setwd('/Users/dangoodall/Documents/PhD/termination_mac/BLAST_results/Ecoli/Replication_Proteins')

# read in protein csv from tblastn program
dnaA  <- read.csv('DnaA_Ecoli_tblastn.csv', 
                  header=FALSE)

dnaB  <- read.csv('DnaB_Ecoli_tblastn.csv', 
                  header=FALSE)

recA  <- read.csv('RecA_Ecoli_tblastn.csv', 
                  header=FALSE)

recB  <- read.csv('RecB_Ecoli_tblastn.csv', 
                  header=FALSE)

recC  <- read.csv('RecC_Ecoli_tblastn.csv', 
                  header=FALSE)

recD  <- read.csv('RecD_Ecoli_tblastn.csv', 
                  header=FALSE)

recG  <- read.csv('RecG_Ecoli_tblastn.csv', 
                  header=FALSE)

rnase  <- read.csv('RNaseHI_Ecoli_tblastn.csv', 
                   header=FALSE)

ruvA  <- read.csv('RuvA_Ecoli_tblastn.csv', 
                  header=FALSE)

ruvB  <- read.csv('RuvB_Ecoli_tblastn.csv', 
                  header=FALSE)

ruvC  <- read.csv('RuvC_Ecoli_tblastn.csv', 
                  header=FALSE)

xerC  <- read.csv('XerC_Ecoli_tblastn.csv', 
                  header=FALSE)

xerD  <- read.csv('XerD_Ecoli_tblastn.csv', 
                  header=FALSE)

# set wd back to the project folder
setwd("/Users/dangoodall/Documents/PhD/termination_mac/")

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

#keep only strain name, start and end
dnaA <- subset(dnaA, select = c(1:3),)
dnaB <- subset(dnaB, select = c(1:3),)
recA <- subset(recA, select = c(1:3),)
recB <- subset(recB, select = c(1:3),)
recC <- subset(recC, select = c(1:3),)
recD <- subset(recD, select = c(1:3),)
recG <- subset(recG, select = c(1:3),)
rnase <- subset(rnase, select = c(1:3),)
ruvA <- subset(ruvA, select = c(1:3),)
ruvB <- subset(ruvB, select = c(1:3),)
ruvC <- subset(ruvC, select = c(1:3),)
xerC <- subset(xerC, select = c(1:3),)
xerD <- subset(xerD, select = c(1:3),)


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

# fourth column which indicates the Protein name variable 
dnaA[4] <- 'dnaA'
dnaB[4] <- 'dnaB'
recA[4] <- 'recA'
recB[4] <- 'recB'
recC[4] <- 'recC'
recD[4] <- 'recD'
recG[4] <- 'recG'
rnase[4] <- 'rnase'
ruvA[4] <- 'ruvA'
ruvB[4] <- 'ruvB'
ruvC[4] <- 'ruvC'
xerC[4] <- 'xerC'
xerD[4] <- 'xerD'


# create concatenated df with all the proteins
total_df <- rbind(dnaA,
                  dnaB,
                  recA,
                  recB,
                  recC,
                  recD,
                  recG,
                  rnase,
                  ruvA,
                  ruvB,
                  ruvC,
                  xerC,
                  xerD)

rownames(total_df) <- NULL
total_df

# send to csv for easy reference laer on
write.csv(total_df, '/Users/dangoodall/Documents/PhD/termination_mac/BLAST_results/Ecoli/Replication_Proteins/Protein_Positions_Ecoli.csv')

# separate by strain name
O42     <- total_df[total_df$V1=='042',]
E11368  <- total_df[total_df$V1=='11368',]
APECO78 <- total_df[total_df$V1=='APECO78',]
BW2952  <- total_df[total_df$V1=='BW2952',]
CE10    <- total_df[total_df$V1=='CE10',]
E2348   <- total_df[total_df$V1=='E2348/69 ',]
EDL933  <- total_df[total_df$V1=='EDL933',]
IAI1    <- total_df[total_df$V1=='IAI1',]
IAI39   <- total_df[total_df$V1=='IAI39',]
MG1655  <- total_df[total_df$V1=='MG1655',]
REL606  <- total_df[total_df$V1=='REL606',]
S88     <- total_df[total_df$V1=='S88',]
Sakai   <- total_df[total_df$V1=='Sakai DNA',]
SMS35   <- total_df[total_df$V1=='SMS-3-5',]
TW14359 <- total_df[total_df$V1=='TW14359',]
UMN026  <- total_df[total_df$V1=='UMN026',]
UTI89   <- total_df[total_df$V1=='UTI89',]

# rename columns
colnames(O42) <- c('strain', 'start', 'end', 'name')
colnames(E11368) <- c('strain', 'start', 'end', 'name')
colnames(APECO78) <- c('strain', 'start', 'end', 'name')
colnames(BW2952) <- c('strain', 'start', 'end', 'name')
colnames(CE10) <- c('strain', 'start', 'end', 'name')
colnames(E2348) <- c('strain', 'start', 'end', 'name')
colnames(EDL933) <- c('strain', 'start', 'end', 'name')
colnames(IAI1) <- c('strain', 'start', 'end', 'name')
colnames(IAI39) <- c('strain', 'start', 'end', 'name')
colnames(MG1655) <- c('strain', 'start', 'end', 'name')
colnames(REL606) <- c('strain', 'start', 'end', 'name')
colnames(S88) <- c('strain', 'start', 'end', 'name')
colnames(Sakai) <- c('strain', 'start', 'end', 'name')
colnames(SMS35) <- c('strain', 'start', 'end', 'name')
colnames(TW14359) <- c('strain', 'start', 'end', 'name')
colnames(UMN026) <- c('strain', 'start', 'end', 'name')
colnames(UTI89) <- c('strain', 'start', 'end', 'name')

# reorder as name, start, end
order <- c('name', 'start', 'end', 'strain')
O42 <- O42[,order]
E11368 <- E11368[,order]
APECO78 <- APECO78[,order]
BW2952 <- BW2952[,order]
CE10 <- CE10[,order]
E2348 <- E2348[,order]
EDL933 <- EDL933[,order]
IAI1 <- IAI1[,order]
IAI39 <- IAI39[,order]
MG1655 <- MG1655[,order]
REL606 <- REL606[,order]
S88 <- S88[,order]
Sakai <- Sakai[,order]
SMS35 <- SMS35[,order]
TW14359 <- TW14359[,order]
UMN026 <- UMN026[,order]
UTI89 <- UTI89[,order]

# drop 4th column
O42 <- O42[-4]
E11368 <- E11368[-4]
APECO78 <- APECO78[-4]
BW2952 <- BW2952[-4]
CE10 <- CE10[-4]
E2348 <- E2348[-4]
EDL933 <- EDL933[-4]
IAI1 <- IAI1[-4]
IAI39 <- IAI39[-4]
MG1655 <- MG1655[-4]
REL606 <- REL606[-4]
S88 <- S88[-4]
Sakai <- Sakai[-4]
SMS35 <- SMS35[-4]
TW14359 <- TW14359[-4]
UMN026 <- UMN026[-4]
UTI89 <- UTI89[-4]

# add strand column with +
O42$strand <- '+'
E11368$strand <- '+'
APECO78$strand <- '+'
BW2952$strand <- '+'
CE10$strand <- '+'
E2348$strand <- '+'
EDL933$strand <- '+'
IAI1$strand <- '+'
IAI39$strand <- '+'
MG1655$strand <- '+'
REL606$strand <- '+'
S88$strand <- '+'
Sakai$strand <- '+'
SMS35$strand <- '+'
TW14359$strand <- '+'
UMN026$strand <- '+'
UTI89$strand <- '+'


# reset index
rownames(O42) <- NULL
rownames(E11368) <- NULL
rownames(APECO78) <- NULL
rownames(BW2952) <- NULL
rownames(CE10) <- NULL
rownames(E2348) <- NULL
rownames(EDL933) <- NULL
rownames(IAI1) <- NULL
rownames(IAI39) <- NULL
rownames(MG1655) <- NULL
rownames(REL606) <- NULL
rownames(S88) <- NULL
rownames(Sakai) <- NULL
rownames(SMS35) <- NULL
rownames(TW14359) <- NULL
rownames(UMN026) <- NULL
rownames(UTI89) <- NULL


# sanity
O42
E11368
APECO78
BW2952
CE10
E2348
EDL933
IAI1
IAI39
MG1655
REL606
S88
Sakai
SMS35
TW14359
UMN026
UTI89



####################################
#  DNA_SEG  OBJECTS WITH GENOPLOTR
###################################
colours <- c('blue','red', 'green', 'orange', 'darkred', 'darkgreen', 'darkblue', 'yellow', 'pink', 'purple', 'gray','black', 'blue')

# O42
O42_dnaseg <- dna_seg(O42)
O42_dnaseg$col <- colours
O42_dnaseg

#E11368
E11368_dnaseg <- dna_seg(E11368)
E11368_dnaseg$col <- colours
E11368_dnaseg

#APECO78
APECO78_dnaseg <- dna_seg(APECO78)
APECO78_dnaseg$col <- colours
APECO78_dnaseg

#BW2952
BW2952_dnaseg <- dna_seg(BW2952)
BW2952_dnaseg$col <- colours
BW2952_dnaseg

#CE10
CE10_dnaseg <- dna_seg(CE10)
CE10_dnaseg$col <- colours
CE10_dnaseg

#E2348
E2348_dnaseg <- dna_seg(E2348)
E2348_dnaseg$col <- colours
E2348_dnaseg

#EDL933
EDL933_dnaseg <- dna_seg(EDL933)
EDL933_dnaseg$col <- colours
EDL933_dnaseg

#IAI1
IAI1_dnaseg <- dna_seg(IAI1)
IAI1_dnaseg$col <- colours
IAI1_dnaseg

#IAI39
IAI39_dnaseg <- dna_seg(IAI39)
IAI39_dnaseg$col <- colours
IAI39_dnaseg

#MG1655
MG1655_dnaseg <- dna_seg(MG1655)
MG1655_dnaseg$col <- colours
MG1655_dnaseg

#REL606
REL606_dnaseg <- dna_seg(REL606)
REL606_dnaseg$col <- colours
REL606_dnaseg

#S88
S88_dnaseg <- dna_seg(S88)
S88_dnaseg$col <- colours
S88_dnaseg

#Sakai
Sakai_dnaseg <- dna_seg(Sakai)
Sakai_dnaseg$col <- colours
Sakai_dnaseg

#SMS35
SMS35_dnaseg <- dna_seg(SMS35)
SMS35_dnaseg$col <- colours
SMS35_dnaseg

#TW14359
TW14359_dnaseg <- dna_seg(TW14359)
TW14359_dnaseg$col <- colours
TW14359_dnaseg

#UMN026
UMN026_dnaseg <- dna_seg(UMN026)
UMN026_dnaseg$col <- colours
UMN026_dnaseg

#UTI89
UTI89_dnaseg <- dna_seg(UTI89)
UTI89_dnaseg$col <- colours
UTI89_dnaseg



#dna_seg list
# add names of the strain
new_list <- list(O42_dnaseg, 
                 E11368_dnaseg,
                 APECO78_dnaseg,
                 BW2952_dnaseg,
                 CE10_dnaseg,
                 E2348_dnaseg,
                 EDL933_dnaseg,
                 IAI1_dnaseg,
                 IAI39_dnaseg,
                 MG1655_dnaseg,
                 REL606_dnaseg,
                 S88_dnaseg,
                 Sakai_dnaseg,
                 SMS35_dnaseg,
                 TW14359_dnaseg,
                 UMN026_dnaseg,
                 UTI89_dnaseg)

names(new_list) <- ecoli_strains
length(new_list)

#mid positions of chromosome
midpos1 <- middle(new_list[[1]])
midpos2 <- middle(new_list[[2]])
midpos3 <- middle(new_list[[3]])
midpos4 <- middle(new_list[[4]])
midpos5 <- middle(new_list[[5]])
midpos6 <- middle(new_list[[6]])
midpos7 <- middle(new_list[[7]])
midpos8 <- middle(new_list[[8]])
midpos9 <- middle(new_list[[9]])
midpos10 <- middle(new_list[[10]])
midpos11 <- middle(new_list[[11]])
midpos12 <- middle(new_list[[12]])
midpos13 <- middle(new_list[[13]])
midpos14 <- middle(new_list[[14]])
midpos15 <- middle(new_list[[15]])
midpos16 <- middle(new_list[[16]])
midpos17 <- middle(new_list[[17]])

# annotations
annot1 <- annotation(x1 = midpos1, text = new_list[[1]]$name, rot = 90, col = 'black')
annot2 <- annotation(x1 = midpos2, text = new_list[[2]]$name, rot = 90, col = 'black')
annot3 <- annotation(x1 = midpos3, text = new_list[[3]]$name, rot = 90, col = 'black')
annot4 <- annotation(x1 = midpos4, text = new_list[[4]]$name, rot = 90, col = 'black')
annot5 <- annotation(x1 = midpos5, text = new_list[[5]]$name, rot = 90, col = 'black')
annot6 <- annotation(x1 = midpos6, text = new_list[[6]]$name, rot = 90, col = 'black')
annot7 <- annotation(x1 = midpos7, text = new_list[[7]]$name, rot = 90, col = 'black')
annot8 <- annotation(x1 = midpos8, text = new_list[[8]]$name, rot = 90, col = 'black')
annot9 <- annotation(x1 = midpos9, text = new_list[[9]]$name, rot = 90, col = 'black')
annot10 <- annotation(x1 = midpos10, text = new_list[[10]]$name, rot = 90, col = 'black')
annot11 <- annotation(x1 = midpos11, text = new_list[[11]]$name, rot = 90, col = 'black')
annot12 <- annotation(x1 = midpos12, text = new_list[[12]]$name, rot = 90, col = 'black')
annot13 <- annotation(x1 = midpos13, text = new_list[[13]]$name, rot = 90, col = 'black')
annot14 <- annotation(x1 = midpos14, text = new_list[[14]]$name, rot = 90, col = 'black')
annot15 <- annotation(x1 = midpos15, text = new_list[[15]]$name, rot = 90, col = 'black')
annot16 <- annotation(x1 = midpos16, text = new_list[[16]]$name, rot = 90, col = 'black')
annot17 <- annotation(x1 = midpos17, text = new_list[[17]]$name, rot = 90, col = 'black')

# annotations list
annots <- list(annot1, annot2, annot3, annot4, annot5, annot6, annot7, annot8, annot9, annot10,
               annot11, annot12, annot13, annot14, annot15, annot16, annot17
               )

# plot
plot_gene_map(dna_segs=new_list, comparisons=NULL,
              annotations = annots, annotation_height = 3, annotation_cex = 0.6,
              main = '13 Replication Proteins - Chromosome Comparison',
              dna_seg_scale=FALSE, gene_type = 'arrows', arrow_head_len	= 120000, dna_seg_label_cex=0.8)






