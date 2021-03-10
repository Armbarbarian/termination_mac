# find Chi sites along genome to find where RecBCD could recognise
# if Chi sites are near ter sites then this would clarify that RecBCD would be able to act near ter

library(Biostrings)
library(seqinr)

# Chi site that blocks RecBCD from digesting too far
chi <- 'CTGGATGCTGGTGGGATAAC'
chi_r <- reverseComplement(DNAString(chi))
chi_core <- 'GCTGGTGG'
chi_core_r <- reverseComplement(DNAString(chi_core))

# load genome MG1655 to test amount of Chi seqs
MG <- read.fasta('MG1655.fasta', seqonly = TRUE)
class(MG)

# turn into DNAString
MGdna <- DNAString(MG[[1]])

# match chi with MG1655 to see how many sequences we get
match_chi <- matchPattern(chi, MGdna, max.mismatch = 3)
match_chi
# found 1 match at:  1680921 (0 mismatches)
# terB MG1655 is at: 1684226
# found 1 match at:  630664
# terI MG1655 is at: 625411
# are these too far away or could RecBCD be cleaving 4000nt long until it reaches Chi?


# reverse chi to find on - strand
match_chi_r <- matchPattern(toString(chi_r), MGdna, max.mismatch = 4)
match_chi_r
# no real matches here

# match just the core sequence to get more possibilities
match_chicore <- matchPattern(chi_core, MGdna)
match_chicore[,1]
# 499 possible Chi sequences
# clean to find if any of these are within the inner termination area

# reverse matching to find - strand chi seqs
match_chicore_r <- matchPattern(toString(chi_core_r), MGdna)
match_chicore_r
# 509 possible Chi sites with the inner core


###############################
#   CHI IN TER AREA
###############################
min <- 1090000
max <- 1690000
chi_df <- as.data.frame(match_chicore[,1])
chi_df

