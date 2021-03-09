# GenoPlotR troubleshoot
# issue from reviewers is that the figures are hard to see
# trying to use the example of Hughey, dewy and lewy to enlarge the ter sites
## possibly try to add the arrows to show the direction.

library(genoPlotR)
library(Biostrings)
library(seqinr)

## simple example
## dna segments
## data.frame with several genes
names1 <- c("terA", "terB", "terC")
starts1 <- c(1341745,
             1684226,
             1609157)
ends1 <- c(starts1 + 10000)
strands1 <- c("+", -1, -1)
cols1 <- c("blue", "grey", "red")


df1 <- data.frame(name=names1, start=starts1, end=ends1,
                  strand=strands1, col='black')
df1

dna_seg1 <- dna_seg(df1)

fill <- c('red', 'blue', 'orange')
dna_seg1$fill <- fill

is.dna_seg(dna_seg1)


## with only one gene, or two, and merging
gene2a <- dna_seg(list(name="terA", start=1341745, end=ends1[1], strand="+", col="black"))
genes2b <- dna_seg(data.frame(name=c("terB", "terC"), start=c(1684226,
                                                              1609157),
                              end=ends1[2:3], strand=c("-", -1),
                              col=c("black", "black")))
dna_seg2 <- c.dna_seg(gene2a, genes2b)
dna_seg2$fill <- fill
is.dna_seg(dna_seg2)

# for dna_seg3
gene3a <- dna_seg(list(name="terA", start=1341745, end=ends1[1], strand="+", col="black"))
genes3b <- dna_seg(data.frame(name=c("terB", "terC"), start=c(1684226,
                                                              1609157),
                              end=ends1[2:3], strand=c("-", -1),
                              col=c("black", "black")))
dna_seg3 <- c.dna_seg(gene3a, genes3b)
dna_seg3$fill <- fill
is.dna_seg(dna_seg3)

## reading from file
#dna_seg3_file <- system.file('extdata/dna_seg3.tab', package = 'genoPlotR')
#dna_seg3 <- read_dna_seg_from_tab(dna_seg3_file)
#dna_seg3
#is.dna_seg(dna_seg3)


## comparison
## from a data.frame
comparison1 <- as.comparison(data.frame(start1=starts1, end1=ends1,
                                        start2=dna_seg2$start,
                                        end2=dna_seg2$end))
comparison2 <- as.comparison(data.frame(start1=starts1, end1=ends1,
                                        start2=dna_seg2$start,
                                        end2=dna_seg3$end))


## from a file
#comparison2_file <- system.file('extdata/comparison2.tab',
                                #package = 'genoPlotR')
#comparison2 <- read_comparison_from_tab(comparison2_file,
                                        #color_scheme="red_blue")
#is.comparison(comparison1)



## plot
plot_gene_map(dna_segs=list(dna_seg1, dna_seg2, dna_seg3),
              comparisons=list(comparison1, comparison2),
              scale = TRUE,
              dna_seg_scale = TRUE, 
              gene_type = 'arrows', 
              arrow_head_len = 10000)






