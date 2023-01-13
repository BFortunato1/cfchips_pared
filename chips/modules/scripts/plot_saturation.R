# plot library saturation data from PICARD MarkDuplicates
# generates a plot of predicted number of fragments sampled at a given sequencing depth
# note: only considers paired reads

# the ESTIMATED_LIBRARY_SIZE comes from the Lander Waterman equation:

#     *     C/X = 1 - exp( -N/X )
#     * where
#     *     X = number of distinct molecules in library
#     *     N = number of read pairs
#     *     C = number of distinct fragments observed in read pairs

# see https://github.com/broadinstitute/picard/blob/master/src/main/java/picard/sam/DuplicationMetrics.java#L115

# Sylvan Baca 3/29/21

library(ggplot2)
library(stringr)

args = commandArgs(trailingOnly=T)

dup.file = args[1]
plot.file = args[2]
complexity.file = args[3]

#---tmp
#dup.file="analysis/align/KamP1u1R/KamP1u1R_dup.metrics.txt"
#plot.file='tmp.pdf'
#complexity.file='tmp.tsv'
#---end tmp
sample.name = dup.file %>% str_replace(".*align/","") %>% str_replace("/.*","")

stats = read.table(dup.file, nrows=1, header=T)
sat = read.table(dup.file, skip=10, header=T)

#simulated sequencing depths
sat$depth = stats$READ_PAIRS_EXAMINED * sat$BIN

#how many fragments are predicted to be sampled at these depths?
sat$sampled_fragments = (stats$READ_PAIRS_EXAMINED-stats$READ_PAIR_DUPLICATES) * sat$VALUE

#how deeply must we sequence to sample 90% of fragments?
s90 = round(min(sat$depth[sat$sampled_fragments >= stats$ESTIMATED_LIBRARY_SIZE * 0.9]), digits=0)

##subset sat to go no further than twice the depth needed to sequence 90% of fragments (to improve display of plot)
##sat = subset(sat, depth < s90)

ggplot(sat, aes(y=sampled_fragments, x=depth)) + 
  geom_line() +
  ylim(0,max(sat$sampled_fragments)) +
  xlim(0,5e7) +
  ggtitle(paste0(sample.name, "\nSequenced depth: ", stats$READ_PAIRS_EXAMINED * 2, 
    "\nEst complexity (Lander-Waterman): ", stats$ESTIMATED_LIBRARY_SIZE,
    "\nDepth to sequence 90% of fragments: ", s90)) +
  theme_classic() +
  theme(plot.title = element_text(size = 8)) +
  geom_vline(xintercept=s90, color="blue", linetype="dashed")

ggsave(plot.file, height=3, width=4)

#write est library complexity to file
out=data.frame(name=sample.name, complexity=stats$ESTIMATED_LIBRARY_SIZE)

write.table(out, file=complexity.file, quote=F, row.names=F, sep="\t")

