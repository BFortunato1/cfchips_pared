#compile summary metrics useful for qc / technical optimization
library(stringr)

args = commandArgs( trailingOnly = TRUE )

peaks_file = args[1]
mapping_file = args[2]
enrich_file = args[3] 
complexity_file = args[4]
fragments_file = args[5]
out_file = args[6]

peaks = read.table(peaks_file, header=T, sep=",")
mapping = read.table(mapping_file, header=T, sep=",")
enrich = read.table(enrich_file, header=T, sep="\t")
complexity = read.table(complexity_file, header=T, sep="\t")
fragments = read.table(fragments_file, header=T, sep="\t")

peaks = peaks[,1:5] #first 5 columns have metrics of interest
colnames(peaks)[1] = "name"
peaks$name = str_replace(peaks$name, ".rep1","")
colnames(mapping)[1] = "name"

out = merge(mapping, peaks, by = "name")
out = merge(out, enrich, by = "name")
out = merge(out, complexity, by = "name")
out = merge(out, fragments, by = "name")

write.table(out, out_file, row.names=F, col.names=T, quote=F, sep="\t")

