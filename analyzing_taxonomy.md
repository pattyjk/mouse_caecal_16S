## Analyze taxonomic data from SINTAX
```
library(ggplot2)
library(plyr)
library(tidyr)
library(reshape2)

#read in OTU table
otu_table<-read.delim("mouse_caecal_16S/UNOISE_files/zotus_table.txt", header=T, row.names=1)

#change OTU table to relative abundanc, can ignor this step if absolute abundance is needed
otu_table<-sweep(otu_table, 2, colSums(otu_table), '/')

#add OTUs as a column
otu_table$X.OTU.ID<-row.names(otu_table)

#read in taxonomic annotations from SINTAX
tax_anno<-read.delim("mouse_caecal_16S/UNOISE_files/zotus_sintax_annotations.txt", header=F)
tax_anno<-tax_anno[,c(-3,-2)]

#add taxonomy to OTU table
otu_table_tax<-merge(otu_table, tax_anno, by.x='X.OTU.ID', by.y='V1')

#split tax annotations by type, some will be missing taxonomic annotations
otu_table_tax<-separate(data = otu_table_tax, col = V4, sep = ',', into = c("Kingdom", 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'))

#for each sample, get total abundance of each taxonimic group
otu_table_m<-melt(otu_table_tax)
otu_table_sum<-ddply(otu_table_m, c("Kingdom", 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'variable'), summarize, total=sum(value))

#data is read to be plotted
```
