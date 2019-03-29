## Analyze taxonomic data from SINTAX
```
library(ggplot2)
library(plyr)
library(tidyr)
library(reshape2)

#read in OTU table
otu_table<-read.delim("mouse_caecal_16S/UNOISE_files/zotus_table.txt", header=T, row.names=1)

#change OTU table to relative abundanc, can ignore this step if absolute abundance is needed
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

#Reshape table
otu_table_m<-melt(otu_table_tax)

#read in metadata
meta<-read.delim("mouse_caecal_16S/mouse_metadata.txt", header=T)

#add metadata to taxonomic file
otu_table_m<-merge(otu_table_m, meta, by.x='variable', by.y='SampleID')


#for each sample, get total abundance of each taxonimic group and take the mean
otu_table_sum<-ddply(otu_table_m, c("Kingdom", 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'variable', 'Sex', 'Diet_group'), summarize, total=sum(value), mean=mean(value))

```

## Plot data
```
```

## Identify differentiall abundant OTUs bewteen male/female mice
```
library(dplyr)
#read in OTU table
otu_table<-read.delim("mouse_caecal_16S/UNOISE_files/zotus_table.txt", header=T, row.names=1)

#read in metadata
meta<-read.delim("mouse_caecal_16S/mouse_metadata.txt", header=T)
row.names(meta)<-meta$SampleID

#remove extra samples (those without enough depth)
otu_table<-otu_table[,names(otu_table) %in% meta$SampleID]


#match the sample IDs in map to columns
otu_table<-otu_table[,row.names(meta),drop=F]

#transpose OTU table
otu_table_t<-as.data.frame(t(otu_table))

#split table into a list of columns
testy2<-split(c(t(otu_table_t)), rep(seq(1, ncol(otu_table_t)), each = 1))

#calculate kruskal-wallis test 
outy<-lapply(testy2, function(x) kruskal.test(x, meta$Sex))

#convert lists to tables
outy2<-lapply(outy, function(x) as.data.frame(t(unlist(x))))

#catenate all tables
outy3<-ldply(outy2, data.frame)

#add OTU information back in
outy3$OTU<-names(otu_table_t)

#read in taxonomic annotations from SINTAX
tax_anno<-read.delim("mouse_caecal_16S/UNOISE_files/zotus_sintax_annotations.txt", header=F)
tax_anno<-tax_anno[,c(-3,-2)]

#add taxonomic annotations to Kruskal-wallis results
outy4<-merge(outy3, tax_anno, by.x='OTU', by.y='V1')

#fix header names
outy4<-outy4[,c(-7,-2)]
names(outy4)<-c("OTU", "KW_Statistic", "DF", "P_value", "Test", "Taxonomy")

#write table to file and read back in
write.table(outy4, 'mouse_caecal_16S/sex_KW_output.txt', row.names = F, quote=F, sep='\t')
sex_kw_output<-read.delim('mouse_caecal_16S/sex_KW_output.txt', header=T)

#correct p-values
sex_kw_output$p.value.fdr<-p.adjust(sex_kw_output$P_value, method='fdr')
sex_kw_output$p.value.hochberg<-p.adjust(sex_kw_output$P_value, method='BH')
sex_kw_output$p.value.bonferonni<-p.adjust(sex_kw_output$P_value, method='bonferroni')
```
