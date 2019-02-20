```
library(ggplot2)
library(vegan)

otu_table<-read.delim("~/Desktop/otu_table_tax.txt", header=T, row.names = 1)
dim(otu_table)
#3294 by 116

meta<-read.delim("~/Desktop/mouse_metadata.txt", header=T)
dim(meta)
# 117 by 10

#remove taxonomy
otu_table<-otu_table[,-116]

#transpose table
otu_table_t<-t(otu_table)

#get sample counts
counts<-rowSums(otu_table_t)
which(rowSums(otu_table_t)<1000)

#remove sample M30
otu_table_t<-otu_table_t[-39,]

#rarefy data
min(rowSums(otu_table_t))
#2893

otu_table_t<-rrarefy(otu_table_t, sample=2893)

mouse_pcoa<-capscale(otu_table_t  ~ 1, distance='bray')

mouse.scores<-scores(mouse_pcoa)
mouse.coords<-as.data.frame(mouse.scores$sites)
mouse.coords$SampleID<-row.names(mouse.coords)
mouse.coords<-merge(mouse.coords, meta, by.x='SampleID', by.y='Sample')

#plot by group
ggplot(mouse.coords, aes(MDS1, MDS2, colour=Diet_group))+
  geom_point(aes(size=1.5))+
  theme_bw()

#plot by sex
ggplot(mouse.coords, aes(MDS1, MDS2, colour=Sex))+
  geom_point(aes(size=1.5))+
  theme_bw()

#plot by sex and group
ggplot(mouse.coords, aes(MDS1, MDS2, colour=Diet_group, shape=Sex))+
  geom_point(aes(size=1.5))+
  theme_bw()


##Subsample OTU table by sex
#split map
meta_split<-split(meta, meta$Sex)

#make sex-specifc tables
mouse.male<-otu_table_t[row.names(otu_table_t) %in% meta_split$M$Sample,]
mouse.female<-otu_table_t[row.names(otu_table_t) %in% meta_split$F$Sample,]



#Make new map
tiss_map<- subset(map, Sample_Type=="Gills" | Sample_Type=="Stomach")
#Get rows that match IDs in map
tiss_df<-weighted_u[rownames(weighted_u) %in% tiss_map$SampleID, ]
#gets columns that match IDs in map
tiss_df1<-tiss_df[ ,colnames(tiss_df) %in% tiss_map$SampleID]
```
