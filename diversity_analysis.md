```
library(ggplot2)
library(vegan)

#read in OTU table
otu_table<-read.delim("~/Desktop/mouse_caecal_16S/otu_table_tax_R.txt", header=T, row.names = 1)
dim(otu_table)
#3294 by 116

#read in meta data
meta<-read.delim("~/Desktop/mouse_caecal_16S/mouse_metadata.txt", header=T)
dim(meta)
# 117 by 10

#remove taxonomy
otu_table<-otu_table[,-116]

#transpose table
otu_table_t<-t(otu_table)

#get sample counts
counts<-rowSums(otu_table_t)
which(rowSums(otu_table_t)<3000)

#remove sample M30 and M21
otu_table_t<-otu_table_t[c(-39, -29),]

#plot rarefaction curve
#rarecurve(otu_table_t)

#rarefy data
min(rowSums(otu_table_t))
#11,510
otu_table_t<-rrarefy(otu_table_t, sample=11510)

#calculate PCoA
mouse_pcoa<-capscale(otu_table_t  ~ 1, distance='bray')

#pull out x/y coordinates
mouse.scores<-scores(mouse_pcoa)

#grab only sample coordinates, write to data frame
mouse.coords<-as.data.frame(mouse.scores$sites)

#create sample names as a column
mouse.coords$SampleID<-row.names(mouse.coords)

#map back meta data
mouse.coords<-merge(mouse.coords, meta, by.x='SampleID', by.y='Sample')

#calculate percent variation explained
100*round(mouse_pcoa$CA$eig[1]/sum(mouse_pcoa$CA$eig), 3)
#10.4
100*round(mouse_pcoa$CA$eig[2]/sum(mouse_pcoa$CA$eig), 3)
#7.5

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


female.pcoa<-capscale(mouse.female ~ 1, distance = 'bray')
male.pcoa<-capscale(mouse.male ~ 1, distance = 'bray')


female.scores<-scores(female.pcoa)
female.coords<-as.data.frame(female.scores$sites)
female.coords$SampleID<-row.names(female.coords)
female.coords<-merge(female.coords, meta, by.x='SampleID', by.y='Sample')

ggplot(female.coords, aes(MDS1, MDS2, colour=Sex))+
  geom_point(aes(size=2))+
  theme_bw()

ggplot(female.coords, aes(MDS1, MDS2, colour=Diet_group, label=SampleID))+
  geom_point(aes(size=2))+
  theme_bw()+
  geom_text(hjust=2)


male.scores<-scores(male.pcoa)
male.coords<-as.data.frame(male.scores$sites)
male.coords$SampleID<-row.names(male.coords)
male.coords<-merge(male.coords, meta, by.x='SampleID', by.y='Sample')

ggplot(male.coords, aes(MDS1, MDS2, colour=Sex))+
  geom_point(aes(size=2))+
  theme_bw()

ggplot(male.coords, aes(MDS1, MDS2, colour=Diet_group))+
  geom_point(aes(size=2))+
  theme_bw()

###########Calculate Alpha diversity
mouse.shan<-diversity(otu_table_t, index='shannon')
mouse.obs<-rowSums(otu_table_t>0)
mouse.even<-(diversity(otu_table_t))/mouse.obs

mouse.div<-as.data.frame(cbind(mouse.shan, mouse.obs, mouse.even))
mouse.div$SampleID<-row.names(mouse.div)

#add metadata
mouse.div<-merge(mouse.div, meta, by.x='SampleID', by.y='Sample')

ggplot(mouse.div, aes(Sex, mouse.shan, fill=Diet_group))+
  geom_boxplot()+
  theme_bw()

ggplot(mouse.div, aes(Sex, mouse.obs, fill=Diet_group))+
  geom_boxplot()+
  theme_bw()

ggplot(mouse.div, aes(Sex, mouse.even, fill=Diet_group))+
  geom_boxplot()+
  theme_bw()

#split by sex
mouse.div.split<-split(mouse.div, mouse.div$Sex)

bartlett.test(mouse.div$mouse.obs ~ mouse.div$Sex)
t.test(mouse.div$mouse.obs ~ mouse.div$Sex)
kruskal.test(mouse.div$mouse.obs ~ mouse.div$Sex)


pairwise.t.test(mouse.div.split$F$mouse.even, mouse.div.split$F$Diet_group, p.adjust.method = 'BH')
pairwise.t.test(mouse.div.split$M$mouse.even, mouse.div.split$M$Diet_group, p.adjust.method = 'BH')



#read in PICRUSt KO predictions
pc_ko<-read.delim("~/Desktop/PICRUSt_files/metaG_predict.txt", header=T, row.names=1)


ko_vitK<-c("K01851","K02552","K02361","K02551","K08680","K02549","K14759","K01911","K14760","K01661","K12073","K19222","K02548","K03183")

#get KO of interest
pc_ko_vitk<-pc_ko[row.names(pc_ko) %in% ko_vitK,]
pc_ko_vitk$KO<-row.names(pc_ko_vitk)


#reshape data
library(reshape2)
library(plyr)

pc_ko_vitk_m<-melt(pc_ko_vitk)
pc_ko_vitk_m<-merge(pc_ko_vitk_m, meta, by.x='variable', by.y='Sample')

pc_ko_vitk_m$log_value<-log(pc_ko_vitk_m$value)

ggplot(pc_ko_vitk_m, aes(Diet_group, KO, fill=log_value))+
  geom_tile()+
  theme_bw()+
  facet_wrap(~Sex, scales='free')+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") 

ko_sum<-ddply(pc_ko_vitk_m, c("Diet_group", 'Sex', 'KO'), summarize, mean=mean(value), sd=sd(value), n=length(value), se=sd/n)


ggplot(ko_sum, aes(Diet_group, mean))+
  geom_point()+
  theme_bw()+
  facet_grid(KO~Sex, scales='free')


pc_ko2<-read.delim("~/Desktop/PICRUSt_files/metaG_predict_L2.txt", header=T)
pc_ko2<-read.delim("~/Desktop/PICRUSt_files/metaG_predict_L1.txt", header=T)
pc_ko2_m<-melt(pc_ko2)
pc_ko2_m<-merge(pc_ko2_m, meta, by.x='variable', by.y='Sample')
pc_ko2_m$log_value<-log(pc_ko2_m$value)

ggplot(pc_ko2_m, aes(OTUID, Diet_group, fill=log_value))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  theme_bw()+
  coord_flip()





```
