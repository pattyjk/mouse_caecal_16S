## Beta diversity

```
library(ggplot2)
library(vegan)
library(gridExtra)
##will use zOTUS for the analysis

#read in zOTU table
otu_table<-read.delim("mouse_caecal_16S/UNOISE_files/zotus_table.txt", header=T, row.names = 1)
dim(otu_table)
#1214 by 115

#read in meta data
meta<-read.delim("mouse_caecal_16S/mouse_metadata.txt", header=T)
dim(meta)
# 117 by 35

#transpose table
otu_table_t<-t(otu_table)

#get sample counts
which(rowSums(otu_table_t)<3000)

#remove sample M30 and M21
otu_table_t<-otu_table_t[c(-39, -29),]

#plot rarefaction curve
#rarecurve(otu_table_t)

#rarefy data
min(rowSums(otu_table_t))
#11,558
otu_table_t<-rrarefy(otu_table_t, sample=11558)

#calculate PCoA
mouse_pcoa<-capscale(otu_table_t  ~ 1, distance='bray')

#pull out x/y coordinates
mouse.scores<-scores(mouse_pcoa)

#grab only sample coordinates, write to data frame
mouse.coords<-as.data.frame(mouse.scores$sites)

#create sample names as a column
mouse.coords$SampleID<-row.names(mouse.coords)

#map back meta data
mouse.coords<-merge(mouse.coords, meta, by.x='SampleID', by.y='SampleID')

#calculate percent variation explained
100*round(mouse_pcoa$CA$eig[1]/sum(mouse_pcoa$CA$eig), 3)
#10.2
100*round(mouse_pcoa$CA$eig[2]/sum(mouse_pcoa$CA$eig), 3)
#7.5


#plot by sex
total<-ggplot(mouse.coords, aes(MDS1, MDS2, colour=Sex))+
  geom_point(aes(size=2))+
  theme_bw()+
  xlab("PC1- 10.2%")+
  ylab("PC2- 7.5%")+
  ggtitle("A) All mice")


## Subsample OTU table by sex
#split map
meta_split<-split(meta, meta$Sex)

#make sex-specifc tables
mouse.male<-otu_table_t[row.names(otu_table_t) %in% meta_split$M$Sample,]
mouse.female<-otu_table_t[row.names(otu_table_t) %in% meta_split$F$Sample,]

#calculate PCoA
female.pcoa<-capscale(mouse.female ~ 1, distance = 'bray')
male.pcoa<-capscale(mouse.male ~ 1, distance = 'bray')

#extract coordinates
female.scores<-scores(female.pcoa)
female.coords<-as.data.frame(female.scores$sites)
female.coords$SampleID<-row.names(female.coords)
female.coords<-merge(female.coords, meta, by.x='SampleID', by.y='SampleID')

male.scores<-scores(male.pcoa)
male.coords<-as.data.frame(male.scores$sites)
male.coords$SampleID<-row.names(male.coords)
male.coords<-merge(male.coords, meta, by.x='SampleID', by.y='SampleID')

#calculate percent variation explained, female mice
100*round(female.pcoa$CA$eig[1]/sum(female.pcoa$CA$eig), 3)
#10.4
100*round(female.pcoa$CA$eig[2]/sum(female.pcoa$CA$eig), 3)
#8.1

#calculate percent variation explained, male mice
100*round(male.pcoa$CA$eig[1]/sum(male.pcoa$CA$eig), 3)
#9.9
100*round(male.pcoa$CA$eig[2]/sum(male.pcoa$CA$eig), 3)
#8.5

#plot female PCoA
fem<-ggplot(female.coords, aes(MDS1, MDS2, colour=Diet_group))+
  geom_point(aes(size=2))+
  theme_bw()+
  #geom_text(hjust=2)+
  ggtitle("B) Female mice")+
  xlab("PC1- 10.4%")+
  ylab("PC2- 8.1%")


#plot male PCoA
men<-ggplot(male.coords, aes(MDS1, MDS2, colour=Diet_group))+
  geom_point(aes(size=2))+
  #geom_text(hjust=2)+
  theme_bw()+
  ggtitle("C) Male mice")+
  xlab("PC1- 9.9%")+
  ylab("PC2- 8.5%")       

#plot on same plot
grid.arrange(total, fem, men, ncol=3)
````

## Test for significance of beta diversity between sexes/diet groups
```
library(vegan)

#read in zOTU table
otu_table<-read.delim("mouse_caecal_16S/UNOISE_files/zotus_table.txt", header=T, row.names = 1)
dim(otu_table)
#1214 by 115

#read in meta data
meta<-read.delim("mouse_caecal_16S/mouse_metadata.txt", header=T)
dim(meta)
# 117 by 35

#transpose table
otu_table_t<-t(otu_table)

#get sample counts
which(rowSums(otu_table_t)<3000)

#remove sample M30 and M21
otu_table_t<-otu_table_t[c(-39, -29),]

#order metadata by OTU table order
row.names(meta)<-meta$SampleID
meta<-meta[row.names(otu_table_t),]

#rarefy data
otu_rare<-rrarefy(otu_table_t, sample=11558)

#calculate adons
adonis(otu_rare ~ Diet_group + Sex, data=meta, permutations =10000, method = 'bray')

#            Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
#Diet_group   5    1.1224 0.22449  1.7820 0.07409 9.999e-05 ***
#Sex          1    0.6742 0.67416  5.3516 0.04450 9.999e-05 ***
#Residuals  106   13.3533 0.12597         0.88141              


#split meta data by sex and create sex-specific OTU tables
meta_split<-split(meta, meta$Sex)

male_otu<-otu_rare[row.names(meta_split$M),]
female_otu<-otu_rare[row.names(meta_split$F),]

#adonis of diet group by for male/female

adonis(male_otu ~ Diet_group, data=meta_split$M, permutations =10000, method = 'bray')
#           Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
#Diet_group  5    1.0651 0.21302  1.7914 0.14939 9.999e-05 ***
#Residuals  51    6.0648 0.11892         0.85061              
#Total      56    7.1299                 1.00000 

adonis(female_otu ~ Diet_group, data=meta_split$F, permutations =10000, method = 'bray')
#           Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
#Diet_group  5    1.2195 0.24389  1.9884 0.16586 9.999e-05 ***
#Residuals  50    6.1330 0.12266         0.83414              
#Total      55    7.3524                 1.00000             

```

## Alpha diversity
````
#Calculate Alpha diversity
mouse.shan<-diversity(otu_table_t, index='shannon')
mouse.obs<-rowSums(otu_table_t>0)
mouse.even<-(diversity(otu_table_t))/mouse.obs

mouse.div<-as.data.frame(cbind(mouse.shan, mouse.obs, mouse.even))
mouse.div$SampleID<-row.names(mouse.div)

#add metadata
mouse.div<-merge(mouse.div, meta, by.x='SampleID', by.y='SampleID')

shannon<-ggplot(mouse.div, aes(Sex, mouse.shan, fill=Diet_group))+
  geom_boxplot()+
  theme_bw()+
  ylab("Shannon Diversity")+
  xlab("")

richness<-ggplot(mouse.div, aes(Sex, mouse.obs, fill=Diet_group))+
  geom_boxplot()+
  theme_bw()+
  ylab("Bacterial Richness")+
  xlab("")

evenness<-ggplot(mouse.div, aes(Sex, mouse.even, fill=Diet_group))+
  geom_boxplot()+
  theme_bw()+
  ylab("Pielou's Evenness")+
  xlab("")

grid.arrange(shannon, richness, evenness, ncol=3)
```

## Test for significance of alpha diversity between sexes
```
#test for differences in variances between sexes

bartlett.test(mouse.div$mouse.obs ~ mouse.div$Sex)
#Bartlett's K-squared = 6.3181, df = 1, p-value = 0.01195

bartlett.test(mouse.div$mouse.shan ~ mouse.div$Sex)
#Bartlett's K-squared = 1.8497, df = 1, p-value = 0.1738

bartlett.test(mouse.div$mouse.even ~ mouse.div$Sex)
#Bartlett's K-squared = 91.56, df = 1, p-value < 2.2e-16

#t-test to test for significance between male and female mice
t.test(mouse.div$mouse.obs ~ mouse.div$Sex)
#t = 3.3878, df = 101.11, p-value = 0.001005

t.test(mouse.div$mouse.even ~ mouse.div$Sex)
#t = -1.7215, df = 61.886, p-value = 0.09016

t.test(mouse.div$mouse.shan ~ mouse.div$Sex)
#t = 3.0707, df = 108.08, p-value = 0.002702
```

## Test for significance of alpha diversity between diet groups
```
#split data by sex
mouse.div.split<-split(mouse.div, mouse.div$Sex)

#calculate pairwise t-test for female mice
pairwise.t.test(mouse.div.split$F$mouse.even, mouse.div.split$F$Diet_group, p.adjust.method = 'BH')
#no significant differences

pairwise.t.test(mouse.div.split$F$mouse.shan, mouse.div.split$F$Diet_group, p.adjust.method = 'BH')
#no significant differences

pairwise.t.test(mouse.div.split$F$mouse.obs, mouse.div.split$F$Diet_group, p.adjust.method = 'BH')
#no significant differences


#calculate pairwise t-test for male mice
pairwise.t.test(mouse.div.split$M$mouse.even, mouse.div.split$M$Diet_group, p.adjust.method = 'BH')
#no significant differences

pairwise.t.test(mouse.div.split$M$mouse.shan, mouse.div.split$M$Diet_group, p.adjust.method = 'BH')
#no significant differences

pairwise.t.test(mouse.div.split$M$mouse.obs, mouse.div.split$M$Diet_group, p.adjust.method = 'BH')
#no significant differences
```
