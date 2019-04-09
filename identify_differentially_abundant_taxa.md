## Identify differentiall abundant OTUs bewteen male/female mice
```
library(dplyr)
library(plyr)
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

## Identify taxa in male and female mice that respond to diets
```
#read in OTU table
otu_table<-read.delim("mouse_caecal_16S/UNOISE_files/zotus_table.txt", header=T, row.names=1)

#read in metadata
meta<-read.delim("mouse_caecal_16S/mouse_metadata.txt", header=T)
row.names(meta)<-meta$SampleID

#remove extra samples (those without enough sequencing depth)
otu_table<-otu_table[,names(otu_table) %in% meta$SampleID]

#match the sample IDs in map to columns
otu_table<-otu_table[,row.names(meta),drop=F]

#transpose OTU table
otu_table_t<-as.data.frame(t(otu_table))

#split map, make male/female OTU tables
meta_split<-split(meta, meta$Sex)
otu_male<-otu_table_t[row.names(otu_table_t) %in% meta_split$M$SampleID,]
otu_female<-otu_table_t[row.names(otu_table_t) %in% meta_split$F$SampleID,]

#split tables into a list of columns
otu_male2<-split(c(t(otu_male)), rep(seq(1, ncol(otu_male)), each = 1))
otu_female2<-split(c(t(otu_female)), rep(seq(1, ncol(otu_female)), each = 1))

#calculate kruskal-wallis test 
otu_male3<-lapply(otu_male2, function(x) kruskal.test(x, meta_split$M$Diet_group))
otu_female3<-lapply(otu_female2, function(x) kruskal.test(x, meta_split$F$Diet_group))

#convert lists to tables
otu_male3<-lapply(otu_male3, function(x) as.data.frame(t(unlist(x))))
otu_female3<-lapply(otu_female3, function(x) as.data.frame(t(unlist(x))))

#catenate all tables
otu_male3<-ldply(otu_male3, data.frame)
otu_female3<-ldply(otu_female3, data.frame)

#add OTU information back in
otu_male3$OTU<-names(otu_male)
otu_female3$OTU<-names(otu_female)

#read in taxonomic annotations from SINTAX
tax_anno<-read.delim("mouse_caecal_16S/UNOISE_files/zotus_sintax_annotations.txt", header=F)
tax_anno<-tax_anno[,c(-3,-2)]

#add taxonomic annotations to Kruskal-wallis results
otu_male3<-merge(otu_male3, tax_anno, by.x='OTU', by.y='V1')
otu_female3<-merge(otu_female3, tax_anno, by.x='OTU', by.y='V1')

#fix header names
otu_male3<-otu_male3[,c(-7,-2)]
otu_female3<-otu_female3[,c(-7,-2)]

names(otu_male3)<-c("OTU", "KW_Statistic", "DF", "P_value", "Test", "Taxonomy")
names(otu_female3)<-c("OTU", "KW_Statistic", "DF", "P_value", "Test", "Taxonomy")

#write table to file and read back in
write.table(otu_male3, 'mouse_caecal_16S/male_diet_KW_output.txt', row.names = F, quote=F, sep='\t')
male_diet_kw_output<-read.delim('mouse_caecal_16S/male_diet_KW_output.txt', header=T)

write.table(otu_female3, 'mouse_caecal_16S/female_diet_KW_output.txt', row.names = F, quote=F, sep='\t')
female_diet_kw_output<-read.delim('mouse_caecal_16S/female_diet_KW_output.txt', header=T)

#correct p-values
female_diet_kw_output$p.value.fdr<-p.adjust(female_diet_kw_output$P_value, method='fdr')
female_diet_kw_output$p.value.hochberg<-p.adjust(female_diet_kw_output$P_value, method='BH')
female_diet_kw_output$p.value.bonferonni<-p.adjust(female_diet_kw_output$P_value, method='bonferroni')

male_diet_kw_output$p.value.fdr<-p.adjust(male_diet_kw_output$P_value, method='fdr')
male_diet_kw_output$p.value.hochberg<-p.adjust(male_diet_kw_output$P_value, method='BH')
male_diet_kw_output$p.value.bonferonni<-p.adjust(male_diet_kw_output$P_value, method='bonferroni')

#write output to file
write.table(male_diet_kw_output, "male_diet_kw_output.txt", row.names=F, quote=F, sep='\t')
write.table(female_diet_kw_output, "female_diet_kw_output.txt", row.names=F, quote=F, sep='\t')
```


## Plot data for male/female
```
library(ggplot2)
library(reshape2)
library(vegan)
library(tidyr)
library(plyr)

#read in metadata
meta<-read.delim("mouse_caecal_16S/mouse_metadata.txt", header=T)

#get taxa that significantly vary between sexes
taxa_sex<-which(sex_kw_output$p.value.fdr <0.05)
sig_taxa<-sex_kw_output[row.names(sex_kw_output) %in% taxa_sex,]

#read in OTU table
otu_table<-read.delim("mouse_caecal_16S/UNOISE_files/zotus_table.txt", header=T, row.names=1)

#rarefy OTU table
otu_table_t<-t(otu_table)
min(rowSums(otu_table_t))
#11558
otu_table_t<-rrarefy(otu_table_t, sample=11558)
otu_table_t<-t(otu_table_t)


#make significant OTU table
otu_table_sig<-otu_table[row.names(otu_table_t) %in% sig_taxa$OTU,]
otu_table_sig$OTU<-row.names(otu_table_sig)

#reshape OTU table, add meta data
otu_sig_m<-melt(otu_table_sig, id.vars = c('OTU'))
otu_sig_m<-merge(otu_sig_m, meta, by.x='variable', by.y='SampleID')

#get mean abundance
otu_sig_sum<-ddply(otu_sig_m, c("Sex", "OTU"), summarize, mean=mean(value), sd=sd(value), n=length(value), se=sd/n)

#add taxonomy to data
tax_anno<-read.delim("mouse_caecal_16S/UNOISE_files/zotus_sintax_annotations.txt", header=F)
tax_anno<-tax_anno[,c(-3,-2)]
otu_sig_sum<-merge(otu_sig_sum, tax_anno, by.x='OTU', by.y='V1')

#split taxonomy into phylum, class, ect...
otu_sig_sum<-separate(otu_sig_sum, V4, sep=',', into = c("Kingdom", 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'))

#get only top taxa with abundance >100 reads/sample
abundant_taxa<-which(otu_sig_sum$mean>100)
abundant_sig_taxa<-otu_sig_sum[row.names(otu_sig_sum) %in% abundant_taxa,]


#plot data
ggplot(abundant_sig_taxa, aes(OTU, mean, fill=Genus))+
  geom_bar(stat='identity')+
  theme_bw()+
  facet_wrap(~Sex)+
  coord_flip()+
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se, width=0.2), stat="identity")+
  ylab("Mean OTU Abundance")+
  xlab("")
  ```
  
  ## Plot data for diet
```
library(ggplot2)
library(reshape2)
library(vegan)
library(tidyr)
library(plyr)

#read in metadata
meta<-read.delim("mouse_caecal_16S/mouse_metadata.txt", header=T)

#read in OTU table
otu_table<-read.delim("mouse_caecal_16S/UNOISE_files/zotus_table.txt", header=T, row.names=1)

#rarefy OTU table
otu_table_t<-t(otu_table)
min(rowSums(otu_table_t))
#11558
otu_table_t<-rrarefy(otu_table_t, sample=11558)
otu_table_t<-t(otu_table_t)

#make significant OTU table for male/female
sig_male<-which(male_diet_kw_output$P_value<0.05)
sig_female<-which(female_diet_kw_output$P_value<0.05)

sig_male<-male_diet_kw_output[row.names(male_diet_kw_output) %in% sig_male,]
sig_female<-female_diet_kw_output[row.names(female_diet_kw_output) %in% sig_female,]

otu_sig_male<-as.data.frame(otu_table_t[row.names(otu_table_t) %in% sig_male$OTU,])
otu_sig_female<-as.data.frame(otu_table_t[row.names(otu_table_t) %in% sig_female$OTU,])

otu_sig_male$OTU<-row.names(otu_sig_male)
otu_sig_female$OTU<-row.names(otu_sig_female)

#reshape OTU table, add meta data
sig_male_m<-melt(otu_sig_male, id.vars = c('OTU'))
sig_female_m<-melt(otu_sig_female, id.vars = c('OTU'))

sig_male_m<-merge(sig_male_m, meta, by.x='variable', by.y='SampleID')
sig_female_m<-merge(sig_female_m, meta, by.x='variable', by.y='SampleID')

#get mean abundance
otu_sig_male_sum<-ddply(sig_male_m, c("Sex", "OTU", "Diet_group"), summarize, mean=mean(value), sd=sd(value), n=length(value), se=sd/n)
otu_sig_female_sum<-ddply(sig_female_m, c("Sex", "OTU", "Diet_group"), summarize, mean=mean(value), sd=sd(value), n=length(value), se=sd/n)

#add taxonomy to data
tax_anno<-read.delim("mouse_caecal_16S/UNOISE_files/zotus_sintax_annotations.txt", header=F)
tax_anno<-tax_anno[,c(-3,-2)]
otu_sig_male_sum<-merge(otu_sig_male_sum, tax_anno, by.x='OTU', by.y='V1')
otu_sig_female_sum<-merge(otu_sig_female_sum, tax_anno, by.x='OTU', by.y='V1')

#split taxonomy into phylum, class, ect...
otu_sig_male_sum<-separate(otu_sig_male_sum, V4, sep=',', into = c("Kingdom", 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'))
otu_sig_female_sum<-separate(otu_sig_female_sum, V4, sep=',', into = c("Kingdom", 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'))

#get only most abundant taxa
male_abun<-which(otu_sig_male_sum$mean>100)
female_unbun<-which(otu_sig_female_sum$mean>100)

abun_otu_sig_male_sum<-otu_sig_male_sum[row.names(otu_sig_male_sum) %in% male_abun,]
abun_otu_sig_female_sum<-otu_sig_female_sum[row.names(otu_sig_female_sum) %in% male_abun,]
  

#plot the data
ggplot(abun_otu_sig_female_sum, aes(OTU, mean, fill=Genus))+
  theme_bw()+
  geom_bar(stat='identity')+
  coord_flip()+
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se, width=0.2), stat="identity")+
  ylab("Mean OTU Abundance")+
  xlab("")+
  facet_wrap(~Diet_group)+
  ggtitle("Female Mice")
  
ggplot(abun_otu_sig_male_sum, aes(OTU, mean, fill=Genus))+
  theme_bw()+
  geom_bar(stat='identity')+
  coord_flip()+
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se, width=0.2), stat="identity")+
  ylab("Mean OTU Abundance")+
  xlab("")+
  ggtitle("Male Mice")+
  facet_wrap(~Diet_group)
```
