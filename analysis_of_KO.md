## Visualize KO
```
#read in PICRUSt KO predictions
pc_ko<-read.delim("mouse_caecal_16S/PICRUSt_files/metaG_predict.txt", header=T, row.names=1)

ko_vitK<-c("K01851","K02552","K02361","K02551","K08680","K02549","K14759","K01911","K14760","K01661","K12073","K19222","K02548","K03183")

#get KO of interest
pc_ko_vitk<-pc_ko[row.names(pc_ko) %in% ko_vitK,]
pc_ko_vitk$KO<-row.names(pc_ko_vitk)


#reshape data
library(reshape2)
library(plyr)

pc_ko_vitk_m<-melt(pc_ko_vitk)
pc_ko_vitk_m<-merge(pc_ko_vitk_m, meta, by.x='variable', by.y='SampleID')

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


pc_ko2<-read.delim("mouse_caecal_16S/PICRUSt_files/metaG_predict_L2.txt", header=T)
pc_ko2<-read.delim("mouse_caecal_16S/PICRUSt_files/metaG_predict_L1.txt", header=T)
pc_ko2_m<-melt(pc_ko2)
pc_ko2_m<-merge(pc_ko2_m, meta, by.x='variable', by.y='Sample')
pc_ko2_m$log_value<-log(pc_ko2_m$value)

ggplot(pc_ko2_m, aes(OTUID, Diet_group, fill=log_value))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  theme_bw()+
  coord_flip()
  ```


## Global DESeq Analysis

```
library(DESeq2)

#read in meta data
meta<-read.delim("mouse_caecal_16S/mouse_metadata.txt", header=T)

#subset to only include diet as a treatment
meta<-meta[,c(1,8,9)]

#read in PICRUSt KO predictions
pc_ko<-read.delim("mouse_caecal_16S/PICRUSt_files/metaG_predict.txt", header=T, row.names=1)

#change numbers to integers by rounding
pc_ko<-as.integer(pc_ko)

#remove extra samples
row.names(meta)<-meta$SampleID
meta<-meta[row.names(meta) %in% names(pc_ko),]

#reorder sample names to match column names in KO data
meta<-meta[names(pc_ko),,drop=F]

#convert to DESEQ2 format
meta_deseq<-DESeqDataSetFromMatrix(countData = pc_ko, colData = meta, design=~Diet_group)

#run deseq2 analysis
dds<-DESeq(meta_deseq)

#write results to a table
deseq_results<-results(dds)
deseq_results<-as.data.frame(deseq_results)
deseq_results$KO<-row.names(deseq_results)
write.table(deseq_results, 'deseq_results_full.txt',  row.names=F, quote=F,, sep='\t')

#create list of desired KO's assocaited with vitamin K metabolism
ko_vitK<-c("K01851","K02552","K02361","K02551","K08680","K02549","K14759","K01911","K14760","K01661","K12073","K19222","K02548","K03183")

#get KO of interest
pc_ko_vitk<-deseq_results[deseq_results$KO %in% ko_vitK,]
pc_ko_vitk$KO<-row.names(pc_ko_vitk)
write.table(pc_ko_vitk, 'deseq_ko_vitK.txt', quote=F, row.names = F, sep='\t')
````

## Sex-specific deseq analysis

```
#split map by sex
meta_split<-split(meta, meta$Sex)

#create meta data and KO tables for each sex
meta_male<-meta_split$M
meta_female<-meta_split$F
row.names(meta_female)<-meta_female$SampleID
row.names(meta_male)<-meta_male$SampleID

#create a KO table for each sex
ko_male<-pc_ko[,names(pc_ko) %in% meta_male$SampleID]
ko_female<-pc_ko[,names(pc_ko) %in% meta_female$SampleID]

#remove missing female samples from meta data
meta_female<-meta_female[names(ko_female),,drop=T]

#reorder sample names
meta_female<-meta_female[names(ko_female),,drop=F]
meta_male<-meta_male[names(ko_male),,drop=F]

#create DESeq2 objects
deseq_male<-DESeqDataSetFromMatrix(countData =ko_male, colData = meta_male, design=~Diet_group)
deseq_female<-DESeqDataSetFromMatrix(countData =ko_female, colData = meta_female, design=~Diet_group)

#run DESEQ2
dds_male<-DESeq(deseq_male)
dds_female<-DESeq(deseq_female)

#write results to a table
deseq_results_male<-results(dds_male)
deseq_results_male<-as.data.frame(deseq_results_male)
deseq_results_male$KO<-row.names(deseq_results_male)
write.table(deseq_results_male, 'deseq_results_full_male.txt',  row.names=F, quote=F,, sep='\t')

#create list of desired KO's assocaited with vitamin K metabolism
ko_vitK<-c("K01851","K02552","K02361","K02551","K08680","K02549","K14759","K01911","K14760","K01661","K12073","K19222","K02548","K03183")

#get KO of interest
pc_ko_vitk_male<-deseq_results_male[deseq_results_male$KO %in% ko_vitK,]
pc_ko_vitk_male$KO<-row.names(pc_ko_vitk_male)
write.table(pc_ko_vitk_male, 'deseq_ko_vitK_male.txt', quote=F, row.names = F, sep='\t')


deseq_results_female<-results(dds_female)
deseq_results_female<-as.data.frame(deseq_results_female)
deseq_results_female$KO<-row.names(deseq_results_female)
write.table(deseq_results_female, 'deseq_results_full_female.txt',  row.names=F, quote=F,, sep='\t')

#create list of desired KO's assocaited with vitamin K metabolism
ko_vitK<-c("K01851","K02552","K02361","K02551","K08680","K02549","K14759","K01911","K14760","K01661","K12073","K19222","K02548","K03183")

#get KO of interest
pc_ko_vitk_female<-deseq_results_female[deseq_results_female$KO %in% ko_vitK,]
pc_ko_vitk_female$KO<-row.names(pc_ko_vitk_female)
write.table(pc_ko_vitk_female, 'deseq_ko_vitK_female.txt', quote=F, row.names = F, sep='\t')
