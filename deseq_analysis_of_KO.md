```
library(DESeq2)

#read in meta data
meta<-read.delim("mouse_caecal_16S/mouse_metadata.txt", header=T)

#subset to only include diet as a treatment
meta<-meta[,c(1,9)]

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
```
