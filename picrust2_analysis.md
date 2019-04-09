## Analysis of zOTUs with PICRUSt2
```
#installation if neccessary 
#conda create -n picrust2 -c bioconda -c conda-forge picrust2
#source activate picrust2

#convert zOTU table to biom format
biom convert -i mouse_caecal_16S/UNOISE_files/zotus_table.txt -o mouse_caecal_16S/UNOISE_files/zotus_table.biom --to-hdf5 --table-type='OTU table'

#run PICRUSt 2 full pipeline
picrust2_pipeline.py -s mouse_caecal_16S/UNOISE_files/zotus.fa -i mouse_caecal_16S/UNOISE_files/zotus_table.biom -o mouse_caecal_16S/picrust2_output --threads 4
```

## Analysis in R
```
library(plyr)
library(ggplot2)
library(reshape2)

#read in data
ko<-read.delim('mouse_caecal_16S/picrust2_output/KO_metagenome_out/pred_metagenome_unstrat.tsv', header=T, row.names=1)

#get KO of interest
ko_vitK<-c("K01851","K02552","K02361","K02551","K08680","K02549","K14759","K01911","K14760","K01661","K12073","K19222","K02548","K03183")

ko_vk<-ko[row.names(ko) %in% ko_vitK,]
ko_vk$KO<-row.names(ko_vk)

#reshape data
ko_vk_m<-melt(ko_vk)

#read in meta data
meta<-read.delim("mouse_caecal_16S/mouse_metadata.txt", header=T)

#add metadata to KO data
ko_vk_m<-merge(ko_vk_m, meta, by.x='variable', by.y='SampleID')

#summarize data
ko_vk_sum<-ddply(ko_vk_m, c('Sex', 'Diet_group', "KO"), summarize, mean=mean(value), sd=sd(value), n=length(value), se=sd/n)

#plot means
ggplot(ko_vk_sum, aes(KO, mean))+
geom_point()+
theme_bw()+
geom_errorbar(aes(ymax=mean+se, ymin=mean-se, width=0.2), stat="identity")+
facet_grid(Diet_group~Sex)+
coord_flip()

#write to table
write.table(ko_vk_sum, "picrust2_ko_vitK_sum.txt", quote=F, sep="\t", row.names=F)

#split by sex
vk_split<-split(ko_vk_m, ko_vk_m$Sex)
vk_male<-vk_split$M
vk_female<-vk_split$F

#split sex by KO
vk_male_split<-split(vk_male, vk_male$KO)
vk_female_split<-split(vk_female, vk_female$KO)

#t-test!
pairwise.t.test(vk_male_split$K01661$value, vk_male_split$K01661$Diet_group, p.adjust.method = 'hochberg')
pairwise.t.test(vk_male_split$K01851$value, vk_male_split$K01661$Diet_group, p.adjust.method = 'hochberg')
pairwise.t.test(vk_male_split$K01911$value, vk_male_split$K01661$Diet_group, p.adjust.method = 'hochberg')
pairwise.t.test(vk_male_split$K02361$value, vk_male_split$K01661$Diet_group, p.adjust.method = 'hochberg')
pairwise.t.test(vk_male_split$K02548$value, vk_male_split$K01661$Diet_group, p.adjust.method = 'hochberg')
pairwise.t.test(vk_male_split$K02549$value, vk_male_split$K01661$Diet_group, p.adjust.method = 'hochberg')
pairwise.t.test(vk_male_split$K02551$value, vk_male_split$K01661$Diet_group, p.adjust.method = 'hochberg')
pairwise.t.test(vk_male_split$K02552$value, vk_male_split$K01661$Diet_group, p.adjust.method = 'hochberg')
pairwise.t.test(vk_male_split$K03183$value, vk_male_split$K01661$Diet_group, p.adjust.method = 'hochberg')
pairwise.t.test(vk_male_split$K08680$value, vk_male_split$K01661$Diet_group, p.adjust.method = 'hochberg')
pairwise.t.test(vk_male_split$K19222$value, vk_male_split$K01661$Diet_group, p.adjust.method = 'hochberg')
#nothing significant

pairwise.t.test(vk_female_split$K01661$value, vk_female_split$K01661$Diet_group, p.adjust.method = 'hochberg')
pairwise.t.test(vk_female_split$K01851$value, vk_female_split$K01661$Diet_group, p.adjust.method = 'hochberg')
pairwise.t.test(vk_female_split$K01911$value, vk_female_split$K01661$Diet_group, p.adjust.method = 'hochberg')
pairwise.t.test(vk_female_split$K02361$value, vk_female_split$K01661$Diet_group, p.adjust.method = 'hochberg')
pairwise.t.test(vk_female_split$K02548$value, vk_female_split$K01661$Diet_group, p.adjust.method = 'hochberg')
pairwise.t.test(vk_female_split$K02549$value, vk_female_split$K01661$Diet_group, p.adjust.method = 'hochberg')
pairwise.t.test(vk_female_split$K02551$value, vk_female_split$K01661$Diet_group, p.adjust.method = 'hochberg')
pairwise.t.test(vk_female_split$K02552$value, vk_female_split$K01661$Diet_group, p.adjust.method = 'hochberg')
pairwise.t.test(vk_female_split$K03183$value, vk_female_split$K01661$Diet_group, p.adjust.method = 'hochberg')
pairwise.t.test(vk_female_split$K08680$value, vk_female_split$K01661$Diet_group, p.adjust.method = 'hochberg')
pairwise.t.test(vk_female_split$K19222$value, vk_female_split$K01661$Diet_group, p.adjust.method = 'hochberg')
#nothing significant

```
