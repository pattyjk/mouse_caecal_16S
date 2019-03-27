###  Cluster Zero Radius OTUs (ZOTUs)
```
#use dereplicated reads
#need the newest version of usearch (11.0.667, 32 bit cause it free)
usearch32 -unoise3 mergedfastq/uniques_combined_merged.fastq -zotus zotus.fa -tabbedout zotus.txt


#count ZOTUS
grep -c '>' zotus.fa
#1214

```

### Map reads at 97% to ZOTUs
```
usearch64 -otutab mergedfastq/merged.fq -otus zotus.fa -uc zotus_mapped.uc -otutabout zotus_table.txt -notmatchedfq zotus_not_matched.fq
```

### Assign taxonomy to zOTUs
```
~/usearch64 -sintax UNOISE_files/zotus.fa -db ~/SILVA_132_SINTAX/SILVA_132-90_SINTAX.udb --tabbedout zotus_sintax_annotations.txt -strand both -sintax_cutoff 0.8
```
