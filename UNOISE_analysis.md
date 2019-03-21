

###  Cluster Zero Radius OTUs (ZOTUs)
```
#use dereplicated reads
#need the newest version of usearch (11.0.667, 32 bit cause it free)
usearch64 -unoise3  -zotus zotus.fa -tabbedout zotus_report.txt
```

### Map reads at 97% to ZOTUs
```
~/usearch64 -otutab mergedfastq/merged.fq -otus zotus.fa -uc zotus_mapped.uc -otutabout zotus_table.txt -notmatchedfq zotus_not_matched.fq
```
