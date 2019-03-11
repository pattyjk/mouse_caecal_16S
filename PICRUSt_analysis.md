## Generating predicted metagenomes wit PICRUSt v. 1

```
#download GreenGrenes 13.5
http://qiime.org/home_static/dataFiles.html
```

## Re-do OTU picking against GreenGenes 13.5
```
#closed reference OTU pick aainst GreenGenes
./usearch64 -usearch_global mergedfastq/denoised_nosigs_uniques_combined_merged.fastq -id 0.97 -strand plus -uc mergedfastq/picrust_ref_seqs.uc -dbmatched /mergedfastq/picrust_closed_reference.fasta -db gg_13_5_otus/rep_set/97_otus.fasta

#map OTUs back to raw reads
./usearch64 -usearch_global mergedfastq/merged.fq'-db mergedfastq/picrust_closed_reference.fasta -strand plus -id 0.97 -uc mergedfastq/picrust_OTU_map.uc -otutabout picrust_OTU_table.txt

#convert to JSON biom
biom convert -i picrust_OTU_table.txt'-o picrust_OTU_table.biom --to-json --table-type="OTU table"

#analysis on Galaxy
http://huttenhower.sph.harvard.edu/galaxy/

#including copy number correction
```
