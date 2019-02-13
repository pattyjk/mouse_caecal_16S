## Processing of 16S rRNA gene amplicons

### Merge Paired End Reads
```
#decompress the reads
gunzip *.gz

mkdir mergedfastq

./usearch64 -fastq_mergepairs *R1*.fastq -relabel @ -fastq_maxdiffs 10 -fastqout mergedfastq/merged.fq -fastq_merge_maxee 1.0 -fastq_minmergelen 250 -fastq_maxmergelen 300

   8443547  Pairs (8.4M)
   7172229  Merged (7.2M, 84.94%)
   3450300  Alignments with zero diffs (40.86%)
   1255657  Too many diffs (> 10) (14.87%)
     14763  No alignment found (0.17%)
         0  Alignment too short (< 16) (0.00%)
       845  Merged too short (< 250)
        53  Merged too long (> 300)
         0  Exp.errs. too high (max=1.0) (0.00%)
      1117  Staggered pairs (0.01%) merged & trimmed
    247.16  Mean alignment length
    252.84  Mean merged length
      0.41  Mean fwd expected errors
      0.78  Mean rev expected errors
     0.12  Mean merged expected errors
```

## Dereplicate sequences
```
./usearch64 -fastx_uniques mergedfastq/merged.fq -fastqout mergedfastq/uniques_combined_merged.fastq -sizeout
```

## Remove Singeltons
```
./usearch64 -sortbysize mergedfastq/uniques_combined_merged.fastq -fastqout mergedfastq/nosigs_uniques_combined_merged.fastq -minsize 2
```

## Precluster Sequences
```
./usearch64 -cluster_fast mergedfastq/nosigs_uniques_combined_merged.fastq -centroids_fastq mergedfastq/denoised_nosigs_uniques_combined_merged.fastq -id 0.9 -maxdiffs 5 -abskew 10 -sizein -sizeout -sort size
```

## Reference-based OTU picking against the Silva v. 1.28
```
#pull Silva and extract it
# wget https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_128_release.tgz

./usearch64 -usearch_global mergedfastq/denoised_nosigs_uniques_combined_merged.fastq -id 0.97 -db ./Silva_128_release/SILVA_128_QIIME_release/rep_set/rep_set_16S_only/97/97_otus_16S.fasta  -strand plus -uc mergedfastq/ref_seqs.uc -dbmatched mergedfastq/closed_reference.fasta -notmatchedfq mergedfastq/failed_closed.fq
```

## Sort by size and then de novo OTU picking on sequences that failed to hit GreenGenes
```
./usearch64 -sortbysize mergedfastq/failed_closed.fq -fastaout mergedfastq/sorted_failed_closed.fq

./usearch64 -cluster_otus mergedfastq/sorted_failed_closed.fq -minsize 2 -otus mergedfastq/denovo_otus.fasta -relabel OTU_dn_ -uparseout denovo_out.up
```

## Combine the rep sets between de novo and reference-based OTU picking
```
cat mergedfastq/closed_reference.fasta mergedfastq/denovo_otus.fasta > full_rep_set.fna
```

## Map rep_set back to pre-dereplicated sequences and make OTU tables
```
./usearch64 -usearch_global mergedfastq/merged.fq -db full_rep_set.fna  -strand plus -id 0.97 -uc OTU_map.uc -otutabout OTU_table.txt
```

## Assign taxonomy with SILVA 
```
assign_taxonomy.py -i full_rep_set.fna -o taxonomy -r Silva_128_release/SILVA_128_QIIME_release/rep_set/rep_set_16S_only/97/97_otus_16S.fasta -t '/media/pattyjk/TOSHIBA EXT/MiSeq Runs/Wolfe_Mouse_cecal_16S_amplicons/Silva_128_release/SILVA_128_QIIME_release/taxonomy/16S_only/97/consensus_taxonomy_all_levels.txt'

#add tax to OTU table
biom convert -i OTU_table.txt -o OTU_table.biom --table-type='OTU table' --to-hdf5

biom add-metadata -i OTU_table.biom -o otu_table_tax.biom --observation-metadata-fp=taxonomy/full_rep_set_tax_assignments.txt --sc-separated=taxonomy --observation-header=OTUID,taxonomy

#summarize OTU table
biom summarize-table -i otu_table_tax.biom -o otu_table_summmary.txt
```


