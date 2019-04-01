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
stuff goes hurr
```
