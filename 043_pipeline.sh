#!/bin/bash

#043 Seagrass amplicon metagenomics analysis pipeline and command lines
#T. R. Allnutt, 2017

#merge and filter reads
usearch -fastq_mergepairs ~/d/043_seagrass/reads/*_R1.fastq -fastqout merged.fastq -relabel @
usearch -fastq_mergepairs ~/d/043_seagrass/reads/*_R1.fastq -fastqout merged2.fastq -relabel @

usearch -fastq_filter merged.fastq -fastq_truncqual 18 -fastq_maxee 1.0 -fastq_minlen 170 -fastaout filtered.fasta
usearch -fastq_filter merged2.fastq -fastq_truncqual 18 -fastq_maxee 1.0 -fastq_minlen 170 -fastaout filtered2.fasta

#replace dashes in sample names
./rename.sh
./rename2.sh

#UPARSE 16S
usearch -derep_fulllength filtered.fasta -sizeout -fastaout uniques.fa
usearch -cluster_otus uniques.fa -minsize 2 -otus otus.fa -relabel Otu
usearch -utax otus.fa -db ~/db/utax/16s.udb -strand both -strand both -fastaout otus_tax.fa -utax_cutoff 0.9 -log utax.log
usearch -usearch_global otus.fa -db ~/db/utax/16s.udb -strand both -id 0.97 -alnout otus_ref.aln -userout otus_ref.user -userfields query+target+id
usearch -usearch_global merged.fastq -db otus_tax.fa -strand both -id 0.97 -log make_otutab.log -otutabout otutab.txt -biomout otutab.biom
biom convert -i otutab.biom -o otutab.txt --to-tsv --header-key taxonomy && sed -i 's/; //g' otutab.txt && sed -i 's/"//g' otutab.txt

#UPARSE ITS
#make Unite db
usearch -makeudb_utax ~/db/utax/refdb.fa -output its.udb -taxconfsin ~/db/utax/utaxref/unite_v7/taxconfs/full_length.tc
#n.b. in practice the file suffix '2' was not used because analysis was carried out in separate directories
search -derep_fulllength filtered2.fasta -sizeout -fastaout uniques2.fa
usearch -cluster_otus uniques2.fa -minsize 2 -otus otus2.fa -relabel Otu
usearch -utax otus2.fa -db ~/db/utax/its.udb -strand both -strand both -fastaout otus_tax2.fa -utax_cutoff 0.9 -log utax.log
usearch -usearch_global otus2.fa -db ~/db/utax/its.udb -strand both -id 0.97 -alnout otus_ref2.aln -userout otus_ref2.user -userfields query+target+id
usearch -usearch_global merged2.fastq -db otus_tax2.fa -strand both -id 0.97 -log make_otutab.log -otutabout otutab2.txt -biomout otutab2.biom
biom convert -i otutab2.biom -o otutab2.txt --to-tsv --header-key taxonomy && sed -i 's/; //g' otutab2.txt && sed -i 's/"//g' otutab2.txt

#Community analysis using QIIME - idenitcal for 16S and ITS
#filter
filter_otus_from_otu_table.py -i merged2otutab.biom -o merged2otutab2.biom -e exclude.txt
filter_otus_from_otu_table.py -i merged2otutab2.biom -o merged2otufilt.biom --min_count_fraction 0.001 -s 2
normalize_table.py -i merged2otufilt.biom -o merged2otunorm.biom
alpha_rarefaction.py -i merged2otufilt.biom -m merged2mapping.txt -o rarefaction2 -n 25 -f -p qiime_parameters1.txt
#compare alpha
compare_alpha_diversity.py -i rarefaction2/alpha_div_collated/simpson.txt -m mapping.txt -c treatment,site -o merged2compare_simpson -t nonparametric -d 30079
compare_alpha_diversity.py -i rarefaction2/alpha_div_collated/shannon.txt -m mapping.txt -c treatment,site -o merged2compare_shannon -t nonparametric  -d 30079
compare_alpha_diversity.py -i rarefaction2/alpha_div_collated/observed_species.txt -m mapping.txt -c treatment,site -o merged2compare_obs -t nonparametric -d 30079
#test groups
group_significance.py -i merged2otunorm.biom -m merged2mapping.txt -s ANOVA -c site -o merged2site-anova.txt
group_significance.py -i merged2otunorm.biom -m merged2mapping.txt -s ANOVA -c treatment -o merged2treatment-anova.txt
group_significance.py -i merged2otunorm.biom -m merged2mapping.txt -s ANOVA -c merged2salinity_descript -o salinity_desc-anova.txt
group_significance.py -i merged2otunorm.biom -m merged2mapping.txt -s ANOVA -c merged2salinity_category -o salinity_cat-anova.txt
group_significance.py -i merged2otunorm.biom -m merged2mapping.txt -s ANOVA -c merged2leaf_length_cat -o merged2leaf_len_cat-anova.txt

#Beta diversity analsyis
#get list of OTUs
cut -f1 otufilt.txt | sed -n '1!p' > otu.list
get-seqs-from-file.py otu.list otus.fa otulist.fa fasta fasta
#align and tree
muscle -in otulist.fa -out otuslist.aln
FastTree -nt otuslist.aln > otuslist.tree
#n.b. use appropriate mapping file
sort_otu_table.py -i otunorm.biom -o otusort.biom -m [16S/ITS]mapping.txt -s Description
#calculate unifrac matrix
beta_diversity.py -i otusort.biom -m weighted_unifrac -o ./ -t otuslist.tree
#PCO
principal_coordinates.py -i weighted_unifrac_otusort.txt -o pco.out

#PICRUSt functional analysis
#replace uparse format headers with '_'
sed -i 's/\./_/g' filtered-clipped.fasta

#QIIME otu picking
pick_closed_reference_otus.py -i filtered-clipped.fasta -o output -r ~/db/gg_13_8_otus/rep_set/97_otus.fasta -t ~/db/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt -f  -p params.txt
biom convert -i output/otu_table.biom -o otu_table.txt --to-tsv

#convert gg codes to taxonomy in order to remove green otus:
gg2tax-chl.py otu_table.txt ~/db/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt otu_tax.txt chlorop,cyan,plant,alga
cut -f1 chl.hits >exclude.txt
#manually fixed gg format taxonomy to uparse format.
biom convert -i otu_tax.txt -o otu_tax.biom --to-json 
filter_otus_from_otu_table.py -i otu_tax.biom -o otufilt.biom -e exclude.txt
filter_otus_from_otu_table.py -i otufilt.biom -o otufilt2.biom --min_count_fraction 0.001 -s 2
biom convert -i otufilt2.biom -o otufilt2.txt --to-tsv --header-key taxonomy && sed -i 's/; //g' otufilt2.txt && sed -i 's/"//g' otufilt2.txt

#picrust
normalize_by_copy_number.py -i otufilt2.biom -o otunorm.c.biom
predict_metagenomes.py -i otunorm.c.biom -o picrust.biom --with_confidence
biom convert -i picrust.biom -o picrust.txt --to-tsv
categorize_by_function.py -i picrust.biom -c KEGG_Pathways -l 3 -o kegg.biom

#all levels
categorize_by_function.py -i picrust.biom -c KEGG_Pathways -l 1 -o kegg-1.biom
categorize_by_function.py -i picrust.biom -c KEGG_Pathways -l 2 -o kegg-2.biom &
categorize_by_function.py -i picrust.biom -c KEGG_Pathways -l 3 -o kegg-3.biom

biom convert -i kegg-1.biom -o kegg-1.txt --to-tsv 
biom convert -i kegg-2.biom -o kegg-2.txt --to-tsv
biom convert -i kegg-3.biom -o kegg-3.txt --to-tsv
#change samples names to underscore, not '-'
biom convert -i kegg-1.txt -o kegg-1.biom --to-json
biom convert -i kegg-2.txt -o kegg-2.biom --to-json
biom convert -i kegg-3.txt -o kegg-3.biom --to-json

biom convert -i kegg.txt -o kegg.biom --to-json

#test significance
group_significance.py -i kegg.biom -m mapping.txt -s ANOVA -c site -o site-anova.txt
group_significance.py -i kegg.biom -m mapping.txt -s ANOVA -c treatment -o treatment-anova.txt
group_significance.py -i kegg.biom -m mapping.txt -s ANOVA -c salinity_descript -o salinity_descript-anova.txt
group_significance.py -i kegg.biom -m mapping.txt -s ANOVA -c salinity_category -o salinity_category.txt
group_significance.py -i kegg.biom -m mapping.txt -s ANOVA -c leaf_length_cat -o leaf_length_cat-anova.txt




