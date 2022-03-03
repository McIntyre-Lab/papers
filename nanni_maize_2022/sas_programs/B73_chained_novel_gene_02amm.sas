


/*

B73 chained classification file from sqanti filtered:  should have 1531 novel genes
    associated_gene column search for "novelGene_"

*/

proc import datafile = "!MCLAB/maize_ozone_FINAL/2018/PacBio/sqanti_post_filter_b73/SQANTI_classification.txt"
out = classFile 
dbms = tab replace ;
guessingrows = MAX ;
run;

proc contents data = classFile ; run;

data classFile2 ;
set classFile ;
if find(associated_gene, "novelGene_") ge 1 then flag_novel_gene = 1 ;
else flag_novel_gene = 0;
keep associated_gene flag_novel_gene ;
run ;
 
proc sort data = classFile2 nodups ;
by _all_;
run;

proc freq data = classFile2 ;
tables flag_novel_gene ;
run;

data list_novel_genes_b73_chain ;
set classFile2 ;
where flag_novel_gene = 1 ;
keep associated_gene ;
run ;

proc export data = list_novel_genes_b73_chain 
outfile = "!MCLAB/maize_ozone_FINAL/pacbio_paper/penultimate_version/draft_figs_tables/list_novel_genes_B73_chained_filtered.tsv"
dbms = tab replace ;
run ;



