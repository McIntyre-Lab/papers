
libname make "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/make_combination_flag_file/sasdata";


/*
flag genes after sqanti filtering

gtf with sqanti post filtered  genes: /McIntyre_Lab/maize_ozone_FINAL/2018/PacBio/sqanti_post_filter_b73/sqanti_b73_filtered_corrected_associated_gene.gtf

        14,877 genes with assembled transcript (flag_assembled_transcript = 1)
        39 assembled genes with fusion transcripts (flag_assembled_fusion_transcript = 1)
        1531 assembled genes with novel transcripts  ( flag_assembled_novel_transcript = 1)
        fusion and novel = 1531 + 39 = 1570

Gtf with all genes:
/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/useful_maize_info/Zea_mays.B73_RefGen_v4.41.gtf

    46,430 unique genes 

gtf fields:
    seqname
    source
    feature
    start 
    end
    score
    strand
    frame
    attribute



fusion annotations = mclab/SHARE/McIntyre_Lab/useful_maize_info/FSM_consolidation_maize_B73_EA_150bp/FSM_consolidation_maize_B73_fusions.tsv

*/
/* include flag for multigene - note we do NOT do coverage counts on multigenes, rRNA, tRNA etc, */
/* list fusions 2 genes */
proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/useful_maize_info/FSM_consolidation_maize_B73_EA_150bp/FSM_consolidation_maize_B73_fusions.tsv"
out = anno 
dbms = tab replace ;
guessingrows = MAX ;
run; /* 402,691 rows  */


data flag_multigene;
set anno ;
keep primary_FBgn flag_multigene ;
rename primary_FBgn = gene ;
run;  /* 48033 obs  */

proc sort data = flag_multigene nodups ;
by _all_ ;
run;

proc sort data = flag_multigene ;
by gene descending flag_multigene;
run;

proc sort data = flag_multigene out =flag_multigene_nodup nodupkey ;
by gene ;
run;  /* 46022 obs */

proc freq data = flag_multigene_nodup ;
tables flag_multigene ;
run;
/*flag_multigene    Frequency
---------------------------
             0       41512
             1        4510  */

proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio/sqanti_post_filter_b73/sqanti_b73_filtered_corrected_associated_gene.gtf"
out = gtf 
dbms = tab replace ;
getnames=no;
guessingrows = MAX ;
run;

data sqanti ;
retain gene gene2 gene3 ;
set gtf ;
gene = scan(var9, 2, ';') ;
gene2 = tranwrd(gene, "gene_id ", "");
gene3 = compress(tranwrd(gene2, '"', ''));
run ;

data sqanti2 ;
set sqanti ;
flag_assembled_transcript = 1 ;
if find(gene3, "novelGene") ge 1 then flag_assembled_novel_transcript = 1 ;
if find(gene3, "novelGene") = 0 and find(gene3, "_") ge 1 then do ;
    flag_assembled_fusion_transcript = 1 ;
end ;
keep gene3 flag_assembled_transcript flag_assembled_fusion_transcript flag_assembled_novel_transcript ;
run ;

data find ;
set sqanti2 ;
if find(gene3, "_") ge 1 ;
run;

proc sort data = sqanti2 nodups ;
by _ALL_ ;
run ;  /*   326870 starting obs
            14877 unique genes */

data sqanti_gene ;
retain gene3 flag_assembled_transcript ;
set sqanti2 ;
rename gene3 = gene ;
if flag_assembled_novel_transcript ne 1 then flag_assembled_novel_transcript = 0 ;
if flag_assembled_transcript ne 1 then flag_assembled_transcript = 0 ;
if flag_assembled_fusion_transcript ne 1 then flag_assembled_fusion_transcript = 0 ;
run;
    /*  14,877 genes with assembled transcript
        39 assembled genes with fusion transcripts
        1531 assembled genes with novel transcripts    */  


proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/useful_maize_info/Zea_mays.B73_RefGen_v4.41.gtf"
out = gtf_all
dbms = tab replace ;
datarow=6;
getnames=no;
guessingrows = MAX ;
run;


data gtf_all2 ;
retain gene gene2 gene3 ;
set gtf_all ;
gene = scan(var9, 1, ';') ;
gene2 = tranwrd(gene, "gene_id ", "");
gene3 = compress(tranwrd(gene2, '"', ''));
run ;  

data gtf_all3 ;
set gtf_all2 ;
keep gene3 ;
run ;

proc sort data = gtf_all3 nodups ;
by _ALL_ ;
run ;   /*  3083980 starting obs
            46430 unique genes */

data gtf_all_gene ;
retain gene3  ;
set gtf_all3 ;
flag_b73_reference_gene = 1 ;
rename gene3 = gene ;
run;

proc sort data = sqanti_gene ;
by gene ;
proc sort data = gtf_all_gene ;
by gene ;
proc sort data = flag_multigene_nodup ;
by gene ;
run ;

data flag_assembled ;
merge sqanti_gene (in=in1) gtf_all_gene (in=in2) ;
by gene ;
run ;  /* 48000 */


data flag_assembled2 ;
set flag_assembled ;
if flag_b73_reference_gene ne 1 then flag_b73_reference_gene = 0;
if flag_assembled_transcript ne 1 then flag_assembled_transcript = 0 ;
if flag_assembled_fusion_transcript ne 1 then flag_assembled_fusion_transcript = 0 ;
if flag_assembled_novel_transcript ne 1 then flag_assembled_novel_transcript = 0;
run ;

proc freq data = flag_assembled2  ;
tables flag_b73_reference_gene * flag_assembled_transcript * flag_assembled_fusion_transcript * flag_assembled_novel_transcript   / out = freq ;
run;
proc print data = freq ; run;
    /*  
flag_b73_
eference_    flag_assembled_     flag_assembled_      flag_assembled_
  gene          transcript      fusion_transcript    novel_transcript    COUNT

    0               1                   0                    1            1531
    0               1                   1                    0              39
    1               0                   0                    0           33123
    1               1                   0                    0           13307

*/

data flag_assembled3 ;
merge flag_assembled2 (in=in1) flag_multigene_nodup (in=in2) ;
by gene ;
run; 

data flag_assembled4 ;
set flag_assembled3 ;
if flag_multigene = . or flag_multigene = 1 then flag_in_cvrg_cnts_bed = 0 ;
else flag_in_cvrg_cnts_bed = 1;
run ;


proc freq data = flag_assembled4 ;
tables flag_multigene * flag_in_cvrg_cnts_bed * flag_assembled_novel_transcript * flag_assembled_fusion_transcript / out = fr_chk;
run;  
proc print data = fr_chk ; run;
/*
               flag_in_
    flag_    cvrg_cnts_     flag_assembled_     flag_assembled_
multigene        bed       novel_transcript    fusion_transcript    COUNT

        .         0                0                   0              408
        .         0                0                   1               39
        .         0                1                   0             1531
        0         1                0                   0            41512
        1         0                0                   0             4510  total = 48,000
*/

data make.flag_assembled_transcripts ;
set flag_assembled4 ;
run;

proc export data = make.flag_assembled_transcripts  
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/make_combination_flag_file/flag_assembled_transcripts.csv"
dbms = csv replace ;
run;






