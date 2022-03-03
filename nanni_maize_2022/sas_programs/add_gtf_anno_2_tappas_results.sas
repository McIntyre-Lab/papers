



/* 
create table with all differentially regulated genes for Lisa
include annotation 

tappas.tappas_results_genes 
tappas.tappas_results_transcripts
*/


/* genes */
data genes ;
set tappas.tappas_results_genes ;
run ;

proc import datafile = "/home/ammorse/TB14/maize_gtf_files/B73_v4_geneIDs_geneNames.txt"
out = anno 
dbms = tab replace ;
guessingrows = MAX ;
getnames = no ;
run;

data anno2 ;
set anno ;
rename var1 = gene ;
rename var2 = geneName ;
run ;

proc sort data = genes ;
by gene ;
proc sort data = anno2 ;
by gene ;
run;

data tappas_results_genes_wAnno ;
merge genes (in=in1) anno2 (in=in2) ;
by gene ;
if in1 then output tappas_results_genes_wAnno ;
run ;

data tappas.tappas_results_genes_wAnno ;
retain gene geneName ;
set tappas_results_genes_wAnno ;
run ;

/* transcripts */
data transcript ;
set tappas.tappas_results_transcripts ;
run;

proc sort data = transcript ;
by gene ;
proc sort data = anno2 ;
by gene ;
run ;

data tappas_results_transcripts_wAnno ;
merge transcript (in=in1) anno2 (in=in2) ;
by gene ;
if in1 then output tappas_results_transcripts_wAnno ;
run ;

data tappas.tappas_results_transcripts_wAnno ;
retain transcript name_description gene geneName  ;
set tappas_results_transcripts_wAnno ;
run ;

proc sort data = transcripts_sig_wanno ;
by geneName ;
run;






