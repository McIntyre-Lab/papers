
libname chiprna 'Z:\SHARE\McIntyre_Lab\Dros_PB_ChIP\sasdata\chipRnaMs';

libname chiprna "!MCLAB/Dros_PB_ChIP/sasdata/chipRnaMs";


data mel_gene_flags_anno1;
length mel_geneID $11;
set dTemp.mel_ttest_flags_with_anno;
mel_geneID = strip(fbgn);
drop fbgn;
run;

/* checking*/
proc sort data= mel_gene_flags_anno1 nodups;
by mel_geneid;
run;

proc contents data=mel_gene_flags_anno1;
run;

proc contents data=mel_2_sim;  /* from import_flags_anno_genome_features8_02amm.sas */
run;


proc sort data=mel_gene_flags_anno1;
by mel_geneID;

proc sort data= mel_2_sim;
by mel_geneID;

data mel_ortho mel_ortho1 mel_gene oops;
merge mel_2_sim (in=in1) mel_gene_flags_anno1(in=in2);
by mel_geneID;
if in1 and in2 then  output mel_ortho;
else if in1 then output mel_ortho1;
else if in2 then output mel_gene;
else output oops;
run;

proc contents data = mel_ortho ; run;

proc contents data = dTemp.sim_ttest_flags_with_anno ;
run ;

data sim ;
set dTemp.sim_ttest_flags_with_anno ;
geneID = strip(fbgn);
run;

proc sql noprint ;
select cats(name, '=', 'sim_',name)
into :list separated by ' '
from dictionary.columns where libname = "WORK" and memname = "SIM";
quit ;

proc datasets  nolist ;
modify sim;
rename &list ;
quit;


data sim_gene_flags_anno1 ;
set sim ;

label sim_num_frags_ttest = "sim_num_frags_ttest";
 label sim_ttest_minpval = "sim_ttest_minpval";
label     sim_num_transcripts = "sim_num_transcripts" ;

drop
sim_PERCENT
sim_XCHR
sim_annotation_ID
sim_organism_abbreviation
sim_secondary_FBgn_s_
sim_secondary_annotation_ID_s_
sim_unb ;
run;





proc sort data=sim_gene_flags_anno1 nodups;
by sim_geneID;

proc sort data= mel_ortho nodups;
by sim_geneID;
run;

data mel_sim_ortho no_ortho sim_gene oops;
merge mel_ortho (in=in1) sim_gene_flags_anno1(in=in2);
by sim_geneID;
if in1 and in2 then  output mel_sim_ortho;
else if in1 then output no_ortho;
else if in2 then output sim_gene;
else output oops;
run;    /*
            mel_sim_ortho 14241 obs
            no_ortho  19 obs
            sim_gene  2154 obs   */

data chiprna.mel_sim_ortho;
set mel_sim_ortho;
run;


PROC EXPORT DATA= chiprna.mel_sim_ortho
            OUTFILE= "!MCLAB/Dros_PB_ChIP/CHIP_RNA_ms/results/data_files/mel_sim_ortho.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;
