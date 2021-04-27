libname tappas "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata/tappas";

/*
import tappas output  
DEA = Differential Expression Analysis
add genotype prefix to all columns except gene or transcript
merge on gene or transcript

*/


/* genes */
%macro genes (genotype) ;

proc import datafile = "!MCLAB/maize_ozone_FINAL/2018/PacBio/RNAseq_rsem_expression_subset_fsm_ism_nic_nnc/corrected_tappas_output/&genotype._tappAS_DEA_Genes.tsv"
out = gene_&genotype. 
dbms = tab replace ;
guessingrows = MAX ;
run;
    
proc sql noprint;
   select cats(name,'=',name,"_&genotype.")
          into :list
          separated by ' '
          from dictionary.columns
          where libname = 'WORK' and memname = "GENE_&genotype.";
quit;

proc datasets library = work nolist;
modify gene_&genotype.;
rename &list;
quit ;

data tappas.gene_tappas_&genotype. ;
set gene_&genotype. ;
rename _gene_&genotype. = gene;
rename Name_Description_&genotype. = Name_Description ;
rename _1___Probability__&genotype. = One_minus_Probability_&genotype.;
run ;

proc sort data = tappas.gene_tappas_&genotype. ;
by gene ;
run;

%mend ;

%genes (B73) ;
%genes (MO17) ;
%genes (HP301) ;
%genes (C123) ;
%genes (NC338) ;



/* transcripts */
%macro transcripts (genotype) ;

proc import datafile = "!MCLAB/maize_ozone_FINAL/2018/PacBio/RNAseq_rsem_expression_subset_fsm_ism_nic_nnc//corrected_tappas_output/&genotype._tappAS_DEA_Transcripts.tsv"
out = transcript_&genotype. 
dbms = tab replace ;
guessingrows = MAX ;
run;
    
proc sql noprint;
   select cats(name,'=',name,"_&genotype.")
          into :list
          separated by ' '
          from dictionary.columns
          where libname = 'WORK' and memname = "TRANSCRIPT_&genotype.";
quit;

proc datasets library = work nolist;
modify transcript_&genotype.;
rename &list;
quit ;

data tappas.transcript_tappas_&genotype. ;
set transcript_&genotype. ;
rename _Transcript_&genotype. = Transcript ;
rename Gene_&genotype. = Gene ;
rename Name_Description_&genotype. = Name_Description ;
rename Gene_Description_&genotype. = Gene_Description ;

rename _1___Probability__&genotype. = One_minus_Probability_&genotype.;
run ;

proc sort data = tappas.transcript_tappas_&genotype. ;
by transcript ;
run;

%mend ;

%transcripts (B73) ;
%transcripts (MO17) ;
%transcripts (HP301) ;
%transcripts (C123) ;
%transcripts (NC338) ;

/* combine all gene results together */

data tappas.tappas_results_genes ;
merge tappas.gene_tappas_B73 tappas.gene_tappas_C123 tappas.gene_tappas_Hp301 tappas.gene_tappas_Mo17 tappas.gene_tappas_NC338 ;
by gene ;
run ;

/* combine all transcript results together */

data tappas.tappas_results_transcripts ;
merge tappas.transcript_tappas_B73 tappas.transcript_tappas_C123 tappas.transcript_tappas_Hp301 tappas.transcript_tappas_Mo17 tappas.transcript_tappas_NC338 ;
by transcript ;
run ;

