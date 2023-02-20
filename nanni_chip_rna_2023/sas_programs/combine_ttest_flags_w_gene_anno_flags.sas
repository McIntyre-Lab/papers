
libname chiprna '!MCLAB/Dros_PB_ChIP/sasdata/chipRnaMs';

libname dTemp '/home/ammorse/TB14/maize_ainsworth/models';

/* assemble ttest flags with species_gene_flags_anno

    inputs: 
        dTemp.mel_manual_model_check  from        mel_frags_ttest_4_de_sex_assemble_flags.sas
        chiprna.mel_gene_flags_anno   from        gene_count_14amm.sas
            

** chiprna.mel_gene_flags_anno ==> this file featureType = fragments and flag_multigene = 0  
*/

%macro flagging (species) ;

data &species._flags ;
set dTemp.&species._manual_model_check ;
where featureType = "fragment" ;
drop featureType ;
run ;  

data &species._anno ;
set chiprna.&species._gene_flags_anno ;  /* this file featureType = fragments and flag_multigene = 0 */
run;

proc sort data = &species._flags ;
by FBgn;
proc sort data = &species._anno ;
by FBgn ;
run; 

data &species._ttest_anno  oops;
merge &species._flags (in=in1) &species._anno  (in=in2) ;
by FBgn ;
if in2 then output &species._ttest_anno ;
else output oops ;  /* 0 in oops */
run ; 

data dtemp.&species._ttest_flags_with_anno ;
set &species._ttest_anno ;
run;

proc export data = dtemp.&species._ttest_flags_with_anno 
outfile = "!MCLAB/Dros_PB_ChIP/RNAseq/model_output/&species._ttest_flags_with_anno.csv"
dbms = csv replace ;
run;

data chiprna.&species._ttest_flags_with_anno ;
set dtemp.&species._ttest_flags_with_anno ;
run ;

%mend ;

%flagging (mel) ;
%flagging (sim) ;

/*  merge kept if in ttest flags only 

  */

