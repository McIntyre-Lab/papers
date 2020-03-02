
libname tappas "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata/tappas";




/* for regulated transcript, what is the overlap between the 5 genotypes? */

data transcript_cnts2 ;
set tappas.tappas_results_transcripts ;
where regulated_B73  ne "" or regulated_Mo17  ne "" or regulated_C123  ne "" or regulated_Hp301  ne "" or regulated_NC338  ne "" ; 
keep name_description regulated_: ;
run;

data transcript_cnts ;
set transcript_cnts2 ;
if regulated_B73 ne ""  then flag_on_B73 = 1 ;     else flag_on_B73 = 0 ;
if regulated_Mo17 ne "" then flag_on_Mo17 = 1 ;   else flag_on_Mo17 = 0 ;
if regulated_C123 ne "" then flag_on_C123 = 1 ;   else flag_on_C123 = 0 ;
if regulated_Hp301 ne "" then flag_on_Hp301 = 1 ; else flag_on_Hp301 = 0 ;
if regulated_NC338 ne "" then flag_on_NC338 = 1 ; else flag_on_NC338 = 0 ;
counts = (flag_on_B73 + flag_on_Mo17 + flag_on_C123 + flag_on_Hp301 + flag_on_NC338) ;
drop regulated_:;
run ;


proc freq data = transcript_cnts ;
tables flag_on_B73 * flag_on_C123 * flag_on_Mo17 * flag_on_Hp301 * flag_on_NC338 / out = transcript_venn ;
run;
    
data transcript_venn_cnts ;
set transcript_venn ;
drop percent ;
sum_on = (flag_on_B73 + flag_on_C123 + flag_on_Mo17 + flag_on_Hp301 + flag_on_NC338) ;
run;

proc sort data = transcript_venn_cnts ;
by sum_on ;
run;

proc means data = transcript_venn_cnts sum ;
var count ;
by sum_on ;
output out = geno_sums sum = sum ;
run;
    /*  regulated in 1 genotype = 3603
                     2 genotype = 1643
                     3 genotype = 1123
                     4 genotype = 1006
                     5 genotype = 81
                    total      7456 transcript at 5% FDR     
*/  
                              
