libname tappas "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata/tappas";
libname pacbio "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata";

/*
identify transcripts that are highly discordant between genotypes using BA flags
*/

proc import datafile = "!MCLAB/maize_ozone_FINAL/2018/PacBio/tappas_results/BA_feature_flags_mean_Amb.tsv"
out = ba 
dbms = tab replace ;
run;

proc contents data = ba ; run;

proc freq data = BA noprint;
tables flag_feature_BA_cooks * flag_feature_BA_dffits * flag_feature_BA_pearson/ out = ba_flags ;
run;

proc freq data = BA noprint;
tables flag_feature_BA_outlier/ out = ba_flags_outlier ;
run;

ods pdf file = "!MCLAB/maize_ozone_FINAL/2018/PacBio/tappas_results/table_BA_flag_cnts_mean_Amb.pdf" ;

title "mean ambient transcript counts discordant across genotypes" ;
proc print data = BA_flags ; run ;
proc print data = ba_flags_outlier ; run;

ods pdf close ;

