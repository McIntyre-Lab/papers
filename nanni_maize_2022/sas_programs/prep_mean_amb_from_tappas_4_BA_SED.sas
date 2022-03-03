libname tappas "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata/tappas";
libname pacbio "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata";

/*
create sbys mean ambient for BA and SED 
*/



%macro prep (geno) ;

data transcript_&geno ;
set tappas.transcript_tappas_&geno ;
keep ambient_meanexplevel_&geno transcript;
rename ambient_meanexplevel_&geno = &geno.;
run;

proc sort data = transcript_&geno ;
by transcript ;
run;

%mend ;

%prep (B73) ;
%prep (C123) ;
%prep (Hp301) ;
%prep (Mo17) ;
%prep (NC338) ;

data transcript_tappas_amb_mean ;
merge transcript_: ;
by transcript ;
run ;

data tappas.transcript_tappas_amb_mean ;
set transcript_tappas_amb_mean ;
run ;

/* create design file for BA and SED */
proc transpose data = tappas.transcript_tappas_amb_mean out = tall ;
by transcript ;
run ;

data df ;
set tall ;
label _name_ = "sampleID";
keep _name_ ;
rename _name_ = sampleID ;
run ;

proc sort data = df nodups ;
by _all_ ;
run;

data tappas.df_genotypes ;
set df ;
run ;

/* export for BA and SED */
proc export data = tappas.df_genotypes 
outfile = "!MCLAB/maize_ozone_FINAL/2018/PacBio/text_data/df_genotypes.tsv"
dbms = tab replace ;
run ;

proc export data = tappas.transcript_tappas_amb_mean 
outfile = "!MCLAB/maize_ozone_FINAL/2018/PacBio/text_data/transcript_tappas_amb_mean.tsv"
dbms = tab replace ;
run ;







