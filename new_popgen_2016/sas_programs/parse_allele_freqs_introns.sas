/* I need to parse the allele frequency output from vcftools so that I can plot the distribution*/

libname dsim "!MCLAB/ethanol/Sim_Pop_Gen/sas_data";
libname fnew "/home/fnew/dsim";

*Import the allele freq file, be careful of the fifth column;

data dsim.allele_freqs_intron;
  set fnew.allele_freq_intron;
  run;

*Set missing values to 0;
data dsim.allele_freqs_intron;
  set dsim.allele_freqs_intron;
  array change _character_;
    do over change;
    if change=' ' then change=0;
    end;
  run;

*How many bi-, tri-allelic loci are there? ;

*Create flags!;
data allele_freq_intron_2;
  set dsim.allele_freqs_intron;
  if N_alleles = 1 then flag_mono=1; else flag_mono=0;
  if N_alleles = 2 then flag_bi=1; else flag_bi=0;
  if N_alleles = 3 then flag_tri=1; else flag_tri=0;
  if N_alleles = 4 then flag_tetra=1; else flag_tetra=0;
  run;

proc freq data=allele_freq_intron_2;
  tables flag_:;
  run;

/* mono:    785,227     96.21%
   bi:      30,843      3.78%
   tri:     125         0.02%
   total:   816195      100%
*/

*Drop all info except allele_freqs;
data mono;
  set allele_freq_intron_2;
  if flag_mono=1;
  run;

data di;
  set allele_freq_intron_2;
  if flag_bi=1;
  run;

data tri;
  set allele_freq_intron_2;
  if flag_tri=1;
  run;


*Parse allele freqs;
data di2;
  set di;
  keep allele_freq allele_freq_2 ;
  run;

data tri2;
  set tri;
  keep allele_freq allele_freq_2 allele_freq_3;
  run;



data di3;
  set di2;
  af1=substr(allele_freq,3);
  af2=substr(allele_freq_2,3);
  keep af1 af2;
  run;

data tri3;
  set tri2;
  af1=substr(allele_freq,3);
  af2=substr(allele_freq_2,3);
  af3=substr(allele_freq_3,3);
  keep af1 af2 af3;
  run;



*I need the minor allele frequencies, get the minimum value in each row;
data minDI;
  set di3;
  array x {*} _char_;
  min=min(of x[*]);
  run;


data minTRI;
  set tri3;
  array x {*} _char_;
  min=min(of x[*]);
  run;



*Export the allele frequencies;
proc export data = minDI
	outfile="!MCLAB/ethanol/Sim_Pop_Gen/output/biallele_frequencies_intron.csv"
	dbms=CSV REPLACE;
	run;


proc export data = minTRI
	outfile="!MCLAB/ethanol/Sim_Pop_Gen/output/triallele_frequencies_intron.csv"
	dbms=CSV REPLACE;
	run;



 data WORK.Rogers_gff    ;
      %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
      infile '/home/fnew/mclab/useful_dsim_data/rogers_2014/GFF/dsim_update.final.gff'
 delimiter='09'x MISSOVER DSD lrecl=32767 ;
         informat VAR1 best32. ;
         informat VAR2 $7. ;
         informat VAR3 $30. ;
         informat VAR4 best32. ;
         informat VAR5 best32. ;
         informat VAR6 best32. ;
         informat VAR7 $1. ;
         informat VAR8 best32. ;
         informat VAR9 $16. ;
         format VAR1 best12. ;
         format VAR2 $7. ;
         format VAR3 $30. ;
         format VAR4 best12. ;
         format VAR5 best12. ;
         format VAR6 best12. ;
         format VAR7 $1. ;
         format VAR8 best12. ;
         format VAR9 $16. ;
      input
                  VAR1
                  VAR2 $
                  VAR3 $
                  VAR4
                  VAR5
                  VAR6
                  VAR7 $
                  VAR8
                  VAR9 $
      ;
      if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
      run;

data rogers_feat;
  set rogers_gff;
  keep var3;
  run;

  proc sort data=rogers_feat nodupkey;
    by var3;
    run;
