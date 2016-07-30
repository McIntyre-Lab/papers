/* Import LD files from R output. Then add flags for bins for LD and distance.
 * Then proc freq to see counts between the flags. Output for plotting in R. */

libname dsim "!MCLAB/ethanol/Sim_Pop_Gen/sas_data";



*import chr4;
DATA  dsim.ld4 ;
INFILE  "/home/fnew/dsim/ld/ld4_dist.txt" 
     DSD 
     LRECL= 22 ;
INPUT
 CHR_A
 R2
 dist
;
RUN;


*Import chrom 2L;
PROC FORMAT;
value CHR_A 
     1 = "2L" 
;

DATA  dsim.ld2l ;
INFILE  "/home/fnew/dsim/ld/ld2l_dist.txt" 
     DSD 
     LRECL= 24 ;
INPUT
 CHR_A
 R2
 dist
;
FORMAT CHR_A CHR_A. ;
RUN;


*Import Chrom 2R;

PROC FORMAT;
value CHR_A 
     1 = "2R" 
;

DATA  dsim.ld2r ;
INFILE  "/home/fnew/dsim/ld/ld2r_dist.txt" 
     DSD 
     LRECL= 24 ;
INPUT
 CHR_A
 R2
 dist
;
FORMAT CHR_A CHR_A. ;
RUN;


*Import Chrom 3L;


PROC FORMAT;
value CHR_A 
     1 = "3L" 
;

DATA  dsim.ld3l ;
INFILE  "/home/fnew/dsim/ld/ld3l_dist.txt" 
     DSD 
     LRECL= 24 ;
INPUT
 CHR_A
 R2
 dist
;
FORMAT CHR_A CHR_A. ;
RUN;

*Import Chrom 3R;


PROC FORMAT;
value CHR_A 
     1 = "3R" 
;

DATA  dsim.ld3r ;
INFILE  "/home/fnew/dsim/ld/ld3r_dist.txt" 
     DSD 
     LRECL= 23 ;
INPUT
 CHR_A
 R2
 dist
;
FORMAT CHR_A CHR_A. ;
RUN;

*Import Chrom X;

DATA  dsim.ldx ;
INFILE  "/home/fnew/dsim/ld/ldx_dist.txt" 
     DSD 
     LRECL= 24 ;
INPUT
 CHR_A
 R2
 dist
;
RUN;


/* Bin the data by distance */
data ld4_bin;
  set ld4;
  if dist LE 250 then flag_bin=1;
   else if dist > 250 and dist LE 500 then flag_bin=2;
   else if dist > 500 and dist LE 1000 then flag_bin=3;
   else if dist > 1000 and dist LE 2000 then flag_bin=4;
   else if dist >2000 and dist LE 2500 then flag_bin=5;
  
  if dist > 2500 then delete;

  if R2 LE 0.1 then flag_ld =1;
   else if R2 > 0.1 and R2 LE 0.2 then flag_ld=2;
   else if R2 > 0.2 and R2 LE 0.3 then flag_ld=3;
   else if R2 > 0.3 and R2 LE 0.4 then flag_ld=4;
   else if R2 > 0.4 and R2 LE 0.5 then flag_ld=5;
   else if R2 > 0.5 and R2 LE 0.6 then flag_ld=6;
   else if R2 > 0.6 and R2 LE 0.7 then flag_ld=7;
   else if R2 > 0.7 and R2 LE 0.8 then flag_ld=8;
   else if R2 > 0.8 and R2 LE 0.9 then flag_ld=9;
   else if R2 > 0.9 and R2 LE   1 then flag_ld=10;

  run;

data ld2l_bin;
  set ld2l;
  if dist LE 250 then flag_bin=1;
   else if dist > 250 and dist LE 500 then flag_bin=2;
   else if dist > 500 and dist LE 1000 then flag_bin=3;
   else if dist > 1000 and dist LE 2000 then flag_bin=4;
   else if dist >2000 and dist LE 2500 then flag_bin=5;
  
  if dist > 2500 then delete;

  if R2 LE 0.1 then flag_ld =1;
   else if R2 > 0.1 and R2 LE 0.2 then flag_ld=2;
   else if R2 > 0.2 and R2 LE 0.3 then flag_ld=3;
   else if R2 > 0.3 and R2 LE 0.4 then flag_ld=4;
   else if R2 > 0.4 and R2 LE 0.5 then flag_ld=5;
   else if R2 > 0.5 and R2 LE 0.6 then flag_ld=6;
   else if R2 > 0.6 and R2 LE 0.7 then flag_ld=7;
   else if R2 > 0.7 and R2 LE 0.8 then flag_ld=8;
   else if R2 > 0.8 and R2 LE 0.9 then flag_ld=9;
   else if R2 > 0.9 and R2 LE   1 then flag_ld=10;

  run;

data ld2r_bin;
  set ld2r;
  if dist LE 250 then flag_bin=1;
   else if dist > 250 and dist LE 500 then flag_bin=2;
   else if dist > 500 and dist LE 1000 then flag_bin=3;
   else if dist > 1000 and dist LE 2000 then flag_bin=4;
   else if dist >2000 and dist LE 2500 then flag_bin=5;
  
  if dist > 2500 then delete;

  if R2 LE 0.1 then flag_ld =1;
   else if R2 > 0.1 and R2 LE 0.2 then flag_ld=2;
   else if R2 > 0.2 and R2 LE 0.3 then flag_ld=3;
   else if R2 > 0.3 and R2 LE 0.4 then flag_ld=4;
   else if R2 > 0.4 and R2 LE 0.5 then flag_ld=5;
   else if R2 > 0.5 and R2 LE 0.6 then flag_ld=6;
   else if R2 > 0.6 and R2 LE 0.7 then flag_ld=7;
   else if R2 > 0.7 and R2 LE 0.8 then flag_ld=8;
   else if R2 > 0.8 and R2 LE 0.9 then flag_ld=9;
   else if R2 > 0.9 and R2 LE   1 then flag_ld=10;

  run;

data ld3l_bin;
  set ld3l;
  if dist LE 250 then flag_bin=1;
   else if dist > 250 and dist LE 500 then flag_bin=2;
   else if dist > 500 and dist LE 1000 then flag_bin=3;
   else if dist > 1000 and dist LE 2000 then flag_bin=4;
   else if dist >2000 and dist LE 2500 then flag_bin=5;
  
  if dist > 2500 then delete;

  if R2 LE 0.1 then flag_ld =1;
   else if R2 > 0.1 and R2 LE 0.2 then flag_ld=2;
   else if R2 > 0.2 and R2 LE 0.3 then flag_ld=3;
   else if R2 > 0.3 and R2 LE 0.4 then flag_ld=4;
   else if R2 > 0.4 and R2 LE 0.5 then flag_ld=5;
   else if R2 > 0.5 and R2 LE 0.6 then flag_ld=6;
   else if R2 > 0.6 and R2 LE 0.7 then flag_ld=7;
   else if R2 > 0.7 and R2 LE 0.8 then flag_ld=8;
   else if R2 > 0.8 and R2 LE 0.9 then flag_ld=9;
   else if R2 > 0.9 and R2 LE   1 then flag_ld=10;

  run;

data ld3r_bin;
  set ld3r;
  if dist LE 250 then flag_bin=1;
   else if dist > 250 and dist LE 500 then flag_bin=2;
   else if dist > 500 and dist LE 1000 then flag_bin=3;
   else if dist > 1000 and dist LE 2000 then flag_bin=4;
   else if dist >2000 and dist LE 2500 then flag_bin=5;
  
  if dist > 2500 then delete;

  if R2 LE 0.1 then flag_ld =1;
   else if R2 > 0.1 and R2 LE 0.2 then flag_ld=2;
   else if R2 > 0.2 and R2 LE 0.3 then flag_ld=3;
   else if R2 > 0.3 and R2 LE 0.4 then flag_ld=4;
   else if R2 > 0.4 and R2 LE 0.5 then flag_ld=5;
   else if R2 > 0.5 and R2 LE 0.6 then flag_ld=6;
   else if R2 > 0.6 and R2 LE 0.7 then flag_ld=7;
   else if R2 > 0.7 and R2 LE 0.8 then flag_ld=8;
   else if R2 > 0.8 and R2 LE 0.9 then flag_ld=9;
   else if R2 > 0.9 and R2 LE   1 then flag_ld=10;

  run;

data ldx_bin;
  set ldx;
  if dist LE 250 then flag_bin=1;
  else if dist > 250 and dist LE 500 then flag_bin=2;
  else if dist > 500 and dist LE 1000 then flag_bin=3;
  else if dist > 1000 and dist LE 2000 then flag_bin=4;
  else if dist >2000 and dist LE 2500 then flag_bin=5;
  if dist > 2500 then delete;

  if R2 LE 0.1 then flag_ld =1;
  else if R2 > 0.1 and R2 LE 0.2 then flag_ld=2;
  else if R2 > 0.2 and R2 LE 0.3 then flag_ld=3;
  else if R2 > 0.3 and R2 LE 0.4 then flag_ld=4;
  else if R2 > 0.4 and R2 LE 0.5 then flag_ld=5;
  else if R2 > 0.5 and R2 LE 0.6 then flag_ld=6;
  else if R2 > 0.6 and R2 LE 0.7 then flag_ld=7;
  else if R2 > 0.7 and R2 LE 0.8 then flag_ld=8;
  else if R2 > 0.8 and R2 LE 0.9 then flag_ld=9;
  else if R2 > 0.9 and R2 LE   1 then flag_ld=10;

  run;

*Save the datasets;
data dsim.chr4_ld_bins


proc freq data=ld4_bin ;
  tables flag_bin*flag_ld / out=ld4_out;
  run;

proc freq data=ld2l_bin ;
  tables flag_bin*flag_ld / out=ld2l_out;
  run;

proc freq data=ld2r_bin ;
  tables flag_bin*flag_ld / out=ld2r_out;
  run;

proc freq data=ld3r_bin ;
  tables flag_bin*flag_ld / out=ld3r_out;
  run;

proc freq data=ld3l_bin ;
  tables flag_bin*flag_ld / out=ld3l_out;
  run;

proc freq data=ldx_bin ;
  tables flag_bin*flag_ld / out=ldx_out;
  run;

* Export the proc freq results;

proc export data=ldx_out
  outfile="!MCLAB/ethanol/Sim_Pop_Gen/output/LD/chrX_LD_proc_freq_by_dist.txt"
  DBMS=TAB REPLACE; run;

proc export data=ld4_out
  outfile="!MCLAB/ethanol/Sim_Pop_Gen/output/LD/chr4_LD_proc_freq_by_dist.txt"
  DBMS=TAB REPLACE; run;

proc export data=ld2l_out
  outfile="!MCLAB/ethanol/Sim_Pop_Gen/output/LD/chr2L_LD_proc_freq_by_dist.txt"
  DBMS=TAB REPLACE; run;

proc export data=ld2r_out
  outfile="!MCLAB/ethanol/Sim_Pop_Gen/output/LD/chr2R_LD_proc_freq_by_dist.txt"
  DBMS=TAB REPLACE; run;

proc export data=ld3r_out
  outfile="!MCLAB/ethanol/Sim_Pop_Gen/output/LD/chr3R_LD_proc_freq_by_dist.txt"
  DBMS=TAB REPLACE; run;

proc export data=ld3l_out
  outfile="!MCLAB/ethanol/Sim_Pop_Gen/output/LD/chr3L_LD_proc_freq_by_dist.txt"
  DBMS=TAB REPLACE; run;





*Now do proc univariate on the data to look at distribution of distance vs ld;
proc sort data=ld4_bin;
  by dist R2;
  run;

proc univariate data=ld4_bin noprint;
  *var R2;
  class flag_bin;
  by dist R2;
   output out=ld4_out;
  run;
