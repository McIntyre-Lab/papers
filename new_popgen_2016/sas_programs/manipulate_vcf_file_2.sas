/* Manipulate vcf file to create binary indicator for het/hom/missing data */

libname sim "/home/fnew/mclab/ethanol/Sim_Pop_Gen/sas_data";

%macro import(chr);
proc import datafile = '/home/fnew/dsim/vcf/filter_10permiss_nolab_nomel_4.recode.vcf'
    out=vcf_4
    dbms=TAB REPLACE;
    delimiter='09'x;
    datarow=2; getnames=yes;
    run;
    %mend;

%import(2L);
%import(2R);
%import(3L);
%import(3R);
%import(4);
%import(X);

proc import datafile = '/home/fnew/dsim/vcf/filter_10permiss_nolab_nomel_2R.recode.vcf'
    out=vcf_2R
    dbms=TAB REPLACE;
    delimiter='09'x;
    datarow=2; getnames=yes;
    run;
proc import datafile = '/home/fnew/dsim/vcf/filter_10permiss_nolab_nomel_2L.recode.vcf'
    out=vcf_2L
    dbms=TAB REPLACE;
    delimiter='09'x;
    datarow=2; getnames=yes;
    run;
proc import datafile = '/home/fnew/dsim/vcf/filter_10permiss_nolab_nomel_3R.recode.vcf'
    out=vcf_3R
    dbms=TAB REPLACE;
    delimiter='09'x;
    datarow=2; getnames=yes;
    run;
proc import datafile = '/home/fnew/dsim/vcf/filter_10permiss_nolab_nomel_3L.recode.vcf'
    out=vcf_3L
    dbms=TAB REPLACE;
    delimiter='09'x;
    datarow=2; getnames=yes;
    run;
proc import datafile = '/home/fnew/dsim/vcf/filter_10permiss_nolab_nomel_X.recode.vcf'
    out=vcf_X
    dbms=TAB REPLACE;
    delimiter='09'x;
    datarow=2; getnames=yes;
    run;

%include '/home/fnew/mclab/ethanol/Sim_Pop_Gen/sas_programs/sample_list.txt';
%put &mylist;
%include '/home/fnew/mclab/macros/iterlist.sas';

*I only want the genotype field;
%macro trimm(line, chr);
    data &line;
    set vcf_&chr ;
    &line = substr(&line,1,3);
    run;
    %mend;

%combine()

proc freq data=sim.vcf_X_b; run;
proc freq data=sim.vcf_4_b;  run;
proc freq data=sim.vcf_3L_b; run;
proc freq data=sim.vcf_3R_b; run;
proc freq data=sim.vcf_2L_b; run;
proc freq data=sim.vcf_2R_b; run;

data vcf_4_b;
  set sim.vcf_4_b;
  cchex=put(_CHROM, hex2.);
  drop _CHROM;
  rename cchex=_CHROM;
  retain _CHROM;
  run;
data vcf_x_b;
  set sim.vcf_x_b;
  length _CHROM $ 2 ;
  run;

*Combine the chromosome counts;
data sim.all_chr;
  set vcf_4_b sim.vcf_3l_b sim.vcf_3r_b vcf_x_b sim.vcf_2r_b sim.vcf_2l_b;
  run;


*Summarize counts...;
proc freq data=all_chr; run;


data sim.vcf_4_b;
  set sim.vcf_4;
  drop ID REF ALT FILTER INFO FORMAT;
  rename _CHROM=CHROM;
  run;

  proc export data=sim.vcf_4_b 
    outfile="/home/fnew/dsim/vcf/vcf_4_for_plot.txt"
    DBMS=TAB REPLACE;
    run;


