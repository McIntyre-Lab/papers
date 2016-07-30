/* For the Swedish R code 'getPopstat.R' I need a bed file with invariant site
 * positions. The only way I can think to get this is the following:
    1) Import a windowed VCF
    2) Make another file that counts 1 by 1 all the positions from vcf_start to vcf_end
        2b) Also, the next column should be col1+1, to give the position of each base
    3) Merge these files and output the sites that are not in the VCF
    4) This result should be like a bed file with all positions that are not in the VCF...
    */
libname chk "!MCLAB/ethanol/Sim_Pop_Gen/sas_data/checks";

%macro import(chr,window);
proc import datafile = ‘!MCLAB/ethanol/Sim_Pop_Gen/output/checks/getPopStat_test/dsim_10kb_chrom4_6.vcf'
    out=vcf
    dbms=TAB REPLACE;
    delimiter='09'x;
    datarow=2; getnames=yes;
    run;
    *208 obs;

*In case I need it later...;
data chk.vcf_&chr._window&window._10kb;
  set vcf;
  run;

/*I need a file that is three columns:
    1) chrom
    2) start
    3) end
  These should be the positions from start to finish in the VCF, with no gaps
  */

proc import datafile ="!MCLAB/ethanol/Sim_Pop_Gen/output/checks/getPopStat_test/chr4_w6_all_sites.csv"
    out = all_sites
    dbms = CSV REPLACE;
    datarow=2; getnames=yes;
    run;

*Merge these two files together and output what is not in1 and in2;
proc sort data=vcf; by pos; run;
proc sort data=all_sites; by pos; run;

data variants invariants;
  merge vcf (in=in1) all_sites (in=in2);
  by pos;
  if in1 and in2 then output variants;
  else if in2 then output invariants;
  run;
* variants has 208 obs;
* invariants has 9,748 obs;


*Export the invariants file;
proc export data=invariants
    outfile="!MCLAB/ethanol/Sim_Pop_Gen/output/checks/getPopStat_test/chr4_w6_invariants_for_R.txt"
    dbms=TAB REPLACE;
    run;


%mend;

%import(3R,204);
%import(2R,31);
%import(4,6);
