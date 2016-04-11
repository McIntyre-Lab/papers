/* I need to parse the allele frequency output from vcftools so that I can plot the distribution*/

libname pop "!MCLAB/ethanol/Sim_Pop_Gen/sas_data";
libname dsim "/home/fnew/dsim";

*Import the allele freq file, be careful of the fifth column;

data dsim.allele_freq;
  set fnew.allele_freq;
  run;

*Set missing values to 0;
data allele_freqs2;
  set allele_freqs;
  array change _character_;
    do over change;
    if change=' ' then change=0;
    end;
  run;

*How many bi-, tri-, and tetra-allelic loci are there? ;

*Create flags;
data allele_freqs3;
  set allele_freqs2;
  if N_alleles = 1 then flag_mono=1; else flag_mono=0;
  if N_alleles = 2 then flag_bi=1; else flag_bi=0;
  if N_alleles = 3 then flag_tri=1; else flag_tri=0;
  if N_alleles = 4 then flag_tetra=1; else flag_tetra=0;
  run;



proc freq data=allele_freqs3;
  tables flag_:;
  run;
/*
bi: 94.66%      5284891
tri: 5.05%      282024
tetra: 0.29%    16085
*/




*Drop all info except allele_freqs;
data mono;
  set allele_freq_2;
  if flag_mono=1;
  run;

data di;
  set allele_freqs3;
  if flag_bi=1;
  run;

data tri;
  set allele_freqs3;
  if flag_tri=1;
  run;

data tetra;
  set allele_freqs3;
  if flag_tetra=1;
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

data tetra2;
  set tetra;
  keep allele_freq allele_freq_2 allele_freq_3 allele_freq_4;
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

data tetra3;
  set tetra2;
  af1=substr(allele_freq,3);
  af2=substr(allele_freq_2,3);
  af3=substr(allele_freq_3,3);
  af4=substr(allele_freq_4,3);
  keep af1 af2 af3 af4;
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

data minTETRA;
  set tetra3;
  array x {*} _char_;
  min=min(of x[*]);
  run;



*Export the allele frequencies;
proc export data = minDI
	outfile="!MCLAB/ethanol/Sim_Pop_Gen/output/dmel_biallele_frequencies.csv"
	dbms=CSV REPLACE;
	run;


proc export data = minTRI
	outfile="!MCLAB/ethanol/Sim_Pop_Gen/output/dmel_triallele_frequencies.csv"
	dbms=CSV REPLACE;
	run;

proc export data = minTETRA
	outfile="!MCLAB/ethanol/Sim_Pop_Gen/output/dmel_tetraallele_frequencies.csv"
	dbms=CSV REPLACE;
	run;





















