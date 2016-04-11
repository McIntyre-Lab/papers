/* I need to parse the allele frequency output from vcftools so that I can plot the distribution*/

libname dsim "!MCLAB/ethanol/Sim_Pop_Gen/sas_data";
libname fnew "/home/fnew/dsim";

*Import the allele freq file, be careful of the fifth column;

data dsim.allele_freqs_exon;
  set fnew.allele_freq_exon;
  run;

*Set missing values to 0;
data dsim.allele_freqs_exon;
  set dsim.allele_freqs_exon;
  array change _character_;
    do over change;
    if change=' ' then change=0;
    end;
  run;

*How many bi-, tri-allelic loci are there? ;

*Create flags!;
data allele_freq_exon_2;
  set dsim.allele_freqs_exon;
  if N_alleles = 1 then flag_mono=1; else flag_mono=0;
  if N_alleles = 2 then flag_bi=1; else flag_bi=0;
  if N_alleles = 3 then flag_tri=1; else flag_tri=0;
  if N_alleles = 4 then flag_tetra=1; else flag_tetra=0;
  run;

proc freq data=allele_freq_exon_2;
  tables flag_:;
  run;

/* mono:    953491      96.20%
   bi:      37,579      3.79%
   tri:     94          0.01%
   tetra:   0           0.00%
*/

*Drop all info except allele_freqs;
data mono;
  set allele_freq_exon_2;
  if flag_mono=1;
  run;

data di;
  set allele_freq_exon_2;
  if flag_bi=1;
  run;

data tri;
  set allele_freq_exon_2;
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
	outfile="!MCLAB/ethanol/Sim_Pop_Gen/output/biallele_frequencies_exon.csv"
	dbms=CSV REPLACE;
	run;


proc export data = minTRI
	outfile="!MCLAB/ethanol/Sim_Pop_Gen/output/triallele_frequencies_exon.csv"
	dbms=CSV REPLACE;
	run;


