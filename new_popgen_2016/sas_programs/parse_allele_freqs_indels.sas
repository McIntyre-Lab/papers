/* I need to parse the allele frequency output from vcftools so that I can plot the distribution*/

libname dsim "!MCLAB/ethanol/Sim_Pop_Gen/sas_data";
libname fnew "/home/fnew/dsim";

*Import the allele freq file, be careful of the fifth column;

data dsim.indel_allele_freq;
  set fnew.indel_allele_freq;
  run;

*Set missing values to 0;
data dsim.indel_allele_freq;
  set dsim.indel_allele_freq;
  array change _character_;
    do over change;
    if change=' ' then change=0;
    end;
  run;

*How many bi-, tri-, and tetra-allelic loci are there? ;

*Create flags!;
data allele_freq_2;
  set dsim.indel_allele_freq;
   if allele_freq_2=0 and allele_freq_3=0 and allele_freq_4=0 then flag_mono=1; else flag_mono=0;
   if allele_freq_2 ne 0 and allele_freq_3=0 and allele_freq_4=0 then flag_di=1; else flag_di=0;
   if allele_freq_3 ne 0 and allele_freq_4=0 then flag_tri=1; else flag_tri=0;
   if allele_freq_4 ne 0 and allele_freq_5=0 and allele_freq_6=0 and allele_freq_7=0 and allele_freq_8=0 then flag_tetra=1; else flag_tetra=0;
   if allele_freq_5 ne 0 and allele_freq_6=0 and allele_freq_7=0 and allele_freq_8=0 then flag_five=1; else flag_five=0;
   if allele_freq_6 ne 0 and allele_freq_7=0 and allele_freq_8=0 then flag_six=1; else flag_six=0;
   if allele_freq_7 ne 0 and allele_freq_8=0 then flag_seven=1; else flag_seven=0;
   if allele_freq_8 ne 0 then flag_eight=1; else flag_eight=0;
  run;

proc freq data=allele_freq_2;
  tables flag_:;
  run;

/* mono:    0
   di:      182674      84.57%
   tri:     25845       11.97%
   tetra:   5479        2.54%
   five:    1611        0.75%
   six:     336
   seven:   50
*/

*Drop all info except allele_freqs;
data mono;
  set allele_freq_2;
  if flag_mono=1;
  run;

data di;
  set allele_freq_2;
  if flag_di=1;
  run;

data tri;
  set allele_freq_2;
  if flag_tri=1;
  run;

data tetra;
  set allele_freq_2;
  if flag_tetra=1;
  run;

data five;
  set allele_freq_2;
  if flag_five=1;
  run;

data six;
  set allele_freq_2;
  if flag_six=1;
  run;

data seven;
  set allele_freq_2;
  if flag_seven=1;
  run;

*Parse allele freqs into individal files;
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

data five2;
  set five;
  keep allele_freq allele_freq_2 allele_freq_3 allele_freq_4 allele_freq_5;
  run;

data six2;
  set six;
  keep allele_freq allele_freq_2 allele_freq_3 allele_freq_4 allele_Freq_5 allele_freq_6;
  run;

data seven2;
  set seven;
  keep allele_freq allele_freq_2 allele_freq_3 allele_freq_4 allele_freq_5 allele_freq_6 allele_freq_7;
  run;



*Remove the allele from the columns, only want the allele freq;
data di3;
  set di2;
  num=index(allele_freq,':');
    af1=substr(allele_freq,num+1);
  num2=index(allele_freq_2, ':');
    af2=substr(allele_freq_2,num2+1);
  keep af1 af2;
  run;

data tri3;
  set tri2;
  num=index(allele_freq,':');
    af1=substr(allele_freq,num+1);
  num2=index(allele_freq_2, ':');
    af2=substr(allele_freq_2,num2+1);
  num3=index(allele_freq_3, ':');
    af3=substr(allele_freq_3,num3+1);
  keep af1 af2 af3;
  run;

data tetra3;
  set tetra2;
  num=index(allele_freq,':');  
    af1=substr(allele_freq,num+1);
  num2=index(allele_freq_2,':');
    af2=substr(allele_freq_2,num2+1);
  num3=index(allele_freq_3,':');
    af3=substr(allele_freq_3,num3+1);
  num4=index(allele_freq_4,':');
    af4=substr(allele_freq_4,num4+1);
  keep af1 af2 af3 af4;
  run;

data five3;
  set five2;
  num=index(allele_freq,':');  
    af1=substr(allele_freq,num+1);
  num2=index(allele_freq_2,':');
    af2=substr(allele_freq_2,num2+1);
  num3=index(allele_freq_3,':');
    af3=substr(allele_freq_3,num3+1);
  num4=index(allele_freq_4,':');
    af4=substr(allele_freq_4,num4+1);
  num5=index(allele_freq_5, ':');
    af5=substr(allele_freq_5,num5+1);
  keep af1 af2 af3 af4 af5;
  run;


data six3;
  set six2;
  num=index(allele_freq,':');  
    af1=substr(allele_freq,num+1);
  num2=index(allele_freq_2,':');
    af2=substr(allele_freq_2,num2+1);
  num3=index(allele_freq_3,':');
    af3=substr(allele_freq_3,num3+1);
  num4=index(allele_freq_4,':');
    af4=substr(allele_freq_4,num4+1);
  num5=index(allele_freq_5, ':');
    af5=substr(allele_freq_5,num5+1);
  num6=index(allele_freq_6, ':');
    af6=substr(allele_freq_6,num6+1);
  keep af1 af2 af3 af4 af5 af6;
  run;

data seven3;
  set seven2;
  num=index(allele_freq,':');  
    af1=substr(allele_freq,num+1);
  num2=index(allele_freq_2,':');
    af2=substr(allele_freq_2,num2+1);
  num3=index(allele_freq_3,':');
    af3=substr(allele_freq_3,num3+1);
  num4=index(allele_freq_4,':');
    af4=substr(allele_freq_4,num4+1);
  num5=index(allele_freq_5, ':');
    af5=substr(allele_freq_5,num5+1);
  num6=index(allele_freq_6, ':');
    af6=substr(allele_freq_6,num6+1);
  num7=index(allele_freq_7, ':');
    af7=substr(allele_freq_7,num7+1);
  keep af1 af2 af3 af4 af5 af6 af7;
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

data minFIVE;
  set five3;
  array x {*} _char_;
  min=min(of x[*]);
  run;


data minSIX;
  set six3;
  array x {*} _char_;
  min=min(of x[*]);
  run;

data minSEVEN;
  set seven3;
  array x {*} _char_;
  min=min(of x[*]);
  run;

*Export the allele frequencies;
proc export data = minDI
	outfile="!MCLAB/ethanol/Sim_Pop_Gen/output/indel_biallele_frequencies.csv"
	dbms=CSV REPLACE;
	run;


proc export data = minTRI
	outfile="!MCLAB/ethanol/Sim_Pop_Gen/output/indel_triallele_frequencies.csv"
	dbms=CSV REPLACE;
	run;

proc export data = minTETRA
	outfile="!MCLAB/ethanol/Sim_Pop_Gen/output/indel_tetraallele_frequencies.csv"
	dbms=CSV REPLACE;
	run;
proc export data = minFIVE
	outfile="!MCLAB/ethanol/Sim_Pop_Gen/output/indel_5allele_frequencies.csv"
	dbms=CSV REPLACE;
	run;


proc export data = minSIX
	outfile="!MCLAB/ethanol/Sim_Pop_Gen/output/indel_6allele_frequencies.csv"
	dbms=CSV REPLACE;
	run;

proc export data = minSEVEN
	outfile="!MCLAB/ethanol/Sim_Pop_Gen/output/indel_7allele_frequencies.csv"
	dbms=CSV REPLACE;
	run;





















