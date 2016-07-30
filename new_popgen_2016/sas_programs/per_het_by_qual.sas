/*Compare the quality score with the percent het at positions */

libname sim "!MCLAB/ethanol/Sim_Pop_Gen/sas_data";

* First on chromosome 4;

data vcf4;
 set sim.vcf_4;
 drop ID FILTER INFO FORMAT;
 run;


*Now need to make homo=0, het=1, missing=.;
*Then get the average across the row for the  percent het;
data vcf4_c;
  set vcf4_b;
  array rows _character_;
    do over rows;
    if rows = './.' then rows = '.';
    else if rows='1/1' or rows='0/0' then rows = 0;
    else if rows = '0/1' or rows = '1/0' then rows = 1;
    per_het=mean(of Sz100 -- Sz99);
    end;
  run;

*Now want to compare % het and qual;
proc plot data=vcf4_c;
  plot qual*per_het;
  title 'Percent Het vs Quality';
  title1 'By Position';
    run;

data het_qual;
  set vcf4_c;
  keep chrom pos qual per_het;
  run;
proc export data=het_qual
    outfile=‘!MCLAB/ethanol/Sim_Pop_Gen/output/checks/het_vs_qual_check/chr4_het_qual.txt'
    DBMS=TAB REPLACE;
    run;


* Make macros;
%macro drop(chrom);
data vcf&chrom;
 set sim.vcf_&chrom;
 drop ID FILTER INFO FORMAT;
 run;
%mend;

%drop(4);
%drop(2L);
%drop(2R);
%drop(3L);
%drop(3R);
%drop(X);


*Now need to make homo=0, het=1, missing=.;
*Then get the average across the row for the  percent het;
%macro avg(chrom);
data vcf&chrom._b;
  set vcf&chrom;
  array rows _character_;
    do over rows;
    if rows = './.' then rows = '.';
    else if rows='1/1' or rows='0/0' then rows = 0;
    else if rows = '0/1' or rows = '1/0' then rows = 1;
    per_het=mean(of Sz100 -- Sz99);
    end;
  run;
%mend;

%avg(4);
%avg(2L);
%avg(2R);
%avg(3L);
%avg(3R);
%avg(X);

*Now want to compare % het and qual;
*Going to export the het and qual columns and plot in R;
%macro outR(chrom);
data het_qual_&chrom;
  set vcf&chrom._b;
  keep chrom pos qual per_het;
  run;
proc export data=het_qual_&chrom
    outfile='/home/fnew/dsim/vcf/chr&chrom._het_qual.txt'
    DBMS=TAB REPLACE;
    run;    
    %mend;
%outR(4);
%outR(2L);
%outR(2R);
%outR(3L);
%outR(3R);
%outR(X);
