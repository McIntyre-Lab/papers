/*
The SelectHapStats.py program requires a specific format 
    Files separated by chromosome
    Column 1 is position, from smallest to largest
    The remaining columns are the nucleotide state for each individual
*/

libname dsim "!MCLAB/ethanol/Sim_Pop_Gen/sas_data";


*There might be a header on the vcf still, check;
proc import datafile = '/home/fnew/dsim/checking_gatk_output/analysis_gvcf_NOmiss_nolabmelindel.recode.vcf'
    out=vcf_x
    dbms=TAB REPLACE;
    delimiter='09'x;
    datarow=2; getnames=yes;
    run;

/* Need to substring the columns for genotype, only want the 0/1, 0/0, or 1/1*/
%macro trm(chrom);
data vcf_&chrom._nomiss;
  set vcf_&chrom;
  drop ID qual filter info format;
  array rows _character_;
    do over rows;
    rows = substr(rows,1,3);
    end;
    run;
%mend;

%trm(x);


/* if sample = 0/0 then sample = ref;
   if sample = 0/1 then sample = 'N';
   if sample = 1/1 then sample = alt;
*/
%macro reformat(chrom);
data vcf&chrom._b;
  set vcf_&chrom._nomiss;
  array rows _character_;
    do over rows;
    if rows = "./." then rows = "N";
    else if rows = "1/1" then rows = alt;
    else if rows = "0/0" then rows = ref;
    else if rows = "0/1" or rows = "1/0" then rows = "N";
    end;
  run;
%mend;

%reformat(x);


/* I only need position and samples for each chromosome, then export */
%macro exp_vcf(chrom);
data vcf&chrom._c;
  set vcf&chrom._b;
  drop _chrom ID ref alt qual filter info format;
  run;
proc sort data=vcf&chrom._c;
  by pos;
  run;
proc export data=vcf&chrom._c
    outfile=“!MCLAB/ethanol/Sim_Pop_Gen/output/checking_gatk_output/gvcf/output/data/chr&chrom._for_h12.csv"
    DBMS=CSV REPLACE;
    run;
%mend;

%exp_vcf(X);

