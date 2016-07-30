/*
The SelectHapStats.py program requires a specific format 
    Files separated by chromosome
    Column 1 is position, from smallest to largest
    The remaining columns are the nucleotide state for each individual
*/

/**** VCFs are not saved on MCLAB — they are on HPC: dsim_pop/Haplotype_Caller_Merged/vcf_by_chrom/*.vcf ****/

libname dsim "!MCLAB/ethanol/Sim_Pop_Gen/sas_data";

proc import datafile = '/home/fnew/dsim/vcf/no_missing/no_missing_biall_chr_nolab_nomel_noindel_2R.recode.vcf'
    out=vcf_2R
    dbms=TAB REPLACE;
    delimiter='09'x;
    datarow=2; getnames=yes;
    run;

proc import datafile = '/home/fnew/dsim/vcf/no_missing/no_missing_biall_chr_nolab_nomel_noindel_2L.recode.vcf'
    out=vcf_2L
    dbms=TAB REPLACE;
    delimiter='09'x;
    datarow=2; getnames=yes;
    run;

proc import datafile = '/home/fnew/dsim/vcf/no_missing/no_missing_biall_chr_nolab_nomel_noindel_3R.recode.vcf'
    out=vcf_3R
    dbms=TAB REPLACE;
    delimiter='09'x;
    datarow=2; getnames=yes;
    run;

proc import datafile = '/home/fnew/dsim/vcf/no_missing/no_missing_biall_chr_nolab_nomel_noindel_3L.recode.vcf'
    out=vcf_3L
    dbms=TAB REPLACE;
    delimiter='09'x;
    datarow=2; getnames=yes;
    run;

proc import datafile = '/home/fnew/dsim/vcf/no_missing/no_missing_biall_chr_nolab_nomel_noindel_4.recode.vcf'
    out=vcf_4
    dbms=TAB REPLACE;
    delimiter='09'x;
    datarow=2; getnames=yes;
    run;

proc import datafile = '/home/fnew/dsim/vcf/no_missing/no_missing_biall_chr_nolab_nomel_noindel_X.recode.vcf'
    out=vcf_X
    dbms=TAB REPLACE;
    delimiter='09'x;
    datarow=2; getnames=yes;
    run;
/* Need to substring the columns for genotype, only want the 0/1, 0/0, or 1/1*/
macro %trm(chrom);
data vcf_&chrom._nomiss;
  set vcf_&chrom;
  drop ID qual filter info format;
  array rows _character_;
    do over rows;
    rows = substr(rows,1,3);
    end;
    run;
%mend;

%trm(4);
%trm(2L);
%trm(2R);
%trm(3L);
%trm(3R);
%trm(X);


/* First on chrom4, if sample = 0/0 then sample = ref;
                    if sample = 0/1 then sample = 'N';
                    if sample = 1/1 then sample = alt;
*/

data vcf4_b;
  set dsim.vcf_4_nomiss;
  array rows _character_;
    do over rows;
    if rows = "./." then rows = "N";
    else if rows = "1/1" then rows = alt;
    else if rows = "0/0" then rows = ref;
    else if rows = "0/1" or rows = "1/0" then rows = "N";
    end;
  run;

*checks;
proc freq data=dsim.vcf_4_nomiss;
  tables Sz101;
  run; *0/1: 35;
proc freq data=vcf4_b;
  tables Sz101;
  run; *N: 35;

*make it a macro;
%macro reformat(chrom);
data vcf&chrom._b;
  set dsim.vcf_&chrom._nomiss;
  array rows _character_;
    do over rows;
    if rows = "./." then rows = "N";
    else if rows = "1/1" then rows = alt;
    else if rows = "0/0" then rows = ref;
    else if rows = "0/1" or rows = "1/0" then rows = "N";
    end;
  run;
%mend;

%reformat(2L);
%reformat(2R);
%reformat(3L);
%reformat(3R);
%reformat(X);

*checks;
proc freq data=dsim.vcf_2l_nomiss;
  tables Sz101;
  run; *0/1 = 7875;
proc freq data=vcf2l_b;
  tables Sz101;
  run; *N = 7875;


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
    outfile="!MCLAB/ethanol/Sim_Pop_Gen/output/vcf_for_h12/chr&chrom._for_h12.csv"
    DBMS=CSV REPLACE;
    run;
%mend;

%exp_vcf(2L);
%exp_vcf(2R);
%exp_vcf(3L);
%exp_vcf(3R);
%exp_vcf(4);
%exp_vcf(X);

/*Save the datasets*/
data dsim.vcf_4_nomiss_4H12;
  set vcf4_c;
  run;
data dsim.vcf_2L_nomiss_4H12;
  set vcf2L_c;
  run;
data dsim.vcf_2R_nomiss_4H12;
  set vcf2R_c;
  run;
data dsim.vcf_3R_nomiss_4H12;
  set vcf3R_c;
  run;
data dsim.vcf_3L_nomiss_4H12;
  set vcf3L_c;
  run;
data dsim.vcf_X_nomiss_4H12;
  set vcfX_c;
  run;

