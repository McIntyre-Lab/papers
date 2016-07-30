/*Manipulate vcfs so that I can get percent heterozygosity by individual and
 * position */

 /*Going to do this locally for now */

 libname loc "/home/fnew/dsim/checking_gatk_output/sas_data";

 *Import chromosome X vcf using import wizard;

 data loc.vcf_chrx;
   set chrx;
   run;

data vcfx;
  set loc.vcf_chrx;
  drop ID FILTER INFO FORMAT;
  run;

*Manipulate file, only need genotype field;
data vcfx_fix;
  set vcfx;
  array rows _character_;
    do over rows;
    rows=substr(rows,1,3);
    end;
    run;

*Now need to make homo=0, het=1, missing=.;
*Then get the average across the row for the percent het by pos;

data vcfx_fix2;
  set vcfx_fix;
  array rows _character_;
    do over rows;
    if rows= './.' then rows='.';
    else if rows = '1/1' or rows = '0/0' then rows = 0;
    else if rows = '1/0' or rows = '0/1' then rows = 1;
    per_het=mean(of Sz100 -- Sz99);
    end;
   run;


*Now want the average het by individual;
proc means data=vcfx_fix2 nway noprint;
  var S: ;
  output out=individ_per_het (drop=_:) mean=average;
  run;

proc freq data=vcfx_fix2; tables S: ;
run;
