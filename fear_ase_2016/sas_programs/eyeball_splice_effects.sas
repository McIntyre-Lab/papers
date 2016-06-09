
data check2;
set sig_anno;
*where symbol_cat = "ps" and line="w47" and mating_status="V";
where fusion_id="F50887_SI" /*(fusion_id="F50897_SI"or)*/ ;
keep fusion_id line sum_line sum_tester start end sum_both sum_total cis_i trans_i mean_apn;
run;



proc print data=check2;
run;

proc gplot data=check2;
where line="w47" and mating_status="V";
plot cis_i*trans_i=line;
run;

*gapdh looks real!;
*prat2 looks real; *check imp?;
*sif looks real;
*maybe yp3;
*ps w47 virgin! looks like large trans effect in one of the two fusions but potential complication in mapping;

proc sort data=check2;
by fusion_

proc sort data=check2;
by line mating_status fusion_id;
run;


proc univariate data=fear.Cis_data_estiamtes normal plot;
var cis_i trans_i;
run;

proc gplot data=fear.Cis_data_estiamtes;
*where cis_i <10 and cis_i >-10 and trans_i <10 and trans_i >-10;
plot cis_i*trans_i;
run;


data ps;
set fear.Fb551_si_fusions_unique_flagged;
where symbol_cat ="ps";
*keep fusion_id start end;
run;

proc sort data=ps;
by fusion_id;
proc sort data=fear.clean_ase_stack;
by fusion_id;

data all_ps;
merge fear.clean_ase_stack ps(in=in1);
by fusion_id;
if in1;
run;

proc sort data=all_ps;
by start;
run;

proc print data=all_ps;
var 

proc export data=all_ps outfile="c:\a1stuff\fear\all_ps.csv" dbms=csv;
run;

