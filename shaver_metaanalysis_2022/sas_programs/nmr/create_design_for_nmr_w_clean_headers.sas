

/* create GA NMR design files with only PD1074 samples
    one with all PD1074 and second with only PD1074 in more than 1 batch
*/

%macro subbing (datain, dump) ;

proc import datafile = "/home/ammorse/sweet16/nmr/Design_file_NMR_&datain._14Aug2020.txt"
out = dsgn_&datain. replace ;
guessingrows = MAX ;
run ;

data dsgn2_&datain. ;
set dsgn_&datain. ;
where genotype = "PD1074" ;
sample = compress(scan(wormgrowth_sample_name, 1, '_'));
newSampleID = compress(tranwrd(sampleID, "_&dump.", ''));
rename sampleID = oldSampleID ;
run;

proc sort data = dsgn2_&datain. ;
by sample ;
run;

data single_&datain. dup_&datain. ;
set dsgn2_&datain. ;
by sample ;
if first.sample and last.sample then output single_&datain. ;
else output dup_&datain.;
run ;

data dsgn_NMR_&datain._mult_PD1074 ;
retain newSampleID ;
set dup_&datain.;
rename newSampleID = sampleID ;
run ;

data dsgn_NMR_&datain._all_PD1074 ;
retain newSampleID ;
set dsgn2_&datain. ;
rename newSampleID = sampleID ;
run ;

proc export data = dsgn_NMR_&datain._mult_PD1074
outfile = "/home/ammorse/sweet16/nmr/dsgn_NMR_&datain._mult_PD1074.txt"
dbms = tab 
replace ;
run ;

proc export data = dsgn_NMR_&datain._all_PD1074
outfile = "/home/ammorse/sweet16/nmr/dsgn_NMR_&datain._all_PD1074.txt"
dbms = tab 
replace ;
run ;

%mend ;

%subbing (CDCl3, NMR_CDCL3) ;

%subbing (D2O, NMR_D2O) ;

