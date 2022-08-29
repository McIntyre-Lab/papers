

/* create GA NMR design files with only PD1074 samples
    one with all PD1074 and second with only PD1074 in more than 1 batch
*/

%macro subbing (datain, dump) ;

proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/data_from_group/nmr_sweet16/Design_file_NMR_&datain._14Aug2020.txt"
out = dsgn_&datain. 
replace ;
guessingrows = MAX ;
run ;

data dsgn2_&datain. ;
set dsgn_&datain. ;
sample = compress(scan(wormgrowth_sample_name, 1, '_'));
newSampleID = compress(tranwrd(sampleID, "_&dump.", ''));
rename sampleID = oldSampleID ;
run;

proc sort data = dsgn2_&datain. ;
by sample ;
run;

data dsgn_NMR_&datain. ;
retain newSampleID ;
set dsgn2_&datain.;
rename newSampleID = sampleID ;
run ;

proc export data =  dsgn_NMR_&datain.
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/data_from_group/nmr_sweet16/dsgn_NMR_&datain._shortNames.txt"
dbms = tab 
replace ;
run ;



%mend ;

%subbing (CDCl3, NMR_CDCL3) ;

%subbing (D2O, NMR_D2O) ;

