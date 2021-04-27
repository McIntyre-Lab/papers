
libname pacbio "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata";
libname seq "!MCLAB/maize_ozone_FINAL/2018/RNAseq_Novogene/sasdata" ;

filename mymacros "!MCLAB/maize_ozone_FINAL/2018/PacBio/sas_programs/macros";
options SASAUTOS=(sasautos mymacros);



/*
import blast results of unmapped
    note that have unmapped for each sample - reference combo == 10 samples * 3 refs = 30
*/

data samples ;
set seq.df_ng_1st_lane_sampleid ;
length trt $3.;
where pacbio = 1 ;
if ozone_trt = "Ele" then trt = "oz" ;
else  trt = "amb" ;
ID = compress(tubeID||'_'||inbred_line||'_'||trt) ;
newID = tranwrd(ID, "Hp30_", "Hp301_") ;
newID2 = tranwrd(newID, "NC33_", "NC338_") ;
drop sampleID ;
rename newID2 = sampleID ;
keep newID2 ;
run;

proc sort data = samples nodups ;
by _all_ ;
run;

%macro adding (ref) ;
data sample_&ref ;
length ref $8.;
set samples ;
ref = "&ref";
run ;
%mend ;
%adding (b73) ;
%adding (mo17_yan) ;
%adding (mo17_cau) ;

data sample_list ;
retain sampleID ref ;
set sample_: ;
run;


ods pdf file = "!MCLAB/maize_ozone_FINAL/2018/PacBio/cnts_unmapped_aln_length_50per_qlen.pdf" ;

%macro importing (sampleid, ref) ;

proc import datafile = "!MCLAB/maize_ozone_FINAL/2018/PacBio/blast_unmapped_PBreads/blast_unmapped_&sampleID._unCollapsed_unmapped_&ref._ref.fa_vs_nt.tsv"
out = &ref._&sampleid. 
dbms = tab replace ;
run;

/* flag if (1) aln length ge 50% of PB length and (2) if query length is less than 500bp  */
data ref2_&ref._&sampleid.  ;
length stitle $200.;
length sampleID $26.;
set &ref._&sampleid. ;
if length ge (.5*qlen) then flag_50per = 1;
else flag_50per = 0;
if qlen ge 500 then flag_qlen_500 = 1;
else flag_qlen_500 = 0 ;
sampleID = "&sampleid._ref_&ref." ;
run;

title "counts where aln length ge 50% of PB length for &sampleid. - &ref. ref" ;
proc freq data = ref2_&ref._&sampleid.  ;
tables flag_50per / out = cnts_&ref._&sampleid.;
run;

title "counts where qlen ge 500  for &sampleid. - &ref. ref" ;
proc freq data = ref2_&ref._&sampleid.  ;
tables flag_qlen_500 / out = cnts_&ref._&sampleid.;
run;

title "counts where qlen ge 500 and aln length ge 50% of PB length for &sampleid. - &ref. ref" ;
proc freq data = ref2_&ref._&sampleid.  ;
tables flag_50per * flag_qlen_500 / out = cnts_&ref._&sampleid.;
run;

title ;

%mend ;

%iterdataset(dataset=sample_list, function=%nrstr(%importing(&sampleid, &ref);));
ods pdf close ;

data pacbio.blast_sqanti_unmapped_reads ;
retain sampleID qseqid stitle;
set ref2_: ;
run ;

proc export data = pacbio.blast_sqanti_unmapped_reads
outfile = "!MCLAB/maize_ozone_FINAL/2018/PacBio/blast_sqanti_unmapped_reads.csv"
dbms = csv replace ;
run;

proc sort data = pacbio.blast_sqanti_unmapped_reads ;
by sampleID ;
run;

proc freq data = pacbio.blast_sqanti_unmapped_reads ;
tables flag_50per * flag_qlen_500 ;
by sampleID ;
run;








