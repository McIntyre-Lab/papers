
libname pacbio "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata";
libname seq "!MCLAB/maize_ozone_FINAL/2018/RNAseq_Novogene/sasdata" ;

filename mymacros "!MCLAB/maize_ozone_FINAL/2018/PacBio/sas_programs/macros";
options SASAUTOS=(sasautos mymacros);



/*
import blast results of non-b73 and non-mo17 PB reads that do NOT map to b73 and do NOT map to mo17

flags created:
    if length ge (.5*qlen) then flag_alnLen_50per = 1;
    if qlen ge 500 then flag_qlen_500 = 1;

    if length ge (.8*qlen) then flag_alnLen_80per = 1;
    if qlen ge 1500 then flag_qlen_1500 = 1;


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

data sample_list ;
set samples ;
if find(sampleID, "Mo17") ge 1 then delete ;
if find(sampleID, "B73") ge 1 then delete ;
run;



ods pdf file = "!MCLAB/maize_ozone_FINAL/2018/PacBio/evidence_genotype_specific_seqs/freqs_evidence_genotype_spec_seqs.pdf" ;

%macro importing (sampleid) ;

proc import datafile = "!MCLAB//maize_ozone_FINAL/2018/PacBio/evidence_genotype_specific_seqs/blast_seqs_no_aln_&sampleid..fa_vs_nt.tsv"
out = s_&sampleid. 
dbms = tab replace ;
guessingrows = MAX ;
run;

title "query count for &sampleid." ;
proc freq data = s_&sampleID ;
tables qseqid / out = qseq_&sampleID ;
run;

/* keep top hit based on bit score */
proc sort data = s_&sampleid. ;
by qseqid descending bitscore ;
run;

data s1_&sampleid.  ;
set s_&sampleid. ;
by qseqid ;
if first.qseqid ;
run ;


/* flag if (1) aln length ge 50% of PB length and (2) if query length is ge to 500bp  */
data s2_&sampleid.  ;
length stitle $200.;
length sampleID $26.;
set s1_&sampleid. ;
if length ge (.5*qlen) then flag_alnLen_50per = 1;
else flag_alnLen_50per = 0;
if qlen ge 500 then flag_qlen_500 = 1;
else flag_qlen_500 = 0 ;

if length ge (.8*qlen) then flag_alnLen_80per = 1;
else flag_alnLen_80per = 0;
if qlen ge 2000 then flag_qlen_2000 = 1;
else flag_qlen_2000 = 0 ;

sampleID = "&sampleid" ;
run;

proc sort data = s2_&sampleid.  ; ;
by descending flag_alnLen_50per flag_qlen_500;
run;

title "counts where qlen ge 500 and aln length ge 50% of PB length for &sampleid." ;
proc freq data = s2_&sampleid.  ;
tables flag_alnLen_50per * flag_qlen_500 / out = cnts_&sampleid.;
run;

title "counts where qlen ge 2000 and aln length ge 80% of PB length for &sampleid." ;
proc freq data = s2_&sampleid.  ;
tables flag_alnLen_80per * flag_qlen_2000 / out = cnts2_&sampleid.;
run;

title ;

%mend ;

%iterdataset(dataset=sample_list, function=%nrstr(%importing(&sampleid);));
ods pdf close ;


data blast_unmapped_reads ;
retain sampleID qseqid stitle;
set s2_: ;
run ;


data  blast_unmapped_reads_filter1 ;
set  blast_unmapped_reads ;
where flag_alnLen_50per = 1 and flag_qlen_500 = 1 ;
run ;  

data  blast_unmapped_reads_filter2 ;
set  blast_unmapped_reads ;
where flag_alnLen_80per = 1 and flag_qlen_2000 = 1 ;
run ;




