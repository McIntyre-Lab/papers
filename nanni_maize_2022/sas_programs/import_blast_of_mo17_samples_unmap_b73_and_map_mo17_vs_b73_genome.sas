
libname pacbio "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata";
libname seq "!MCLAB/maize_ozone_FINAL/2018/RNAseq_Novogene/sasdata" ;

filename mymacros "!MCLAB/maize_ozone_FINAL/2018/PacBio/sas_programs/macros";
options SASAUTOS=(sasautos mymacros);



/*
import blast results of mo17 reads that do NOT map to b73 but DO map to mo17 vs b73 genome

flags created:
    if length ge (.5*qlen) then flag_alnLen_50per = 1;
    if qlen ge 500 then flag_qlen_500 = 1;

    if length ge (.8*qlen) then flag_alnLen_80per = 1;
    if qlen ge 1500 then flag_qlen_1500 = 1;


*/


data samples ;
set pacbio.design_pacbio ;
newGeno1 = tranwrd(genotype, "hp301", "Hp301") ;
newGeno2 = tranwrd(newGeno1, "b73", "B73") ;
newGeno3 = tranwrd(newGeno2, "mo17", "Mo17") ;
newGeno4 = tranwrd(newGeno3, "c123", "C123") ;
Genotype2 = tranwrd(newGeno4, "nc338", "NC338") ;
drop newGeno: genotype sampleID sampleID2 otherID;
rename Genotype2 = genotype ;
sampleID3 = compress(ID||'_'||genotype2||'_'||trt) ;
rename sampleID3 = sampleID2 ;
sampleID4 = tranwrd(sampleID3, "21-2", "21_2") ;
rename sampleID4 = sampleID ;
run;


proc sort data = samples nodups ;
by _all_ ;
run;

data sample_list ;
retain sampleID sampleID2 ;
set samples ;
if find(sampleID, "Mo17") ge 1 ;
run;

proc sort data = sample_list nodups ;
by _all_ ;
run;


ods pdf file = "!MCLAB/maize_ozone_FINAL/2018/PacBio/evidence_genotype_specific_seqs/freqs_evidence_genotype_spec_seqs_mo17_samples_vs_b73_genome.pdf" ;

%macro importing (sampleid, sampleid2) ;

proc import datafile = "!MCLAB/maize_ozone_FINAL/2018/PacBio/evidence_genotype_specific_seqs/blast_seqs_yes_aln_&sampleid2..fa_vs_b73_genome.tsv"
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

data cnts3_&sampleID. ;
set cnts2_&sampleID ;
sampleID = "&sampleID";
run;

data s3_&sampleID ;
set s2_&sampleID ;
sampleID = "&sampleID";
run;


title ;

%mend ;

%importing(21_2_Mo17_oz, 21-2_Mo17_oz );

%iterdataset(dataset=sample_list, function=%nrstr(%importing(&sampleid, &sampleid2);));
ods pdf close ;


data blast_mo17_reads ;
retain sampleID qseqid stitle;
set s3_: ;
run ;

proc freq data = blast_mo17_reads ;
tables sampleID;
run;  /* 19 = 2
        21-1 = 16
        21 = 2  */

proc freq data = blast_mo17_reads ;
tables flag_alnLen_80per * flag_qlen_2000 ;
run;
       /* 
flag_alnLen_80per
          flag_qlen_2000

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |     14 |      3 |     17
         |  70.00 |  15.00 |  85.00
         |  82.35 |  17.65 |
         |  82.35 | 100.00 |
---------+--------+--------+
       1 |      3 |      0 |      3
         |  15.00 |   0.00 |  15.00
         | 100.00 |   0.00 |
         |  17.65 |   0.00 |
---------+--------+--------+
Total          17        3       20
            85.00    15.00   100.00

most (14 out of 20) of the transcripts are less than 2kb AND less than 80% of alnLen
None are ge 2kb AND more than 80% of alnLen
*/
proc freq data = blast_mo17_reads ;
tables flag_alnLen_50per * flag_qlen_500 ;
run;
    /*
flag_alnLen_50per     flag_qlen_500

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |      1 |     13 |     14
         |   5.00 |  65.00 |  70.00
         |   7.14 |  92.86 |
         |  50.00 |  72.22 |
---------+--------+--------+
       1 |      1 |      5 |      6
         |   5.00 |  25.00 |  30.00
         |  16.67 |  83.33 |
         |  50.00 |  27.78 |
---------+--------+--------+
Total           2       18       20
            10.00    90.00   100.00

18 of the transcripts are ge 500bp - of these 18, most (13) are less than 50% of alnLen

*/
data blast_cnts_mo17_reads ;
length sampleID $16.;
set cnts3_: ;
run ;


data  blast_mo17_reads_filter1 ;
set  blast_mo17_reads ;
where flag_alnLen_50per = 1 and flag_qlen_500 = 1 ;
run ;  

data  blast_mo17_reads_filter2 ;
set  blast_mo17_reads ;
where flag_alnLen_80per = 1 and flag_qlen_2000 = 1 ;
run ;




