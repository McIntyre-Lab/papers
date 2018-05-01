ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/splicing/sas_data';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';

/* Import iReckon results and process 

   iReckon output is a GTF file consisting of transcripts and exons
   Each transcript has the following attributes:
         gene_id	iReckon assigns this, so the actual "geneID" is not useful
         transcript_id	Gleaned from annotation, else assigned by iReckon
         RPKM		RPKM for transcript
         frac		Not sure? Fraction of gene expression?
         conf_lo	?
         conf_hi	?
         frac		?
         cov		?	

   Each exon has the following additional attribute:
         exon_number	number of exon in isoform
*/

proc import datafile="/mnt/store/event_sandbox/ireckon/NSC1/result.gtf"
     out=NSC1_ireckon dbms=tab replace;
     getnames=no; guessingrows=170000;
run;

/* Parse results GTF */

data NSC1_ireckon2;
  set NSC1_ireckon;
  length gene_id $100.;
  length transcript_id $100.;
  format exon_number best32. ;
  format rpkm best32. ;
  format frac1 best32. ;
  format conf_lo best32. ;
  format conf_hi best32. ;
  format frac2 best32. ;
  format cov best32. ;
  if VAR3="transcript" then do;
    gene_id=compress(scan(scan(VAR9,1,";"),2,'"'));
    transcript_id=compress(scan(scan(VAR9,2,";"),2,'"'));
    rpkm=compress(scan(scan(VAR9,3,";"),2,'"'))+0;
    frac1=compress(scan(scan(VAR9,4,";"),2,'"'))+0;
    conf_lo=compress(scan(scan(VAR9,5,";"),2,'"'))+0;
    conf_hi=compress(scan(scan(VAR9,6,";"),2,'"'))+0;
    frac2=compress(scan(scan(VAR9,7,";"),2,'"'))+0;
    cov=compress(scan(scan(VAR9,8,";"),2,'"'))+0;
  end;

  else if VAR3="exon" then do;
    gene_id=compress(scan(scan(VAR9,1,";"),2,'"'));
    transcript_id=compress(scan(scan(VAR9,2,";"),2,'"'));
    exon_number=compress(scan(scan(VAR9,3,";"),2,'"'))+0;
    rpkm=compress(scan(scan(VAR9,4,";"),2,'"'))+0;
    frac1=compress(scan(scan(VAR9,5,";"),2,'"'))+0;
    conf_lo=compress(scan(scan(VAR9,6,";"),2,'"'))+0;
    conf_hi=compress(scan(scan(VAR9,7,";"),2,'"'))+0;
    frac2=compress(scan(scan(VAR9,8,";"),2,'"'))+0;
    cov=compress(scan(scan(VAR9,9,";"),2,'"'))+0;
  end;
  else delete;
  drop VAR9;
  rename VAR1=chr VAR2=source VAR3=feature_type VAR4=start VAR5=stop VAR6=score
         VAR7=strand VAR8=frame;
  run;

/* Look at the distribution of RPKM */

proc univariate data=NSC1_ireckon2;
   where feature_type="transcript";
   var RPKM;
run;

/*
  N                       14584    Sum Weights              14584
  Mean               2.63823822    Sum Observations    38476.0661
  Std Deviation      46.3263206    Variance            2146.12798
  Skewness           64.0394938    Kurtosis            5416.07969
  Uncorrected SS     31398493.4    Corrected SS        31296984.4
  Coeff Variation     1755.9567    Std Error Mean      0.38360961

 Quantile          Estimate

 100% Max       4.31841E+03
 99%            3.65364E+01
 95%            3.06849E+00
 90%            8.67036E-01
 75% Q3         1.89917E-01
 50% Median     6.00721E-02
 25% Q1         1.61778E-02
 10%            2.77361E-04
 5%             1.17781E-04
 1%             3.75831E-05
 0% Min         1.00671E-05
*/

/* flag novel transcripts, bin all transcripts, and count */

data flag_and_bin;
  set NSC1_ireckon2;
  where feature_type="transcript";
  /* Flag transcripts */
  if index(transcript_id,"unspliced") > 0 then flag_novel_unspliced=1; else flaG_novel_unspliced=0;
  if index(transcript_id,"novel") > 0 then flag_novel_isoform=1; else flaG_novel_isoform=0;
  if index(transcript_id,"Intron") > 0 then flag_novel_IR=1; else flaG_novel_IR=0;
  if index(transcript_id,"NM_") > 0 
     or index(transcript_id,"NR_") > 0
     or index(transcript_id,"XM_") > 0
     or index(transcript_id,"XR_") > 0
     then flag_known_isoform=1; else flag_known_isoform=0;
  /* Bin by log RPKM : 0-1, 1-2, 2-4, 4+ */
  log_rpkm=log(rpkm+1);
  if log_rpkm = 0 then isoform_bin=0;
  else if log_rpkm < 1 then isoform_bin=1;
  else if log_rpkm < 2 then isoform_bin=2;
  else if log_rpkm < 4 then isoform_bin=3;
  else isoform_bin=4;
run;

proc freq data=flag_and_bin noprint;
  tables flaG_known_isoform*flag_novel_isoform*flag_novel_unspliced*flag_novel_IR / out=iso_by_type;
  tables isoform_bin / out=iso_by_bin;
  tables isoform_bin*flaG_known_isoform*flag_novel_isoform*flag_novel_unspliced*flag_novel_IR / out=iso_by_type_and_bin;
run;

proc print data=iso_by_type;
run;

/*
 flag_      flag_
 known_     novel_    flag_novel_      flag_
isoform    isoform     unspliced     novel_IR    COUNT

   0          0            1             0        3792
   0          1            0             0        8121
   0          1            0             1        1552
   1          0            0             0         952
   1          0            0             1         167

Most transcripts are apparently novel...
*/

proc print data=iso_by_bin;
run;

/*
 isoform_
    bin      COUNT    PERCENT

     1       13595    93.2186
     2         522     3.5793
     3         357     2.4479
     4         110     0.7543

Most transcripts have low expression
*/

proc print data=iso_by_type_and_bin;
run;


/*
             flag_      flag_
isoform_     known_     novel_    flag_novel_      flag_
   bin      isoform    isoform     unspliced     novel_IR    COUNT    PERCENT

    1          0          0            1             0        3175    21.7704
    1          0          1            0             0        7873    53.9838
    1          0          1            0             1        1525    10.4567
    1          1          0            0             0         862     5.9106
    1          1          0            0             1         160     1.0971

    2          0          0            1             0         374     2.5645
    2          0          1            0             0          96     0.6583
    2          0          1            0             1          16     0.1097
    2          1          0            0             0          33     0.2263
    2          1          0            0             1           3     0.0206

    3          0          0            1             0         189     1.2959
    3          0          1            0             0         122     0.8365
    3          0          1            0             1           4     0.0274
    3          1          0            0             0          39     0.2674
    3          1          0            0             1           3     0.0206

    4          0          0            1             0          54     0.3703
    4          0          1            0             0          30     0.2057
    4          0          1            0             1           7     0.0480
    4          1          0            0             0          18     0.1234
    4          1          0            0             1           1     0.0069

*/

