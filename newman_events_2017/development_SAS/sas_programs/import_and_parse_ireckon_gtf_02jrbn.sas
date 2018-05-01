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

%macro importIR(sample);

proc import datafile="/mnt/store/event_sandbox/ireckon/&sample./result.gtf"
     out=&sample._ireckon dbms=tab replace;
     getnames=no; guessingrows=170000;
run;

/* Parse results GTF */

data &sample._ireckon2;
  set &sample._ireckon;
  length sample_id $4.;
  length gene_id $100.;
  length transcript_id $100.;
  format exon_number best32. ;
  format rpkm best32. ;
  format frac1 best32. ;
  format conf_lo best32. ;
  format conf_hi best32. ;
  format frac2 best32. ;
  format cov best32. ;
  sample_id="&sample.";
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

%mend;
%importIR(NSC1);
%importIR(NSC2);

/* Stack and flag  novel transcripts, bin all transcripts, and count */

data flag_xs;
  set NSC1_ireckon2 NSC2_ireckon2;
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
run;

data ir_type;
  set flag_xs;
  keep gene_id  transcript_id flag_novel_unspliced flag_novel_isoform
       flag_novel_ir flag_known_isoform;
run;

data npc1 npc2;
  set flag_xs;
  if sample_id="NSC1" then output npc1;
  else if sample_id="NSC2" then output npc2;
  keep gene_id transcript_id;
run;

proc sort data=ir_type nodup;
  by gene_id transcript_id;
proc sort data=npc1 nodup;
  by gene_id transcript_id;
proc sort data=npc2 nodup;
  by gene_id transcript_id;
run;

data ir_npc1_2;
  merge npc1 (in=in1) npc2 (in=in2);
  by gene_id transcript_id;
  if in1 then flag_in_npc1=1; else flag_in_npc1=0;
  if in2 then flag_in_npc2=1; else flag_in_npc2=0;
run;

data ir_npc1_2_flags;
  merge ir_npc1_2 (in=in1) ir_type (in=in2);
  by gene_id transcript_id;
run;

proc freq data=ir_npc1_2_flags noprint;
  tables flag_in_npc1*flag_in_npc2*flaG_known_isoform*flag_novel_isoform*
         flag_novel_unspliced*flag_novel_IR / out=iso_by_type;
run;

proc print data=iso_by_type;
run;


/*

                         flag_      flag_
flag_in_    flag_in_     known_     novel_    flag_novel_      flag_
  npc1        npc2      isoform    isoform     unspliced     novel_IR    COUNT

    0           1          0          0            1             0        1423
    0           1          0          1            0             0        6051
    0           1          0          1            0             1        1040
    0           1          1          0            0             0         524
    0           1          1          0            0             1          82
    1           0          0          0            1             0        1605
    1           0          0          1            0             0        7896
    1           0          0          1            0             1        1520
    1           0          1          0            0             0         697
    1           0          1          0            0             1         138
    1           1          0          0            1             0        2105
    1           1          0          1            0             0         184
    1           1          0          1            0             1          23
    1           1          1          0            0             0         140
    1           1          1          0            0             1          10


*/
