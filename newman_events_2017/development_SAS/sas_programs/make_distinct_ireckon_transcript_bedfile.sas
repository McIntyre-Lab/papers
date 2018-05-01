/* import bedfile for iReckon transcripts and make a set of distinct isoforms to estimate */

    data WORK.IR_NSC1_BED    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile '!MCLAB/event_analysis/references/ireckon_isoforms_NSC1.bed' delimiter='09'x
MISSOVER DSD lrecl=32767 ;
       informat chr $5. ;
       informat start best32. ;
       informat stop best32. ;
       informat transcript_id $53. ;
       informat score $1. ;
       informat strand $1. ;
       informat start2 best32. ;
       informat stop2 best32. ;
       informat color $7. ;
       informat blocks best32. ;
       informat block_len $306. ;
       informat block_start $451. ;
       format chr $5. ;
       format start best12. ;
       format stop best12. ;
       format transcript_id $53. ;
       format score $1. ;
       format strand $1. ;
       format start2 best12. ;
       format stop2 best12. ;
       format color $7. ;
       format blocks best12. ;
       format block_len $306. ;
       format block_start $451. ;
    input
                chr $
                start
                stop
                transcript_id $
                score $
                strand $
                start2
                stop2
                color $
                blocks
                block_len $
                block_start $
    ;
    run;


    data WORK.IR_NSC2_BED    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile '!MCLAB/event_analysis/references/ireckon_isoforms_NSC2.bed' delimiter='09'x
MISSOVER DSD lrecl=32767 ;
       informat chr $5. ;
       informat start best32. ;
       informat stop best32. ;
       informat transcript_id $53. ;
       informat score $1. ;
       informat strand $1. ;
       informat start2 best32. ;
       informat stop2 best32. ;
       informat color $7. ;
       informat blocks best32. ;
       informat block_len $306. ;
       informat block_start $451. ;
       format chr $5. ;
       format start best12. ;
       format stop best12. ;
       format transcript_id $53. ;
       format score $1. ;
       format strand $1. ;
       format start2 best12. ;
       format stop2 best12. ;
       format color $7. ;
       format blocks best12. ;
       format block_len $306. ;
       format block_start $451. ;
    input
                chr $
                start
                stop
                transcript_id $
                score $
                strand $
                start2
                stop2
                color $
                blocks
                block_len $
                block_start $
    ;
    run;

/* Stack and see if there are duplicates */


data ir_nsc1_bed2;
  set ir_nsc1_bed;
  drop strand;
run;
data ir_nsc2_bed2;
  set ir_nsc2_bed;
  drop strand;
run;

proc sort data=ir_nsc1_bed2 nodup;
   by _all_;
proc sort data=ir_nsc2_bed2 nodup;
   by _all_;
run;

data ir_bed2;
  merge ir_nsc1_bed2 (in=in1) ir_nsc2_bed2 (in=in2) ;
  by _all_;
  if in1 and in2 then sample=3;
  else if in1 then sample=1;
  else sample=2;
run;

/* Add strand back in */

data ir_nsc12_bed;
  set ir_nsc1_bed ir_nsc2_bed;
  if strand = "" then delete;
run;

proc sort data=ir_bed2 nodup;
  by chr start stop transcript_id score start2 stop2 color blocks block_len block_start;
proc sort data=ir_nsc12_bed nodup;
  by chr start stop transcript_id score start2 stop2 color blocks block_len block_start;
run;

data ir_nsc12_bed2;
  set ir_nsc12_bed;
  by chr start stop transcript_id score start2 stop2 color blocks block_len block_start;
  if first.block_start then output;
run;


data ir_bed3;
  merge ir_bed2 (in=in1) ir_nsc12_bed2 (in=in2);
  by chr start stop transcript_id score start2 stop2 color blocks block_len block_start;
  if in1;
run;

data ir_bed4;
  set ir_bed3;
  length sample_transcript_id $100.;
  if sample=3 then sample_transcript_id=catx("_","Both",transcript_id);
  if sample=1 then sample_transcript_id=catx("_","NPC1",transcript_id);
  if sample=2 then sample_transcript_id=catx("_","NPC2",transcript_id);
  if score="" then score = ".";
  if strand="" then strand = "+";
  drop sample transcript_id;
run;


data ir_bed5;
 retain chr start stop sample_transcript_id score strand start2 stop2 color blocks block_len block_start;
 set ir_bed4;
run;

proc export data=ir_bed5 outfile="!MCLAB/event_analysis/references/ireckon_transcript_bed_union_npc1_npc2.bed"
  dbms=tab replace;
  putnames=no;
run;

data ir_xs2gene;
  length gene_id $100.;
  set ir_bed5;
  if index(sample_transcript_id,"NM_") > 0 or index(sample_transcript_id,"NR_") > 0 or
     index(sample_transcript_id,"XM_") > 0 or index(sample_transcript_id,"XR_") > 0 then
     gene_id=catx("_",scan(sample_transcript_id,2,"_"),scan(sample_transcript_id,3,"_"));
  else gene_id=scan(sample_transcript_id,2,"_");
  keep gene_id sample_transcript_id;
run;

proc export data=ir_xs2gene outfile="!MCLAB/event_analysis/references/ireckon_transcript_bed_union_npc1_npc2_gene2xs.txt"
  dbms=tab replace;
  putnames=no;
run;


/* extract sequence with bedtools, reimport and collapse sequences 


bedtools getfasta -fi ~/McLab/useful_mouse_data/mm10/fasta/mm10_for_bedtools_v2.fa -bed ireckon_transcript_bed_union_npc1_npc2.bed -fo ireckon_transcript_bed_union_npc1_npc2.tsv -tab -s -split -name
jrbnewman@ufgi-immigrans ~/McLab/event_analysis/references $ bedtools getfasta -fi ~/McLab/useful_mouse_data/mm10/fasta/mm10_for_bedtools_v2.fa -bed ireckon_transcript_bed_union_npc1_npc2.bed -fo ireckon_transcript_bed_union_npc1_npc2.tsv -tab -s -split -name

*/
