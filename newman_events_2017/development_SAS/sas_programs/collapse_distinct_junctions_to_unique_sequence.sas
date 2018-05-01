/* Libraries */

libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';

/* Import distinct junctions FASTA and collapse on sequence */

    data WORK.JUNC_FA    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile '!MCLAB/event_analysis/references/mm10_refseq_distinct_junc_coord.fa'
delimiter='09'x MISSOVER DSD lrecl=32767 ;
       informat junction_id $27. ;
       informat seq $81. ;
       format junction_id $27. ;
       format seq $81. ;
    input
                junction_id $
                seq $
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;

/* Replace lowercase characters with uppercase */

data junc_fa2;
  set junc_fa;
  length seq2 $81.;
  seq2=upcase(seq);
run;

proc sort data=junc_fa2;
   by seq2;
proc freq data=junc_fa2 noprint;
   tables seq2 / out=seq_count;
proc sort data=seq_count;
  by descending count;
run;
*2951862 distinct coordinates, 2927315 unique sequences (strand specific);

data junc_fa_check;
  set junc_fa2;
  if seq2="GTGGTGTTGCTGACTTGCCCTAGGTCACCTGAGTTAAAAGGGTTGCATTCCCTTTCCAGTATGGTCCATGGAAAACATGC";
run;


/* if sequence is on the minus strand, flip sequence and re-collapse */

data junc_fa3;
   set junc_fa2;
   length strand $1.;
   length seq3 $81.;
   length seq2_rev $81.;
   length seq3_tmp1 $81.;
   strand=scan(junction_id,4,":");
   if strand = "-" then do;
       seq2_rev = input(seq2,$reverj81.);
       seq3_tmp1 = compress(tranwrd(tranwrd(tranwrd(tranwrd(seq2_rev,"A","W"),"C","X"),"G","Y"),"T","Z"));
       seq3 = compress(tranwrd(tranwrd(tranwrd(tranwrd(seq3_tmp1,"W","T"),"X","G"),"Y","C"),"Z","A"));
       end;
   else seq3=seq2;
run;

  /* now count number of unique sequences (non-strand specific)*/
proc sort data=junc_fa3;
   by seq2 strand ;
proc freq data=junc_fa3 noprint;
   tables seq2*strand / out=seq_count;
proc sort data=seq_count;
  by descending count;
run;

proc freq data=seq_count noprint;
  tables seq2 / out=seq_count2;
proc sort data=seq_count2;
  by descending count;
run;

/* Okay, I DON'T need to flip the sequence. I am going to export the set of 2927315 unique junction sequences
   and map to these, forcing ONLY forward alignments */

proc sort data=junc_fa3;
   by seq2 ;
proc freq data=junc_fa3 noprint;
   tables seq2 / out=junc_uniq_seq;
run;

data junc_uniq_seq2;
  length unique_junc_id $25.;
  set junc_uniq_seq;
  unique_junc_id=catt("junction_",_n_);
  keep seq2 unique_junc_id;
run;

proc sort data=junc_uniq_seq2;
   by seq2;
proc sort data=junc_fa3;
   by seq2;
run;

data junc2uniq_seq;
  merge junc_fa3 (in=in1) junc_uniq_seq2 (in=in2);
  by seq2;
  if in1 and in2;
run;

/* Make permenant and output seq fa */

data evspl.mm10_refseq_junc2uniq_seq;
   set junc2uniq_seq;
   keep junction_id seq2 unique_junc_id;
run;

data junc_uniq_seq3;
  length seq_name $26.;
  set junc_uniq_seq2;
  seq_name=catt(">",unique_junc_id);
  keep seq_name seq2;
run;

proc export data=junc_uniq_seq3
     outfile="!MCLAB/event_analysis/references/mm10_unique_junction_seq.tsv"
     dbms=tab replace;
     putnames=no;
run;

  


