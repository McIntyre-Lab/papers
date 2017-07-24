ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';


/* Import PacBio BLAST results for detected junctions from expressed genes without
   multigene components */

    data WORK.PACBIO_JUNC    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile
'!MCLAB/event_analysis/analysis_output/blast_output/blast_events_to_pacbio_min_length_50_nomulti.tsv'
delimiter='09'x MISSOVER DSD lrecl=32767 ;
       informat event_id $173. ;
       informat pacbio_id $40. ;
       informat perc_identity best32. ;
       informat length best32. ;
       informat mismatch best32. ;
       informat gapopen best32. ;
       informat query_start best32. ;
       informat query_stop best32. ;
       informat ref_start best32. ;
       informat ref_stop best32. ;
       informat evalue best32. ;
       informat bitscore best32. ;
       format event_id $173. ;
       format pacbio_id $40. ;
       format perc_identity best12. ;
       format length best12. ;
       format mismatch best12. ;
       format gapopen best12. ;
       format query_start best12. ;
       format query_stop best12. ;
       format ref_start best12. ;
       format ref_stop best12. ;
       format evalue best12. ;
       format bitscore best12. ;
    input
                event_id $
                pacbio_id $
                perc_identity
                length
                mismatch
                gapopen
                query_start
                query_stop
                ref_start
                ref_stop
                evalue
                bitscore
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;

/* Add in annotations -- need flags for annot junc and IR, event size, gene_id */

data event_annot;
  set evspl.splicing_events_annot_refseq;
  keep event_id flag_junction_annotated flag_intron_retention gene_id event_size;
run;

proc sort data=event_annot;
  by event_id;
proc sort data=pacbio_junc;
  by event_id;
run;

data pacbio_blast_w_annot;
   merge pacbio_junc (in=in1) event_annot (in=in2);
   by event_id;
   if in1 and in2;
run;

data pacbio_blast_w_annot2;
   set pacbio_blast_w_annot;
   length pacbio_gene_id $10.;
   pacbio_gene_id=catt("PB.",scan(pacbio_id,2,"."));
run;

data pb2rs_gene;
   set event.pacbio2refseq_gene;
   keep gene_id pacbio_gene_id;
run;

proc sort data=pb2rs_gene nodup;
  by pacbio_gene_id;
proc freq data=pb2rs_gene noprint;
  tables pacbio_gene_id / out=gene_count;
proc sort data=gene_count;
  by descending count;
proc print data=gene_count (obs=1);
run; *15 RefSeq genes per PacBio gene max;

data cat_gene;
  array gene[15] $15.;
  retain gene1-gene15;
  set pb2rs_gene;
  by pacbio_gene_id;
  if first.pacbio_gene_id then do;
   call missing(of gene1-gene15);
   records=0;
  end;
  records + 1;
  gene[records] = gene_id;
  if last.pacbio_gene_id then output;
run;

data cat_gene2;
  set cat_gene;
  length refseq_gene_id $200.;
  refseq_gene_id=catx("|", OF gene1-gene2);
  keep pacbio_gene_id refseq_gene_id;
run;


proc sort data=pacbio_blast_w_annot2;
  by pacbio_gene_id;
run;

data pacbio_blast_w_annot3;
  merge pacbio_blast_w_annot2 (in=in1) cat_gene2 (in=in2);
  by pacbio_gene_id;
  if in1 and in2 then do;
        flag_pb_gene_has_refseq=1;
        output; end;
  else if in1 then do;
        flag_pb_gene_has_refseq=0;
        output; end;
run;

/* Check to make sure that none of the unannotated features are only going to genes without refseq IDs */

proc freq data=pacbio_blast_w_annot3;
   tables flag_pb_gene_has_refseq*flag_junction_annotated;
run;

/*
  flag_pb_gene_has_refseq
            flag_junction_annotated

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |   1825 |   5192 |   7017
           |   1.17 |   3.33 |   4.50
           |  26.01 |  73.99 |
           |  15.13 |   3.61 |
  ---------+--------+--------+
         1 |  10235 | 138582 | 148817
           |   6.57 |  88.93 |  95.50
           |   6.88 |  93.12 |
           |  84.87 |  96.39 |
  ---------+--------+--------+
  Total       12060   143774   155834
               7.74    92.26   100.00

Nope, good!
*/

/* Remove hits with gaps or mismatches */


data drop_partial;
  set pacbio_blast_w_annot3;
  if perc_identity < 95 then delete;
  if mismatch > 0 then delete;
  if gapopen > 0 then delete;
run;

proc freq data=drop_partial;
   tables flag_pb_gene_has_refseq*flag_junction_annotated;
run;


/*
 flag_pb_gene_has_refseq
           flag_junction_annotated

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |    897 |   4816 |   5713
          |   0.60 |   3.24 |   3.84
          |  15.70 |  84.30 |
          |  11.45 |   3.42 |
 ---------+--------+--------+
        1 |   6936 | 136193 | 143129
          |   4.66 |  91.50 |  96.16
          |   4.85 |  95.15 |
          |  88.55 |  96.58 |
 ---------+--------+--------+
 Total        7833   141009   148842
              5.26    94.74   100.00


So far so good!!
*/

/* Make permenant so I can process next */

data event.blast_dtct_junc2pb_nomult;
  set drop_partial;
run;


