ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* RefSeq hits for unannotated junctions : 
Need to figure out if hits represent redundant sequences, or a bug in the creation code

If event and target transcript are on different chromosomes, then flag the sequence as redundant
If on the same chromosome, check if the gene of the event and the gene of the target transcript overlap
If so, and the event is an IR event, check if the complete IR coordinate overlaps with an exonic region
If the event is a junction, then check if the event sequence overlap with junctions from the target gene
*/


/* Check if events and target transcripts are on the same chromosome */



data refseq_hits;
   set event.refseq_blast_hits_w_annot;
   where flag_junction_annotated=0 and flag_unannotated_blast_hit=1;
   keep event_id blast_only_transcript_id;
run;

data refseq_hits2;
   length transcript_id $15.;
   set refseq_hits;
    do i=1 by 1 while(scan(blast_only_transcript_id,i,"|") ^= "");
         transcript_id=scan(blast_only_transcript_id,i,"|");
         output;
         end;
   keep event_id transcript_id;
run;

data event_chrom;
  set evspl.splicing_events_annot_refseq;
  if transcript_id ne "" then flag_junction_annotated=1;
  keep event_id gene_id chr strand feature1_start feature1_stop feature2_start feature2_stop
       flag_junction_annotated flag_intron_retention;
run;

data xs_chrom;
   length transcript_id2 $15.;
   set mm10.mm10_exons_w_info;
   do i=1 by 1 while(scan(transcript_id,i,"|") ^= "");
     transcript_id2=scan(transcript_id,i,"|");
     output;
     end;
   keep transcript_id2 gene_id chrom;
   rename transcript_id2=transcript_id chrom=transcript_chr gene_id=transcript_gene;
run;

proc sort data=event_chrom;
   by event_id;
proc sort data=refseq_hits2 nodup;
   by event_id;
run;

data refseq_hits_w_event_chr;
  merge refseq_hits2 (in=in1) event_chrom (in=in2);
   by event_id ;
   if in1 and in2;
run;

proc sort data=refseq_hits_w_event_chr nodup;
   by transcript_id;
proc sort data=xs_chrom nodup;
   by transcript_id;
run;

data refseq_hits_w_chroms;
   merge refseq_hits_w_event_chr (in=in1) xs_chrom (in=in2);
   by transcript_id;
   if in1 and in2;
run;

data flag_chrom;
   set refseq_hits_w_chroms;
   if chr=transcript_chr then flag_same_chrom=1;
   else flag_same_chrom=0;
run;

proc freq data=flag_chrom;
   tables flag_same_chrom;
run;

/*
                                              Cumulative    Cumulative
  flag_same_chrom    Frequency     Percent     Frequency      Percent
  --------------------------------------------------------------------
                0         418       27.66           418        27.66
                1        1093       72.34          1511       100.00

*/



/* Gene check: are events and target transcripts from the same gene? */



data flag_gene;
  set flag_chrom;
  if transcript_gene=gene_id then flag_same_gene=1;
  else flag_same_gene=0;
run;


proc freq data=flag_gene;
   tables flag_same_gene flag_same_gene*flag_same_chrom;
run;

/*
                              The FREQ Procedure

                                                Cumulative    Cumulative
     flag_same_gene    Frequency     Percent     Frequency      Percent
     -------------------------------------------------------------------
                  0        1325       87.69          1325        87.69
                  1         186       12.31          1511       100.00

           flag_same_gene     flag_same_chrom

           Frequency|
           Percent  |
           Row Pct  |
           Col Pct  |       0|       1|  Total
           ---------+--------+--------+
                  0 |    418 |    907 |   1325
                    |  27.66 |  60.03 |  87.69
                    |  31.55 |  68.45 |
                    | 100.00 |  82.98 |
           ---------+--------+--------+
                  1 |      0 |    186 |    186
                    |   0.00 |  12.31 |  12.31
                    |   0.00 | 100.00 |
                    |   0.00 |  17.02 |
           ---------+--------+--------+
           Total         418     1093     1511
                       27.66    72.34   100.00

*/


/* Make permenant */

data event.rfsq_unannot_hits_flag_chr_gene;
   set flag_gene;
run;

