ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';


/* Get quantification levels of Cav1 transcripts for each restricted transcriptome
   and the full PacBio transcriptome and make a table

Cav1: RefSeqID is 12389
      PacBioID is PB.5785
*/


/* Get Cav1 transcripts */

data cav1_xs;
   set event.feature2xs2gene_exp_only_nomulti;
   where gene_id="12389";
   keep transcript_id;
run;


/* All RefSeq */

data refseq_all;
   set event.rsem_refseq_all;
   mean_tpm_refseq=mean(tpm_nsc1,tpm_nsc2);
   keep transcript_id mean_tpm_refseq;
run;   

/* Any EA */

data event_any;
   set event.rsem_events_exp_any;
   mean_tpm_ea_any=mean(tpm_nsc1,tpm_nsc2);
   keep transcript_id mean_tpm_ea_any;
run;   

/* EA: 100% APN>0 */

data event_100p_apn0;
   set event.rsem_events_exp_100perc;
   mean_tpm_ea_100p_apn0=mean(tpm_nsc1,tpm_nsc2);
   keep transcript_id mean_tpm_ea_100p_apn0;
run; 

/* EA: 75% APN>0 */

data event_75p_apn0;
   set event.rsem_events_exp_75perc;
   mean_tpm_ea_75p_apn0=mean(tpm_nsc1,tpm_nsc2);
   keep transcript_id mean_tpm_ea_75p_apn0;
run; 

/* EA: 50% APN>0 */

data event_50p_apn0;
   set event.rsem_events_exp_50perc;
   mean_tpm_ea_50p_apn0=mean(tpm_nsc1,tpm_nsc2);
   keep transcript_id mean_tpm_ea_50p_apn0;
run; 

/* EA: 100% APN>5 */

data event_100p_apn5;
   set event.rsem_events_exp_100perc_apn5;
   mean_tpm_ea_100p_apn5=mean(tpm_nsc1,tpm_nsc2);
   keep transcript_id mean_tpm_ea_100p_apn5;
run; 

/* EA: 75% APN>5 */

data event_75p_apn5;
   set event.rsem_events_exp_75perc_apn5;
   mean_tpm_ea_75p_apn5=mean(tpm_nsc1,tpm_nsc2);
   keep transcript_id mean_tpm_ea_75p_apn5;
run; 

/* EA: 50% APN>5 */

data event_50p_apn5;
   set event.rsem_events_exp_50perc_apn5;
   mean_tpm_ea_50p_apn5=mean(tpm_nsc1,tpm_nsc2);
   keep transcript_id mean_tpm_ea_50p_apn5;
run; 

/* Merge EA and RefSeq data */

proc sort data=cav1_xs nodup;
   by transcript_id;
proc sort data=refseq_all;
   by transcript_id;
proc sort data=event_any;
   by transcript_id;
proc sort data=event_100p_apn0;
   by transcript_id;
proc sort data=event_75p_apn0;
   by transcript_id;
proc sort data=event_50p_apn0;
   by transcript_id;
proc sort data=event_100p_apn5;
   by transcript_id;
proc sort data=event_75p_apn5;
   by transcript_id;
proc sort data=event_50p_apn5;
   by transcript_id;
run;

data cav1_est_refseq;
   merge cav1_xs (in=in1) refseq_all event_any event_100p_apn0
          event_75p_apn0  event_50p_apn0 event_100p_apn5
          event_75p_apn5  event_50p_apn5;
   by transcript_id;
   if in1;
run;

/* PacBio transcriptome */

data pacbio_counts;
   set event.rsem_pacbio_all;
   where transcript_id ? "PB.5785.";
   mean_tpm_pacbio=mean(tpm_nsc1,tpm_nsc2);
   keep transcript_id mean_tpm_pacbio;
   rename transcript_id=pacbio_id;
run;

/* PacBio-to-Refseq */

data pb2refseq;
   set event.pacbio2refseq_id;
   keep pacbio_id transcript_id;
run;

proc sort data=pacbio_counts;
   by pacbio_id;
proc sort data=pb2refseq;
  by pacbio_id;
run;

data pacbio_w_Rs;
  merge pacbio_counts (in=in1) pb2refseq (in=in2);
   by pacbio_id;
  if in1;
run;

/* Merge RefSeq/EA and PAcBio */

proc sort data=pacbio_w_rs;
   by transcript_id;
proc sort data=cav1_est_refseq;
   by transcript_id;
run;

data cav1_est_all;
  merge pacbio_w_rs cav1_est_refseq;
  by transcript_id;
run;

/* Print table */

proc print data=cav1_est_all;
run;

/*

                                                                    mean_     mean_tpm_
                         mean_tpm_    transcript_     mean_tpm_    tpm_ea_     ea_100p_
    Obs    pacbio_id       pacbio     id                refseq       any         apn0

     1     PB.5785.3        0.865                         .           .           .
     2     PB.5785.5       14.260                         .           .           .
     3     PB.5785.6        0.215     NM_001243064      10.600      12.64       12.930
     4     PB.5785.1       29.765     NM_007616          8.390      10.53       11.125
     5     PB.5785.2        4.320     XM_006504974      13.970      16.84       18.090
     6     PB.5785.4        0.340     XM_006504975       0.415       0.49        0.540
     7                       .        XM_006504976       0.345       0.45        0.545

           mean_tpm_    mean_tpm_    mean_tpm_    mean_tpm_    mean_tpm_
            ea_75p_      ea_50p_      ea_100p_     ea_75p_      ea_50p_
    Obs       apn0         apn0         apn5         apn5         apn5

     1         .            .            .            .            .
     2         .            .            .            .            .
     3       12.775       12.720       14.200       13.070       12.790
     4       10.805       10.705       13.075       11.765       10.790
     5       17.080       16.955       20.200       17.890       17.505
     6        0.510        0.505         .            .           0.525
     7        0.470        0.455         .            .           0.495



*/
