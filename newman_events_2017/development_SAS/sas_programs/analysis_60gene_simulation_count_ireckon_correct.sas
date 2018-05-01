/* Libraries */
libname event  '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* 60 gene simulation, iReckon check:
   what does it get right in terms of transcript retention? 

   I am going to use the "union" iReckon transcriptome here.
   I am not going to use the iReckon-to-Refseq BLAST results here, just those IDs
   that iReckon correctly assigns
*/

data ir;
  set event.polyester_60gene_sim_ireckon;
run;

/*   Because of the way that iReckon assigns gene IDs when none are provided (yet, I can't find
   info on this in the IR documentation...), I am going to assign RefSeq gene IDs to the iReckon IDs.
   Then assign the Refseq gene ID to each iR transcript with the same iR gene ID */

data ir_genes;
  set ir;
  length refseq_id $100.;
  if index(gene_id,"NM_") > 0 or index(gene_id,"NR_") > 0
     or index(gene_id,"XM_") > 0 or index(gene_id,"XR_") > 0
     then refseq_id=tranwrd(gene_id,"IntronRetention","");
  else refseq_id=tranwrd(transcript_id,"IntronRetention","");
  keep sample_id gene_id transcript_id refseq_id strand;
run;

data sim_list;
   set event.polyester_xs_list_60genes;
   keep transcript_id;
   rename transcript_id=refseq_id;
run;

data xs2gene;
   set event.feature2xs2gene;
   keep gene_id transcript_id;
   rename transcript_id=refseq_id gene_id=refseq_gene_id;
run;

proc sort data=ir_genes;
  by refseq_id;
proc sort data=xs2gene nodup;
  by refseq_id;
proc sort data=sim_list nodup;
  by refseq_id;
run;

data ir_gene_w_refseq;
  merge ir_genes (in=in1) xs2gene (in=in2) sim_list (in=in3);
  by refseq_id;
  if in3 then flag_simulated=1; else flag_simulated=0;
  if in1 and in2;
  keep sample_id gene_id refseq_gene_id flag_simulated strand;
run;

proc sort data=ir_gene_w_refseq nodup;
  by sample_id  strand gene_id refseq_gene_id;
run;

proc sort data=ir;
  by sample_id  strand gene_id;
run;

data ir_w_refseq_gene_id;
  merge ir (in=in1) ir_gene_w_refseq (in=in2);
  by sample_id strand gene_id;
  if in2 then flag_refseq_gene=1; else flag_refseq_gene=0;
  if in1 then output;
run;

/* Collapse transcripts -- some transcripts have the same RefSeq ID, but different iReckon gene_ids */

data ir_w_refseq_gene_id2;
  set ir_w_refseq_gene_id;
  length transcript_id2 $100.;
  if index(transcript_id,"NM_") > 0 or index(transcript_id,"NR_") > 0 or 
     index(transcript_id,"XM_") > 0 or index(transcript_id,"XR_") > 0 
     then transcript_id2=transcript_id;
  else transcript_id2=catx("|",gene_id, transcript_id);
  run;

proc sort data=ir_w_refseq_gene_id2;
  by transcript_id2 ;
run;

proc means data=ir_w_refseq_gene_id2 noprint;
  by transcript_id2;
  var flag_refseq_gene flag_simulated;
  output out=ir_collapse max=;
run;

data sim_list;
  length transcript_id $100.;
  format transcript_id $100.;
  informat transcript_id $100.;
  set event.polyester_xs_list_60genes;
  rename transcript_id=transcript_id2;
run;

proc sort data=sim_list;
  by transcript_id2;
proc sort data=ir_collapse;
  by transcript_id2;
run;

data ir_collapse_flag_60gn;
  merge ir_collapse (in=in1) sim_list (in=in2);
  by transcript_id2;
  if in1 then flag_detected=1; else flag_detected=0;
  if in2 then flag_simulated_xs=1; else flag_simulated_xs=0;
  if (index(transcript_id2,"NM_") > 0 or index(transcript_id2,"NR_") > 0 or 
     index(transcript_id2,"XM_") > 0 or index(transcript_id2,"XR_") > 0 )
     and index(transcript_id2,"Intron") = 0 and index(transcript_id2,"novel") = 0
          and index(transcript_id2,"unspliced") = 0
     then flag_refseq_xs=1; else flag_refseq_xs=0;
run;

proc freq data=ir_collapse_flag_60gn noprint;
  tables flag_detected*flag_simulated*flag_simulated_xs*flag_refseq_gene*flag_refseq_xs / out=ireckon_60gn_count;
run;

proc print data=ireckon_60gn_count;
run;

/*
                            flag_       flag_      flag_
  flag_       flag_      simulated_    refseq_    refseq_
detected    simulated        xs          gene        xs      COUNT

    0           .             1           .          1         76
    1           .             0           0          0         19
    1           1             0           1          0        546
    1           1             1           1          1        391


CORRECT			391		transcripts that were selected for simulation
RELATED_RS		0		RefSeq transcripts that were also estimated, but not from simulation.
						As all transcripts were simulated, this should be 0
RELATED_OTHER	546		Non-RefSeq transcripts estimated and belong to same genes as a Refseq gene
UNRELATED		19		All other transcripts
MISSING			76		transcripts simulated that were not estimated
*/

