/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* Import NIC junction to catalog BLAST and pull out hits for each transcript */

proc import datafile="!MCLAB/event_analysis/analysis_output/blast_output/NIC_junctions_to_check_BLAST_megablast.tsv"
     out=blast_results dbms=tab replace;
     getnames=no;
     guessingrows=max;
run;

data blast_results2;
  set blast_results;
  rename VAR1=junction_id
         VAR2=reference_junc
         VAR3=perc_identity
         VAR4=hit_length
         VAR5=num_mismatch
         VAR6=num_gapopen
         VAR7=query_start
         VAR8=query_stop
         VAR9=ref_start
         VAR10=ref_stop
         VAR11=evalue
         VAR12=bitscore;
run;

/* Drop hits that are shorter than the read length (56bp) */

data blast_length;
  set blast_results2;
  query_match_len=query_stop-query_start+1;
  ref_match_len=ref_stop-ref_start+1;
run;

data blast_len_gt56;
  set blast_length;
  if query_match_len < 56 then delete;
  if ref_match_len < 56 then delete;
  if hit_length < 56 then delete;
run;

data flag_hits;
  set blast_len_gt56;
  if junction_id=reference_junc then flag_hit_sameID=1; else flag_hit_sameID=0;
run;

proc freq data=flag_hits;
  tables flag_hit_sameID;
run;

/*

                                            Cumulative    Cumulative
flag_hit_sameID    Frequency     Percent     Frequency      Percent
--------------------------------------------------------------------
              1          19      100.00            19       100.00

Junctions don't map to other junctions in catalog, so is not due to sequence redundancy in alignments to catalog

Might need to look at alignments instead, and see exactly which reads are mapping to these junctions,
  and then check the sequence of the read against the sequence of the junctions...
*/


