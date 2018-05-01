/* libraries */

libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
ods listing; ods html close;

/* Import ASTAlavista results GTF */

proc import datafile="!MCLAB/event_analysis/analysis_output/astalavista/asta_test_complete.gtf"
     out=asta_complete dbms=tab replace;
     getnames=no; guessingrows=100;
run;

data asta_complete2;
   set asta_complete;
   length transcript_id $100.;
   length gene_id $100.;
   length flanks $30.;
   length structure $30.;
   length splice_chain $30.;
   transcript_id=compress(scan(scan(VAR9,1,';'),2,'"'));
   gene_id=compress(scan(scan(VAR9,2,';'),2,'"'));
   flanks=compress(scan(scan(VAR9,3,';'),2,'"'));
   structure=compress(scan(scan(VAR9,4,';'),2,'"'));
   splice_chain=compress(scan(scan(VAR9,5,';'),2,'"'));
   drop VAR2 VAR3 VAR6 VAR8 VAR9;
   rename VAR1=chr VAR4=start VAR5=stop VAR7=strand;
run;

/* Merge with catalog junctions */

data juncs;
  set evspl.splicing_events_annot_refseq;
  where event_type="exon_junction";
    if gene_id in ('21419','217378','234725','24058','435802','666528',
                 '70375','71893','74140','74190','12564','12722');
  keep chr event_id gene_id transcript_id strand feature1_stop feature2_start flag_junction_annotated
        feature1_start feature2_stop
       flag_exonskip flag_alt_donor flag_alt_acceptor;
  rename transcript_id=transcript_cat gene_id=gene_cat ;
run;

proc sort data=asta_complete2;
   by chr start stop strand;
proc sort data=juncs;
   by chr feature1_stop feature2_start strand;
run;


data asta2cat asta_nocat cat_noasta;
  merge asta_complete2 (in=in1) juncs2 (in=in2);
  by chr start stop strand;
  if in1 and in2 then output asta2cat;
  else if in1 then output asta_nocat;
  else output cat_noasta;
run;

