ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';


/* Test calculating PSI for gene Rbm7 (geneID 67010) and compare to MISO output */


* Get all junctions in gene ;

data junc;
  set evspl.splicing_events_annot_refseq;
  where gene_id="71276" and flag_intron_retention=0;
  keep event_id gene_id flag_exonskip feature1_id feature2_id event_size;
run;

data exonskip;
  set evspl.skipped_exon_list_rfsq;
run;

data junc_apn;
   set refseq.rfsq_counts_by_splicing;
run;

 data splicing_on;
    set event.splicing_on_apn_gt0;
    keep event_id flag_splicing_on;
 run;

* Get all fragments in Rbm7 ;

data exon2frag;
   set mm10.mm10_exon_fragment_flagged;
   where gene_id ? "71276";
   keep fragment_id exon_id;
run;

data frag_on;
   set event.fragments_on_apn_gt0;
   keep fragment_id flag_fragment_on;
run;

data frag_apn;
   set event.mm10_refseq_fragment_counts;
   where sample_id ? "NSC";
   keep sample_id fragment_id region_depth apn;
run;

/* Exon fragment prep */

data exon2frag2;
   length exon_id2 $36.;
   set exon2frag;
   do i=1 by 1 while(scan(exon_id,i,"|") ^= "");
       exon_id2=scan(exon_id,i,"|");
       output;
   end;
   keep fragment_id exon_id2;
   rename exon_id2=exon_id;
run;

proc sort data=frag_apn;
   by fragment_id;
proc means data=frag_apn noprint;
   by fragment_id;
   var region_depth apn;
   output out=mean_frag_apn mean=;
run;

proc sort data=mean_frag_apn;
  by fragment_id;
proc sort data=exon2frag2;
  by fragment_id;
proc sort data=frag_on;
  by fragment_id;
run;

data frag_info;
  merge exon2frag2 (in=in1) mean_frag_apn frag_on;
  by fragment_id;
  if in1;
run;

/* Junction fragment prep */

proc sort data=junc_apn;
   by event_id;
proc means data=junc_apn noprint;
   by event_id;
   var apn;
   output out=mean_junc_apn mean=;
run;

proc sort data=splicing_on;
   by event_id;
proc sort data=junc;
   by event_id;
proc sort data=mean_junc_apn;
   by event_id;
run;

data junc_info;
  merge junc (in=in1) splicing_on mean_junc_apn;
  by event_id;
  if in1 ;
run;

* relate donor/acceptor features to exon and to exon skipping annotations

data junc_info2;
  set junc_info;
  length donors $10.;
  length acceptors $10.;
  donors=scan(feature1_id,2,":");
  acceptors=scan(feature2_id,2,":");
  keep event_id gene_id event_size flag_splicing_on apn donors acceptors;
run;

data junc_donors;
   length donor_num $3.;
   length donor_exon_id $100.;
   set junc_info2;
   do i=1 by 1 while(scan(donors,i,"|") ^= "");
     donor_num= scan(donors,i,"|");
     donor_exon_id=catx(":",gene_id,donor_num);
     output;
   end;
   drop i donor_num donors;
run;

data junc_acceptors;
   length acceptor_num $3.;
   length acceptor_exon_id $100.;
   set junc_donors;
   do i=1 by 1 while(scan(acceptors,i,"|") ^= "");
     acceptor_num= scan(acceptors,i,"|");
     acceptor_exon_id=catx(":",gene_id,acceptor_num);
     output;
   end;
   drop i acceptor_num acceptors;
run;

data junc_w_id;
  length junction_id $100.;
  set junc_acceptors;
  junction_id=catx("|",donor_exon_id,acceptor_exon_id);
  region_depth=apn*event_size;
run;

proc sort data=junc_w_id;
   by junction_id;
proc sort data=exonskip;
   by junction_id;
run;

data junc_w_exonskip_annot;
   merge junc_w_id (in=in1) exonskip;
   by junction_id;
   if in1;
run;

/* Parse junction data */

data donor;
   set junc_w_exonskip_annot;
   keep donor_exon_id event_id apn region_depth flag_splicing_on;
   rename donor_exon_id=exon_id event_id=fragment_id  flag_splicing_on=flag_fragment_on;
run;

data acceptor;
   set junc_w_exonskip_annot;
   keep acceptor_exon_id event_id apn region_depth flag_splicing_on;
   rename acceptor_exon_id=exon_id event_id=fragment_id  flag_splicing_on=flag_fragment_on;
run;

data exonskipping;
   set junc_w_exonskip_annot;
   where flag_exonskip=1;
   keep skipped_exon_id event_id apn region_depth flag_splicing_on;
   rename skipped_exon_id=exon_id event_id=fragment_id flag_splicing_on=flag_fragment_on;
run;

proc sort data=donor nodup;
   by fragment_id exon_id;
proc sort data=acceptor nodup;
   by fragment_id exon_id;
proc sort data=exonskipping nodup;
   by fragment_id exon_id;
run;

/* Stack all data */

data stack_features;
   length component $2.;
   set donor (in=in1) acceptor (in=in2) exonskipping (in=in3) frag_info (in=in4);
   if in1 then component="B1";
   if in2 then component="B2";
   if in3 then component="C";
   if in4 then component="A";
   drop _TYPE_ _FREQ_;
run;

/* For each exon, sum all in A, B1, B2, and C */

proc sort data=stack_features;
  by exon_id component;
proc means data=stack_features noprint;
  by exon_id component;
  var apn region_depth;
  output out=sum_components sum=;
run;

/* Transpose */

proc sort data=sum_components;
   by exon_id component;
proc transpose data=sum_components out=exon_by_comp_apn;
   by exon_id;
   var apn;
   id component;
run;

proc transpose data=sum_components out=exon_by_comp_depth;
   by exon_id;
   var region_depth;
   id component;
run;

/* Calc PSI */

data psi_from_apn;
   set exon_by_comp_apn;
   if A=. or B1=. or B2=. or C=. then flag_calc_psi=0;
   else do;
     flag_calc_psi=1;
     psi_apn=(A+B1+B2)/(A+B1+B2+C)*100;
     psi_apn_junc=(B1+B2)/(B1+B2+C)*100;
     end;
run;


data psi_from_depth;
   set exon_by_comp_depth;
   if A=. or B1=. or B2=. or C=. then flag_calc_psi=0;
   else do;
     flag_calc_psi=1;
     psi_depth=(A+B1+B2)/(A+B1+B2+C)*100;
     psi_depth_junc=(B1+B2)/(B1+B2+C)*100;
     end;
run;




