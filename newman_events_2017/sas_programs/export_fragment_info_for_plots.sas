libname events '!MCLAB/event_analysis/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';

/* Merge in fragment annotations. I want to do the following:
  1. count number of fragments detected:
	a. all
	b. unique
	c. common, single gene
	d. constitutive, single gene
	e. common, multigene
	f. constitutive, multigene
  2. Distribution between APN and length: is there one?
  3. Distribution of number detected by fragment length
  4. Distribution of APN by categories from (1)
  5. Distribution of fragment length to commonality group
  6. Distribution between fragment length to number of transcripts
  7. Distribution between fragment length to number of exons

I think I can just output a single CSV that has the following (for NSCs and OLDs separately):
fragment_id, apn, fragment_length, detection flag, commonality_group, number of transcripts

Note: do the same for splicing events, but incorporate event type!!

*/

data fragment_annot;
  set mm10.mm10_exon_fragment_flagged;
  length commonality_group $30.;
  if flag_unique=1 then commonality_group="unique - single gene";
  else if flag_common=1 and flag_multigene=0 then commonality_group="common - single gene";
  else if flag_constitutive=1 and flag_multigene=0 then commonality_group="constitutive - single gene";
  else if flag_common=1 and flag_multigene=1 then commonality_group="common - multiple genes";
  else if flag_constitutive=1 and flag_multigene=1 then commonality_group="constitutive - multiple genes";
  else commonality_group="unknown - check annot";
  fragment_length=fragment_end-fragment_start;
    keep fragment_id num_exons_per_fragment num_xscripts_per_fragment fragment_length commonality_group;
run;

data flag_on;
  set events.flag_fragment_on;
run;


data frag_apn;
  length cell_type $3.;
  set events.mm10_refseq_fragment_counts;
  if sample_id in ('NSC1','NSC2') then cell_type="NSC";
  if sample_id in ('OLD1','OLD2') then cell_type="OLD";
  keep sample_id cell_type fragment_id apn;
run;

proc sort data=frag_apn;
  by cell_type fragment_id;
proc means data=frag_apn noprint;
  by cell_type fragment_id;
  var apn;
  output out=mean_apn_by_frag mean=;
run;

data nsc_apn old_apn;
  set mean_apn_by_frag;
  if cell_type="NSC" then output nsc_apn;
  if cell_type="OLD" then output old_apn;
run;

data nsc_apn2;
  set nsc_apn;
  rename apn=mean_apn_nsc;
run;

data old_apn2;
  set old_apn;
  rename apn=mean_apn_old;
run;

proc sort data=fragment_annot;
  by fragment_id;
proc sort data=flag_on;
  by fragment_id;
proc sort data=nsc_apn2;
  by fragment_id;
proc sort data=old_apn2;
  by fragment_id;
run;

data events.fragment_counts_w_info;
  merge fragment_annot (in=in1) flag_on (in=in2) nsc_apn2 (in=in3) old_apn2 (in=in4);
  by fragment_id;
  if in1 and in2 and in3 and in4;
  drop cell_type _TYPE_ _FREQ_;
run;

data nsc_data;
   set events.fragment_counts_w_info;
   drop mean_apn_old flag_fragment_old_on;
run;

data old_data;
   set events.fragment_counts_w_info;
   drop mean_apn_nsc flag_fragment_nsc_on;
run;

proc export data=nsc_data outfile="!MCLAB/event_analysis/analysis_output/annot_by_fragment_nsc.csv"
    dbms=csv replace;
run;


proc export data=old_data outfile="!MCLAB/event_analysis/analysis_output/annot_by_fragment_old.csv"
    dbms=csv replace;
run;



