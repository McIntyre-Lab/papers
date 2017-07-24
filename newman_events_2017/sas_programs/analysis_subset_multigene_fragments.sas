libname event '!MCLAB/event_analysis/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';

 /* Extract the list of mutligene fragments that are detected. I am going to limit myself to fragments that are assigned to two genes so that I can make easier predictions later */

data frags_on;
  set event.flag_fragment_on;
  where flag_fragment_nsc_on=1;
  keep fragment_id;
run;

* I am going to pull only the fragments assigned to 2 genes for simplicity;
data multigene_frags;
  set mm10.mm10_exon_fragment_flagged;
  where num_genes_per_fragment=2;
  keep fragment_id exon_id transcript_id gene_id flag_unique flag_common flag_constitutive;
run;

proc sort data=frags_on;
  by fragment_id;
proc sort data=multigene_frags;
  by fragment_id;
run;

data event.multigene_frags_on;
  merge frags_on (in=in1) multigene_frags (in=in2);
  by fragment_id;
  if in1 and in2;
run;

*6958 2-gene fragments detected;

