libname events '!MCLAB/event_analysis/sas_data';

/* Subset 10000, 20000, 30000 random Refseq transcripts to use to
   simulate reads using Polyester for testing the
   accuracy of Event Analysis
   For each set I am going to have two groups, 6 replicates per group

   My only requirement is to select from genes without multigene regions */

data xscripts;
   set event.feature2xs2gene_nomulti;
   keep transcript_id;
run;

proc sort data=xscripts nodup;
  by transcript_id;
run;

/* Add random number */

data random;
   set xscripts;
   random_num1= rand("Uniform");
   random_num2= rand("Uniform");
   random_num3= rand("Uniform");
run;

proc sort data=random;
   by random_num1;
run;

data xs_list1;
  set random;
  if _n_ le 10000;
  keep transcript_id;
run;

proc sort data=random;
   by random_num2;
run;

data xs_list2;
  set random;
  if _n_ le 20000;
  keep transcript_id;
run;

proc sort data=random;
   by random_num3;
run;

data xs_list3;
  set random;
  if _n_ le 30000;
  keep transcript_id;
run;

/* Make permenant */

data event.polyester_xs_list_10k;
   set xs_list1;
run;

data event.polyester_xs_list_20k;
   set xs_list2;
run;

data event.polyester_xs_list_30k;
   set xs_list3;
run;

/* Export */

proc export data=xs_list1
     outfile="!MCLAB/event_analysis/references/mm10_refseq_10000_random_xscripts.txt"
     dmbs=tab replace;
     putnames=no;
run;

proc export data=xs_list2
     outfile="!MCLAB/event_analysis/references/mm10_refseq_20000_random_xscripts.txt"
     dmbs=tab replace;
     putnames=no;
run;

proc export data=xs_list3
     outfile="!MCLAB/event_analysis/references/mm10_refseq_30000_random_xscripts.txt"
     dmbs=tab replace;
     putnames=no;
run;


