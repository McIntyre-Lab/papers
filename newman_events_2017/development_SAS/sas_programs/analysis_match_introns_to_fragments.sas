ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* For each intron, find the 5' and 3' flanking fusions and fragments

   I already have the 5' and 3' flanking fusions in the introns info dataset,
   so now I want to find the corresponding fragment */

data gene_exp;
  set event.flag_gene_expressed;
  where flag_gene_expressed=1;
  keep gene_id;
run;

data intron_coord;
   set mm10.mm10_introns_from_fusions_si;
run;

proc sort data=gene_exp;
  by gene_id;
proc sort data=intron_coord;
  by gene_id;
run;

data intron_coord2;
   merge intron_coord (in=in1) gene_exp (in=in2);
   by gene_id;
   if in1 and in2;
   drop gene_id;
run; *154834 introns from expressed genes;


data frag_coord;
   set mm10.mm10_exon_fragment_flagged;
   keep fragment_id fragment_start fragment_end fusion_id chr;
run;

data frag_3prime2;
   set frag_coord;
   intron_stop=fragment_start-1;
   rename fusion_id=fusion_id_3prime fragment_id=fragment_id_3prime
   fragment_start=fragment_start_3prime fragment_end=fragment_end_3prime;
run;

data frag_5prime2;
   set frag_coord;
   intron_start=fragment_end;
   rename fusion_id=fusion_id_5prime fragment_id=fragment_id_5prime
   fragment_start=fragment_start_5prime fragment_end=fragment_end_5prime;
run;


proc sort data=frag_5prime2;
  by fusion_id_5prime intron_start;
proc sort data=intron_coord2;
  by fusion_id_5prime intron_start;
run;

data intron_coord_w_frag_5prime no_5prime;
   merge intron_coord2 (in=in1) frag_5prime2 (in=in2);
   by fusion_id_5prime intron_start;
   if in1 and in2 then output intron_coord_w_frag_5prime;
   else if in1 then output no_5prime;
run;

proc sort data=frag_3prime2;
  by fusion_id_3prime intron_stop;
proc sort data=intron_coord_w_frag_5prime;
  by fusion_id_3prime intron_stop;
run;

data intron_coord_w_frag_3prime no_3prime;
   merge intron_coord_w_frag_5prime (in=in1) frag_3prime2 (in=in2);
   by fusion_id_3prime intron_stop;
   if in1 and in2 then output intron_coord_w_frag_3prime;
   else if in1 then output no_3prime;
run;

/* Check to make sure that the 5' and 3' intron coordinates match fragments */

data intron_check;
   set intron_coord_w_frag_3prime;
   if fragment_end_5prime=intron_start then flag_intron_frag5_okay=1; else flag_intron_frag5_okay=0;
   if fragment_start_3prime=intron_stop+1 then flag_intron_frag3_okay=1; else flag_intron_frag3_okay=0;
run;

proc freq data=intron_check;
   tables flag_intron_frag5_okay flag_intron_frag3_okay;
run;
*all is okay!;

/* Make permenant */

data event.mm10_introns_to_fragment_fusion;
  set intron_coord_w_frag_3prime;
  keep fusion_id_5prime fusion_id_3prime intron_id chr intron_start intron_stop
       fragment_id_5prime fragment_id_3prime;
run;



