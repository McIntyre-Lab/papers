ods listing; ods html close;
libname hg19 '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';


/* Exon fragment annotations:
   Get genes, transcripts, coordinations
   Number of genes, transcripts

   Add number of fragments to exons and fusions */

/* I want to make exon fragments relative to fusions AND exons, and base the fragment ID off fusion_id
   e.g. F10000_SI:1, to denote fusion "F10000_SI" fragment 1, since I think this makes the most sense
   since exons are grouped into fusions, which are then split into fragments. This also helps with
   relatability between our internal annotations */

* Relate fragment to fusion -- fusionID will become the "exon fragment group" identifier ;

data fus2exon;
  set hg19.hg19_aceview_fusions_si_info;
  keep fusion_id exon_id;
run;

data fus_coord;
  set hg19.hg19_aceview_fusions_si_bed;
  keep fusion_id fusion_start fusion_stop chr;
run;

proc sort data=fus2exon;
  by fusion_id;
proc sort data=fus_coord;
  by fusion_id;
run;

data fus2exon2;
  merge fus2exon (in=in1) fus_coord (in=in2);
  by fusion_id;
  if in1 and in2;
run;

data frags;
   set hg19.hg19_aceview_exon_fragments;
   length exon_id2 $39.;  
   if exon_id="" then output;
   else do i=1 by 1 while(scan(exon_id,i,"|") ^= " ");
       exon_id2=scan(exon_id,i,"|");
       output;
       end;
   drop exon_id;
   rename chrom=chr exon_id2=exon_id;
run;

proc sort data=fus2exon2;
   by chr exon_id;
proc sort data=frags;
   by chr exon_id;
run;

data frag2exon2fus nofus nofrag;
  merge frags (in=in1) fus2exon2 (in=in2);
  by chr exon_id;
  if in1 and in2 then output frag2exon2fus;
  else if in1 then output nofus; *0 obs!;
  else output nofrag; *0 obs!;
run;

* Check: are the fusion coordinates the same as the group coordinates?;

data flag_coord;
  set frag2exon2fus;
  if group_start=fusion_start then flag_same_start=1; else flag_same_start=0;
  if group_stop=fusion_stop then flag_same_stop=1; else flag_same_stop=0;
run;

proc freq data=flag_coord;
  tables flag_same_start*flag_same_stop;
run;
* Coordinates all the same, so this is good;

data frag2fus;
  set frag2exon2fus;
  keep fusion_id fragment_start fragment_end;
run;


proc sort data=frag2fus nodup;
  by fusion_id fragment_start fragment_end;
run;

data frag_add_id;
  length fragment_id $20.;
  set frag2fus;
  by fusion_id;
  if first.fusion_id then i=1;
  else i+1;
  fragment_id=catx(":",fusion_id,i);
  drop i;
run;

/* Merge with fus2exon2frag */

proc sort data=frag_add_id;
  by fusion_id fragment_start fragment_end;
proc sort data=frag2exon2fus;
  by fusion_id fragment_start fragment_end;
run;

data frag2exon2fus_w_id;
  merge frag_add_id (in=in1) frag2exon2fus (in=in2);
  by fusion_id fragment_start fragment_end;
  if in1 and in2;
run;

/* Merge in exon info to calculate number of exons, transcripts, genes per fragment */

data exon_info;
   set hg19.hg19_aceview_exons;
    keep chrom exon_id gene_id transcript_id;
   rename chrom=chr;
run;

proc sort data=exon_info;
  by chr exon_id;
proc sort data=frag2exon2fus_w_id;
  by chr exon_id;
run;

data frag2exon_info;
  merge frag2exon2fus_w_id (in=in1) exon_info (in=in2);
  by chr exon_id;
  if in1 and in2;
run;

data frag2exon_info2;
   set frag2exon_info;
   length transcript_id2 $60.;  
   if transcript_id="" then output;
   else do i=1 by 1 while(scan(transcript_id,i,"|") ^= " ");
       transcript_id2=scan(transcript_id,i,"|");
       output;
       end;
    keep chr fragment_id exon_id gene_id transcript_id2;
   rename transcript_id2=transcript_id;
run;

* Count number of exons per fragment and cat;

data frag2ex;
   set frag2exon_info2;
   keep fragment_id exon_id;
run;

proc sort data=frag2ex nodup;
   by fragment_id exon_id;
proc freq data=frag2ex noprint;
   tables fragment_id / out=exons_per_frag;
proc sort data=exons_per_frag;
   by descending count;
proc print data=exons_per_frag (obs=1);
run; *47 exons max;

data cat_ex; 
  array exon[47] $ 39.;
  retain exon1-exon47;
  set frag2ex;
  by fragment_id;
  if first.fragment_id then do;
     call missing(of exon1-exon47);
     records = 0;
  end;
  records + 1;
  exon[records]=exon_id;
  if last.fragment_id then output;
run;

data cat_ex2;
  set cat_ex;
  length cat_exon_id $1880.;
  cat_exon_id= catx("|", OF exon1-exon47);
  keep fragment_id records cat_exon_id;
  rename records= num_exons_per_fragment
         cat_exon_id=exon_id;
  run;

* Count number of xscripts per fragment and cat;

data frag2xs;
   set frag2exon_info2;
   keep fragment_id transcript_id;
run;

proc sort data=frag2xs nodup;
   by fragment_id transcript_id;
proc freq data=frag2xs noprint;
   tables fragment_id / out=xs_per_frag;
proc sort data=xs_per_frag;
   by descending count;
proc print data=xs_per_frag (obs=1);
run; *76 xscripts max;

data cat_xs; 
  array xs[76] $ 60.;
  retain xs1-xs76;
  set frag2xs;
  by fragment_id;
  if first.fragment_id then do;
     call missing(of xs1-xs76);
     records = 0;
  end;
  records + 1;
  xs[records]=transcript_id;
  if last.fragment_id then output;
run;

data cat_xs2;
  set cat_xs;
  length cat_transcript_id $2000.;
  cat_transcript_id= catx("|", OF xs1-xs76);
  keep fragment_id records cat_transcript_id;
  rename records= num_xscripts_per_fragment
         cat_transcript_id=transcript_id;
  run;


* Count number of genes per fragment and cat;

data frag2gene;
   set frag2exon_info2;
   keep fragment_id gene_id;
run;

proc sort data=frag2gene nodup;
   by fragment_id gene_id;
proc freq data=frag2gene noprint;
   tables fragment_id / out=genes_per_frag;
proc sort data=genes_per_frag;
   by descending count;
proc print data=genes_per_frag (obs=1);
run; *3 genes max;

data cat_gene; 
  array gene[3] $ 36.;
  retain gene1-gene3;
  set frag2gene;
  by fragment_id;
  if first.fragment_id then do;
     call missing(of gene1-gene3);
     records = 0;
  end;
  records + 1;
  gene[records]=gene_id;
  if last.fragment_id then output;
run;

data cat_gene2;
  set cat_gene;
  length cat_gene_id $120.;
  cat_gene_id= catx("|", OF gene1-gene3);
  keep fragment_id records cat_gene_id;
  rename records= num_genes_per_fragment
         cat_gene_id=gene_id;
  run;

/* Merge together */

data frag_info;
  set frag2exon2fus_w_id;
  keep fragment_id fragment_start fragment_end chr fusion_id;
run;


proc sort data=frag_info nodup;
  by fragment_id;
proc sort data=cat_ex2;
  by fragment_id;
proc sort data=cat_xs2;
  by fragment_id;
proc sort data=cat_gene2;
  by fragment_id;
run;

data hg19.hg19_aceview_exon_fragment_info;
   merge frag_info (in=in1) cat_ex2 (in=in2) cat_xs2 (in=in3) cat_gene2 (in=in4);
   by fragment_id;
   if in1 and in2 and in3 and in4;
run;


