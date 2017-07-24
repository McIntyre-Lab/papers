
ods listing; ods html close;
libname hg19 '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';

/* Count xscripts per exon fragment, flag unique/common/constit, create gene-level flags for gene level analysis */

/* If fragment goes to one gene:
	Unique: 1 transcript
	Constit: all transcripts
	Common: many transcripts
   If fragment goes to multiple genes:
	Common: many transcripts
	Constit: all transcripts for all genes
*/

/* Calculate the total number of transcripts per gene */


data xs2gene;
  set hg19.hg19_aceview_exons;
  length transcript_id2 $60.;  
    do i=1 by 1 while(scan(transcript_id,i,"|") ^= " ");
       transcript_id2=scan(transcript_id,i,"|");
       output;
       end;
  keep transcript_id2 gene_id;
  rename transcript_id2=transcript_id;
run;

proc sort data=xs2gene nodup;
  by gene_id transcript_id;
proc freq data=xs2gene noprint;
  tables gene_id / out=num_xs_per_gene;
run;

data frag2gene;
  set hg19.hg19_aceview_exon_fragment_info;
  length gene_id2 $36.;  
    do i=1 by 1 while(scan(gene_id,i,"|") ^= " ");
       gene_id2=scan(gene_id,i,"|");
       output;
       end;
  keep fragment_id gene_id2;
  rename gene_id2=gene_id;
run;

proc sort data=frag2gene;
  by gene_id;
proc sort data=num_xs_per_gene;
  by gene_id;
run;

data frags_to_flag;
  merge frag2gene (in=in1) num_xs_per_gene (in=in2);
  by gene_id;
  if in1 and in2;
run;

proc sort data=frags_to_flag nodup;
   by fragment_id;
proc means data=frags_to_flag noprint;
   by fragment_id;
   var count;
   output out=frag_max_xs sum(count)=max_xs_per_frag;
run;

data frag_info;
  set hg19.hg19_aceview_exon_fragment_info;
  keep fragment_id num_xscripts_per_fragment num_genes_per_fragment ;
run;

proc sort data=frag_info;
   by fragment_id;
proc sort data=frag_max_xs;
   by fragment_id;
run;

data frag_info2;
  merge frag_info (in=in1) frag_max_xs (in=in2);
  by fragment_id;
  if in1 ;
run;


/* If fragment goes to one gene:
	Unique: 1 transcript
	Constit: all transcripts
	Common: many transcripts
   If fragment goes to multiple genes:
	Common: many transcripts
	Constit: all transcripts for all genes
*/

data hg19_exon_fragment_flagged;
  set frag_info2;
  if num_genes_per_fragment = 1 then do;
    if num_xscripts_per_fragment = 1 then flag_unique=1;
    else flag_unique=0;
    if num_xscripts_per_fragment = max_xs_per_frag then flag_constitutive=1;
    else flag_constitutive=0;
    if num_xscripts_per_fragment > 1 and num_xscripts_per_fragment < max_xs_per_frag then flag_common=1;
    else flag_common=0;
    end;
  if num_genes_per_fragment > 1 then do;
    flag_unique=0;
    if num_xscripts_per_fragment = max_xs_per_frag then flag_constitutive=1;
    else flag_constitutive=0;
    if num_xscripts_per_fragment > 1 and num_xscripts_per_fragment < max_xs_per_frag then flag_common=1;
    else flag_common=0;
    end;
  keep fragment_id flag_unique flag_common flag_constitutive;
run;

proc freq data=hg19_exon_fragment_flagged noprint ;
   tables flag_unique*flag_common*flag_constitutive / out=check_flags;
run;

/* 20,697 fragments constitutive
   454,727 fragments common
   192,009 fragments unique
    77,096 fragments unique and constitutive (ie, single transcript genes);
*/

proc sort data=hg19.hg19_aceview_exon_fragment_info;
   by fragment_id;
proc sort data=hg19_exon_fragment_flagged;
   by fragment_id;
run;

data hg19.hg19_exon_fragment_flagged;
  merge hg19.hg19_aceview_exon_fragment_info (in=in1) hg19_exon_fragment_flagged (in=in2);
  by fragment_id;
  if in1 and in2;
  if num_genes_per_fragment gt 1 then flag_multigene=1;
  else flag_multigene=0;
run;



proc freq data=hg19.hg19_exon_fragment_flagged noprint ;
   tables flag_multigene*flag_unique*flag_common*flag_constitutive / out=check_flags;
run;

proc print data=check_flags;
run;

/* 
  17,893 fragments constitutive, single gene
  403,918 fragments common, single gene
  192,009 fragments unique, single gene
   77,096 fragments unique and constitutive, single gene (single-isoform genes)
    2,804 fragments constitutive, multiple genes
   50,809 fragments common, multiple genes
*/

