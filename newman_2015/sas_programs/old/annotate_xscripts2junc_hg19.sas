/* annotating genes to non-annotated junctions */

libname splice '/home/jrbnewman/McLab/junction_annotations/sas_data/';

libname splice2 '/home/jrbnewman/Desktop/sas_temp/';


/* split junctions in annotated and non-annotated */

data junctions_w_annot junctions_wo_annot oops;
   set splice.junctions_annotated_hg19;
   if flag_junction_annotated=1 then output junctions_w_annot;
   else if flag_junction_annotated=0 then output junctions_wo_annot;
   else output oops; *0 obs!;
run;

/* split non-annotated into junction_id and exonA/exonB*/

data junctions_exonA;
   set junctions_wo_annot;
   keep junction_pos exonA;
   rename exonA=exon_id;
run;

data junctions_exonB;
   set junctions_wo_annot;
   keep junction_pos exonB;
   rename exonB=exon_id;
run;

/* catback together together, remove duplicate observations */

data junction_exons;
   set junctions_exonA junctions_exonB;
run;

proc sort data=junction_exons nodup;
   by junction_pos exon_id;
run;

/* merge in gene */

data exon2gene;
   set splice.donor_exons_hg19;
   keep exonA geneA_id;
   rename exonA=exon_id;
   rename geneA_id=gene_id;
run;

proc sort data=junction_exons;
   by exon_id;
run;

proc sort data=exon2gene;
   by exon_id;
run;

data junc_wo_annot_w_gene nojunc noexon;
   merge junction_exons (in=in1) exon2gene (in=in2);
   by exon_id;
   if in1 and in2 then output junc_wo_annot_w_gene;
   else if in1 then output noexon; *0 obs!;
   else output nojunc; *456425 obs, this is fine;
run;

/* Cat genes by junctions */

/* first remove duplicate observations */

data junc_wo_annot_w_gene2;
   set junc_wo_annot_w_gene;
   drop exon_id;
run;

proc sort data=junc_wo_annot_w_gene2 nodups;
   by junction_pos gene_id;
run;

*2164031 duplicates dropped, left with 1311853;

/* get counts first */

proc freq noprint data=junc_wo_annot_w_gene2;
   tables junction_pos / out=junc_wo_annot_w_gene_cnt;
run;

proc sort data=junc_wo_annot_w_gene_cnt;
  by descending count;
run;
*max=2 genes per junction; *flag multigene later!;


data junction_wo_annot_gene_cat; 
  array genes[2] $ 36;

  retain genes1-genes2;

  set junc_wo_annot_w_gene2;
  by junction_pos;
  
  if first.junction_pos then do;
     call missing(of genes1-genes2);
     records = 0;
  end;

  records + 1;
  genes[records]=gene_id;
  if last.junction_pos then output;
run;

  *clean up the output file;

data junction_wo_annot_gene_cat2;
  set junction_wo_annot_gene_cat;
  length genes_cat2 $ 75;
  rename records= num_genes2;
         genes_cat2= catx("|", OF genes1-genes2);
  drop genes1-genes2 gene_id ;
  run;

/* merge back with non-annotated junctions */

proc sort data=junctions_wo_annot;
   by junction_pos;
run;

proc sort data=junction_wo_annot_gene_cat2;
   by junction_pos;
run;

data junc_wo_annot_w_cat_gene oops1 oops2;
    merge junctions_wo_annot (in=in1) junction_wo_annot_gene_cat2 (in=in2);
    by junction_pos;
    if in1 and in2 then output junc_wo_annot_w_cat_gene;
    else if in1 then output oops1; *0 obs!;
    else output oops2; *0 obs!;
run;

/* replace cat_genes and num_genes */

data junctions_wo_annot2;
    set junc_wo_annot_w_cat_gene;
    num_xscripts=0;
    num_genes=num_genes2;
    genes_cat=genes_cat2;
    if num_genes2>1 then flag_multigene=1;
    else flag_multigene=0;
    drop num_genes2 genes_cat2;
run;

/* cat annotated and non-annotated back together */

data junctions_annotated2;
   set junctions_w_annot junctions_wo_annot2;
run;

/* make permenant */

data splice2.junctions_annotated_all;
   set junctions_annotated2;
run;

/* Clean up */
    proc datasets nolist;
        delete junctions_w_annot oops junctions_wo_annot junctions_exonA junctions_exonB
        junction_exons exon2gene junc_wo_annot_w_gene nojunc noexon junc_wo_annot_w_gene2
        junc_wo_annot_w_gene_cnt junction_wo_annot_gene_cat junction_wo_annot_gene_cat2
        junc_wo_annot_w_cat_gene oops1 oops2 junctions_wo_annot2 junctions_annotated2
        ;
        run;
        quit;

