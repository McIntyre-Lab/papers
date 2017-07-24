ods listing; ods html close;
libname hg19 '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';

/* Make an exon fragment annotation file for isoform pipeline.
   I need (NO CONCATENTATIONS):
   fragment_id = identifier of fragment
   chrom = chromosome
   start = start position of fragment
   end = end position of fragment
   exon_id = exon associated with fragment (this will be stacked)
   gene_id = gene associated with exon
   gene_name = gene symbol associated with exon (make this the same as gene_id for now)
   flag_multigene = fragment associated with multiple genes
*/

data frags;
   set hg19.hg19_aceview_exon_fragment_info;
   length exon_id2 $39.;
   if num_genes_per_fragment > 1 then flag_multigene=1;
   else flag_multigene=0;
   do i=1 by 1 while(scan(exon_id,i,"|") ^= " ");
      exon_id2=scan(exon_id,i,"|");
      output;
      end;
   keep fragment_id exon_id2 fragment_start fragment_end chr flag_multigene;
   rename exon_id2=exon_id;
run;

data exon2gene;
   set hg19.hg19_aceview_exons;
   keep gene_id exon_id;
run;

proc sort data=frags;
   by exon_id;
proc sort data=exon2gene;
   by exon_id;
run;

data frag2exon2gene;
  merge frags (in=in1) exon2gene (in=in2);
  by exon_id;
  if in1 and in2;
run;

data hg19.hg19_fragment2exon2gene;
  retain gene_id gene_name fragment_id exon_id flag_multigene chr fragment_start fragment_end;
  length gene_name $11.;
  set frag2exon2gene;
  gene_name=gene_id;
  rename chr=chrom fragment_start=start fragment_end=end;
run;

