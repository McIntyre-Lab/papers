ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Add donor and acceptor fusions to unannotated junctions */

*subset only the expressed genes -- this will make this run quicker;
data exp_genes;
  set event.flag_gene_expressed;
  where flag_gene_expressed=1;
  keep gene_id;
run;

data junc2exons;
   set evspl.splicing_events_annot_refseq;
   where flag_junction_annotated=0 and flag_intron_retention=0;
   keep event_id feature1_id feature2_id gene_id;
run;

proc sort data=exp_genes;
  by gene_id;
proc sort data=junc2exons;
  by gene_id;
run;

data junc2exons_exp;
  merge junc2exons (in=in1) exp_genes (in=in2);
  by gene_id;
  if in1 and in2;
run;

*stack feature IDs;

data junc2exons_exp2;
   length donor_num $250.;
   length acceptor_num $250.;
   set junc2exons_exp;
   donor_num=scan(feature1_id,2,":");
   acceptor_num=scan(feature2_id,2,":");
run;

data stack_donors;
   length donor_exon_id $15.;
   set junc2exons_exp2;
   do i=1 by 1 while(scan(donor_num,i,"|") ^= "");
        donor_exon_id=catx(":",gene_id,scan(donor_num,i,"|"));
        output;
        end;
run;


data stack_acceptors;
   length acceptor_exon_id $15.;
   set stack_donors;
   do i=1 by 1 while(scan(acceptor_num,i,"|") ^= "");
        acceptor_exon_id=catx(":",gene_id,scan(acceptor_num,i,"|"));
        output;
        end;
   keep event_id gene_id donor_exon_id acceptor_exon_id;
run;

proc sort data=stack_acceptors nodup;
   by event_id donor_exon_id acceptor_exon_id;
run;

data donor_fus;
   set mm10.mm10_refseq_fusion_si_info_v2;
   keep fusion_id exon_id;
   rename fusion_id=donor_fusion_id exon_id=donor_exon_id;
run;

data acceptor_fus;
   set mm10.mm10_refseq_fusion_si_info_v2;
   keep fusion_id exon_id;
   rename fusion_id=acceptor_fusion_id exon_id=acceptor_exon_id;
run;

proc sort data=stack_acceptors;
   by donor_exon_id;
proc sort data=donor_fus;
   by donor_exon_id;
run;

data junc_donor_fus;
  merge stack_acceptors (in=in1) donor_fus (in=in2);
  by donor_exon_id;
  if in1 and in2;
run;

proc sort data=junc_donor_fus;
   by acceptor_exon_id;
proc sort data=acceptor_fus;
   by acceptor_exon_id;
run;

data junc_acceptor_fus;
   merge junc_donor_fus (in=in1) acceptor_fus (in=in2);
   by acceptor_exon_id;
   if in1 and in2;
run;

/* Make permenant -- next step will be flagging if fusions and events are on/off */

data unannot_jnc_w_flanking_fus;
  set junc_acceptor_fus;
  keep event_id gene_id donor_fusion_id acceptor_fusion_id;
run;

proc sort data=unannot_jnc_w_flanking_fus nodup;
   by event_id gene_id donor_fusion_id acceptor_fusion_id;
run;

data event.unannot_jnc_w_flanking_fus;
  set unannot_jnc_w_flanking_fus;
run;
