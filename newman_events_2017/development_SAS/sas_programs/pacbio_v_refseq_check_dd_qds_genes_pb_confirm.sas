
ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';

/* Of the 93 DD and QDS mouse genes, how many are have DD in PacBio? */

data gene2keep;
   set event.mm10_flag_gene_dd_ds_exons_apn5;
   if flag_cell_by_fus_fdr05=1 and flag_gene_exon_dd=1;
   keep gene_id;
run;

data pb_flag_dd_iso;
   set event.pacbio_lr_read_counts;
   mean_NPC=mean((reads_FL_NPC1+reads_nonFL_NPC1),(reads_FL_NPC2+reads_nonFL_NPC2));
   mean_OPC=mean((reads_FL_OPC1+reads_nonFL_OPC1),(reads_FL_OPC2+reads_nonFL_OPC2));

   if mean_NPC>mean_OPC then flag_pb_higher_in_npc=1; else flag_pb_higher_in_npc=0;

   if mean_NPC<mean_OPC then flag_pb_higher_in_opc=1; else flag_pb_higher_in_opc=0;

   if mean_NPC=0 and mean_OPC=0 then do; flag_exp=0; flag_pacbio_dd=.; end;

   else do;
       flag_exp=1;
       if mean_NPC=0 or mean_OPC=0 then flag_pacbio_dd=1;
       else if abs(mean_NPC-mean_NPC) > 5 then flag_pacbio_dd=1;
       else flag_pacbio_dd=0;
       end;
run;

data add_pb_gene;
   set pb_flag_dd_iso;
   length pacbio_gene_id $15.;
   pacbio_gene_id=catt("PB.",scan(pacbio_id,2,"."));
run;

proc sort data=add_pb_gene;
   by pacbio_gene_id;
proc means data=add_pb_gene noprint;
   by pacbio_gene_id;
   var flag_exp flag_pacbio_dd flag_pb_higher_in_npc flag_pb_higher_in_opc;
   output out=pb_gene_flag_dd max=;
run;

data pb2refseq;
   set event.pacbio2refseq_gene_nomulti;
   keep gene_id pacbio_gene_id;
run;

proc sort data=pb2refseq nodup;
  by pacbio_gene_id gene_id;
proc sort data=pb_gene_flag_dd;
   by pacbio_gene_id;
run;

data pb2rs_flag_dd;
   merge pb2refseq (in=in1) pb_gene_flag_dd (in=in2);
   by pacbio_gene_id;
   if in1 and in2;
run;

proc sort data=pb2rs_flag_dd;
   by gene_id;
proc means data=pb2rs_flag_dd  noprint;
   by gene_id;
   var flag_exp flag_pacbio_dd flag_pb_higher_in_npc flag_pb_higher_in_opc;
   output out=pb2rs_gene_flag_dd max=;
run;

proc sort data=pb2rs_gene_flag_dd;
  by gene_id;
proc sort data=gene2keep;
  by gene_id;
run;

data gene2keep2;
  merge gene2keep (in=in1) pb2rs_gene_flag_dd (in=in2);
  by gene_id;
  if in1 and in2;
run; 

proc freq data=gene2keep2 noprint;
   tables flag_exp*flag_pacbio_dd*flag_pb_higher_in_npc*flag_pb_higher_in_opc /out=gene_count;
run;

proc print data=gene_count;
run; *64 genes in PB are also confirmed;


/*
                flag_     flag_pb_    flag_pb_
               pacbio_     higher_     higher_
   flag_exp       dd       in_npc      in_opc     COUNT

       1          0           0           1         12
       1          0           1           0         17
       1          0           1           1          4

       1          1           0           1         10
       1          1           1           0         11
       1          1           1           1         10

*/

