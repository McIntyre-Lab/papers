/* Counts needed:
(1) PacBio NNC and NIC junctions
(2) Novel ID rate:
	Unannot EA vs STAR novel
	Unannot+Border EA vs STAR novel
(3) Redo counts, include ALL PB novel in PB-detected set
(4) Novel counts: NIC: STAR/EA vs PB
			STAR/EA: annotated as "unannot" in event catalog, vs NIC PB
                  NNC: STAR vs PB
			STAR: not in event catalog, vs NIC PB
(5) Annotated junctions:
	Restrict to ONLY PB junctions, and compare to PB junctions
		Also to genes with PB transcripts
	Count junctions detected by each method
	STAR/events vs PB
	STAR vs Events
*/

libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/splicing/sas_data';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';

/* Count PacBio NNC and NIC junctions */

data pb_junc;
   set event.catalog_pacbio_star_junctions;
   where flag_in_pacbio=1;
   if flag_junction_annotated=0 then flag_pb_junc_nic=1; else flag_pb_junc_nic=0;
   if flag_junction_annotated=. then flag_pb_junc_nnc=1; else flag_pb_junc_nnc=0;
run;

proc freq data=pb_junc noprint;
   tables  flag_pb_junc_nic*flag_pb_junc_nnc*flag_junction_annotated / out=pb_junc_type_cnt;
run;

proc print data=pb_junc_type_cnt;
run;


/*


 flag_pb_    flag_pb_    flag_junction_
 junc_nic    junc_nnc       annotated      COUNT    PERCENT

     0           0              1          68393    99.0557
     0           1              .           3691      .
     1           0              0            652     0.9443

*/


proc freq data=pb_junc noprint;
   tables flag_pb_junc_nic*flag_pb_junc_nnc*flag_donor_in_catalog*flag_acceptor_in_catalog
      /out=pb_junc_count;
run;

proc print data=pb_junc_count;
run;

/*

 flag_pb_    flag_pb_    flag_donor_    flag_acceptor_
 junc_nic    junc_nnc     in_catalog      in_catalog      COUNT

     0           0            1                1          68393
     0           1            0                0           1763
     0           1            0                1            962
     0           1            1                0            928
     0           1            1                1             38
     1           0            1                1            652

As expected, all 652 NIC junctions can be derived from existing annotations
38 NNC junctions have both their donor and acceptor site in the catalog, so what are these?
*/

data catalog_donor;
  set evspl.splicing_events_annot_refseq;
  where feature1_id ^? "intron";
  keep chr strand feature1_id feature1_stop;
  rename feature1_id=donor_id feature1_stop=donor_stop;
run;

data catalog_acceptor;
  set evspl.splicing_events_annot_refseq;
  where feature2_id ^? "intron";
  keep chr strand feature2_id feature2_start;
  rename feature2_id=acceptor_id feature2_start=acceptor_start;
run;

data pb_nnc_junc_check;
  set pb_junc;
  where flaG_pb_junc_nnc=1 and flag_donor_in_catalog=1 and flag_acceptor_in_catalog=1;
  keep chr donor_Stop acceptor_start strand;
run;


proc sort data=catalog_donor nodup;
   by chr strand donor_stop donor_id;
proc sort data=pb_nnc_junc_check nodup;
   by chr strand donor_stop;
run;

data pb_nnc_check2;
  merge pb_nnc_junc_check (in=in1) catalog_donor (in=in2);
  by chr strand donor_stop;
  if in1 and in2;
run;

proc sort data=catalog_acceptor nodup;
   by chr strand acceptor_start acceptor_id;
proc sort data=pb_nnc_check2 nodup;
   by chr strand acceptor_start;
run;

data pb_nnc_check3;
  merge pb_nnc_check2 (in=in1) catalog_acceptor (in=in2);
  by chr strand acceptor_start;
  if in1 and in2; *many-to-many merge, but is okay for this;
run;

data pb_nnc_check4;
  set pb_nnc_check3;
  donor_gene=scan(donor_id,1,":") + 0;
  acceptor_gene=scan(acceptor_id,1,":") + 0;
  if donor_gene=acceptor_gene then flaG_same_gene=1;
  else flag_same_gene=0;
run;

proc freq data=pb_nnc_check4;
  tables flag_same_gene;
run;

/*
                      
                                            Cumulative    Cumulative
 flaG_same_gene    Frequency     Percent     Frequency      Percent
 -------------------------------------------------------------------
              0          42      100.00            42       100.00        

*/



