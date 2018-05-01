ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';


/* Match PabBio gene to RefSeq gene */


data pb2refseq;
  length pacbio_gene_id $10.;
   set event.pacbio2refseq_id;
  pacbio_gene_id=catt("PB.",scan(pacbio_id,2,"."));
  keep pacbio_id pacbio_gene_id transcript_id match_type;
run;

data xs2gene;
  set event.feature2xs2gene;
  keep gene_id transcript_id;
run;

/* Assign Refseq geneID to pacbio transcripts */

proc sort data=xs2gene nodup;
  by transcript_id gene_id;
proc sort data=pb2refseq nodup;
  by transcript_id pacbio_id;
run;

data gene2pb;
  merge xs2gene (in=in1) pb2refseq (in=in2);
  by transcript_id;
  if in1 and in2 then output;
  else if in1 then delete;
  else if in2 then do;
      gene_id="NovelGene";
      output; end;
run;

/* Make permenant */

data event.pacbio2refseq_gene;
  set gene2pb;
run;

/* Import list of PB transcripts with type */

proc import datafile="!MCLAB/event_analysis/references/pacbio_transcript_list.txt"
    out=pb_iso dbms=tab replace;
    getnames=no;
run;

data pb_iso2;
  set pb_iso;
  pacbio_gene_id=catt("PB.",scan(VAR1,2,"."));
  keep pacbio_gene_id VAR1 VAR2 VAR3;
  rename VAR1=pacbio_id VAR2=type VAR3=class;
run;

proc sort data=pb_iso2;
  by pacbio_gene_id class;
proc freq data=pb_iso2 noprint;
  by pacbio_gene_id;
  tables class / out=pb_iso_by_class;
run;

proc transpose data=pb_iso_by_class out=pb_iso_class_sbys;
   by pacbio_gene_id;
   id class;
   var count;
run;

/* Merge in Refseq ID */

data pb2rs;
   set gene2pb;
   keep gene_id pacbio_gene_id;
run;

proc sort data=pb2rs nodup;
  by pacbio_gene_id gene_id;
proc sort data=pb_iso_class_sbys;
  by pacbio_gene_id;
run;

data pb_iso_w_refseq;
  merge pb2rs (in=in1) pb_iso_class_sbys (in=in2) ;
  by pacbio_gene_id;
  if gene_id="" then gene_id="NovelGene";
run;

data pb_iso_w_refseq2;
   set pb_iso_w_refseq;
   array change _numeric_;
            do over change;
            if change=. then change=0;
            end;
   drop _NAME_ _LABEL_;
run;

proc sort data=pb_iso_w_refseq2 nodup;
   by pacbio_gene_id;
run;

/* Make permenant */

data event.pacbio_gene_to_refseq;
   set pb_iso_w_refseq2;
run;



