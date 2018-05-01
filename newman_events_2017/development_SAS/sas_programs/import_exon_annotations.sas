/* Import exon chunks, match to fusions, assign ID, save and export a BED */

libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname splicing '!MCLAB/conesa_pacbio/sas_data/splicing';
libname res '!MCLAB/isoform_reconstruction/sas_data/';

/* Import fusion info */
proc import datafile='!MCLAB/conesa_pacbio/created_files/conesa_pacbio_mm10_fusions_si.bed'
            out=fusions_bed
            dbms=tab replace; guessingrows=73221; getnames=no;
run;

data fusions_bed2;
   set fusions_bed;
   rename VAR1=chrom
          VAR2=fusion_start
          VAR3=fusion_end
          VAR4=fusion_id
          VAR5=score
          VAR6=strand
   ;
run;


proc import datafile='!MCLAB/conesa_pacbio/created_files/conesa_pacbio_mm10_fusions_si.tsv'
            out=fusions_info
            dbms=tab replace; guessingrows=91610;
run;

/* Import chromosome list */

proc import datafile='!MCLAB/conesa_pacbio/chrom_list.txt'
            out=chrom_list
            dbms=csv replace; getnames=no;
run;


data chrom_list2;
   set chrom_list;
   rename VAR1=chrom;
run;


/* Import gene and transcript info */
%include '!MCLAB/conesa_pacbio/sas_programs/macros/iterdataset.sas';

%macro import_genes(chrom);

proc import datafile="!MCLAB/conesa_pacbio/created_files/conesa_pacbio_mm10_gene_list_&chrom..csv"
            out=genes
            dbms=csv replace; guessingrows=1500;
run;

data gene_list_&chrom.;
   length chrom $9.;
   set genes;
   chrom="&chrom.";
run;
%mend;

%iterdataset(dataset=chrom_list2, function=%nrstr(%import_genes(&chrom);));

data conesa.pacbio_gene_list;
   set gene_list_chr: ;
run;


/* Import exon chunks */
proc import datafile='!MCLAB/conesa_pacbio/created_files/conesa_pacbio_mm10_exon_chunks.csv'
            out=exon_chunks
            dbms=csv replace; guessingrows=92833;
run;

/* Assembe fusion info dataset */

proc sort data=fusions_bed2;
   by fusion_id;
proc sort data=fusions_info;
   by fusion_id;
run;

data fusions_si_info;
   merge fusions_bed2 (in=in1) fusions_info (in=in2);
   by fusion_id;
   if in1 and in2;
run;

proc sort data=splicing.pacbio_exons;
   by chrom exon_id;
proc sort data=fusions_si_info;
   by chrom exon_id;
run;

data conesa.fusions_si_info;
   merge fusions_si_info (in=in1) splicing.pacbio_exons (in=in2);
   by chrom exon_id;
   if in1 and in2;
run;

/* Merge fusion info with exon chunks */

data fusions_data;
   set conesa.fusions_si_info;
   keep fusion_id chrom fusion_start fusion_end exon_id;
run;

data exon_chunks2;
   length exon_id2 $10.;
   length chunk_coord $25.;
   set exon_chunks;
   chunk_coord=catx(":",chrom,chunk_start,chunk_end);
   do i=1 by 1 while(scan(exon_id,i,'|') ^=' ');
   exon_id2=scan(exon_id,i,'|');
   output;
   end;
   drop i group_start group_stop;
   rename exon_id=exon_cat exon_id2=exon_id;
run;


proc sort data=fusions_data nodup;
     by chrom exon_id ;
proc sort data=exon_chunks2;
     by chrom exon_id ;
run;

data chunk2fusion no_chunk no_fus;
   merge fusions_data (in=in1) exon_chunks2 (in=in2);
   by chrom exon_id ;
   if in1 and in2 then output chunk2fusion;
   else if in1 then output no_chunk;
   else output no_fus;
run;

data chunk2fusion2;
   set chunk2fusion;
   drop exon_id;
   rename exon_cat=exon_id;
run;

proc sort data=chunk2fusion2 nodup;
   by chunk_coord;
run;

*92832;




proc sort data=fusions_data nodup;
     by chrom exon_id ;
proc sort data=exon_chunks2;
     by chrom exon_id ;
run;

data chunk2fusion no_chunk no_fus;
   merge fusions_data (in=in1) exon_chunks2 (in=in2);
   by chrom exon_id ;
   if in1 and in2 then output chunk2fusion;
   else if in1 then output no_chunk;
   else output no_fus;
run;


/* Merge fusion info with exon chunks */

data fusions_data;
   set conesa.fusions_si_info;
   keep fusion_id chrom fusion_start fusion_end;
   rename fusion_start=group_start fusion_end=group_stop;
run;

proc sort data=fusions_data nodup;
   by chrom group_start group_stop;
proc sort data=exon_chunks;
   by chrom group_start group_stop;
run;


data chunk2fusion no_chunk no_fus;
   merge fusions_data (in=in1) exon_chunks (in=in2);
   by chrom group_start group_stop;
   if in1 and in2 then output chunk2fusion;
   else if in1 then output no_chunk;
   else output no_fus;
run;

proc sort data=chunk2fusion;
   by fusion_id group_start group_stop;
run;

data chunk_add_id;
  retain fus_chunk_num;
  length exonchunk_id $15.;
  set chunk2fusion;
  by fusion_id;
  if first.fusion_id then fus_chunk_num=1;
  else fus_chunk_num=fus_chunk_num+1;
  exonchunk_id=catx(":",fusion_id,fus_chunk_num);
run;

/* Make permenant */

data conesa.exon_chunk2fusion;
   set chunk_add_id;
run;

/* Export BED file for exon chunks */

data chunk_bed;
  retain chrom chunk_start chunk_end exonchunk_id;
  set conesa.exon_chunk2fusion;
  keep exonchunk_id chunk_start chunk_end chrom;
run;

proc sort data=chunk_bed;
   by chrom chunk_start chunk_end exonchunk_id;
run;

proc export data=chunk_bed outfile='!MCLAB/conesa_pacbio/created_files/conesa_pacbio_mm10_exon_chunks.bed' dbms=tab replace; putnames=no;
run;


