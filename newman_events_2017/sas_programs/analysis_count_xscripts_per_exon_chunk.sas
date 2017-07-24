/* For each exon chunk, I want to count the number of exons, transcripts and genes*/

libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname splicing '!MCLAB/conesa_pacbio/sas_data/splicing';

/* Will have two files:
   exon chunks with all exons/transcripts/genes stack
   and the same but catted together */
   
data exon_info;
   set splicing.pacbio_exons;
   keep chrom exon_id transcript_id gene_id;
   run;
   
data exon_chunks;
   set conesa.exon_chunk2fusion;
   keep exonchunk_id chrom chunk_start chunk_end exon_id;
   run;
   
*stack exon_ids;

data exon_chunks2;
  set exon_chunks;
  length exon_id2 $15.;
  do i=1 by 1 while(scan(exon_id,i,'|') ^= ' ');
  exon_id2=scan(exon_id,i,'|');
  drop exon_id i;
  output;
  end;
  rename exon_id2=exon_id;
run;
  
proc sort data=exon_chunks2;
   by chrom exon_id;
proc sort data=exon_info;
   by chrom exon_id;
run;

data chunk2info;
  merge exon_chunks2 (in=in1) exon_info (in=in2);
  by chrom exon_id;
  if in1 and in2;
  run;
  
  
*stack transcript_ids;

data chunk2info2;
  set chunk2info;
  length transcript_id2 $15.;
  if transcript_id=' ' then do;
       transcript_id2=transcript_id;
       output;
       end;
  else 
    do i=1 by 1 while(scan(transcript_id,i,'|') ^= ' ');
  transcript_id2=scan(transcript_id,i,'|');
  drop transcript_id i;
  output;
  end;
  rename transcript_id2=transcript_id;
run;

%macro count_feat(feature);

data &feature.2chunk;
  set chunk2info2;
  keep exonchunk_id &feature._id;
  run;

proc sort data=&feature.2chunk nodup;
  by exonchunk_id &feature._id;
proc freq data=&feature.2chunk noprint;
  tables exonchunk_id / out=&feature.s_per_chunk;
run;


data &feature.s_per_chunk2;
   set &feature.s_per_chunk;
   rename count=&feature._count;
   run;

proc sort data=&feature.s_per_chunk2;
   by exonchunk_id;
   run;
%mend;      
  
* count number of exons per chunk;
%count_feat(exon);

* count number of transcripts per chunk;
%count_feat(transcript);

* count number of genes per chunk ;
%count_feat(gene);

proc sort data=chunk2info2;
   by exonchunk_id;
run;


data exon_chunk_w_info;
  merge chunk2info2 (in=in1) exons_per_chunk2 transcripts_per_chunk2 genes_per_chunk2;
  by exonchunk_id;
  if exon_count=1 then flag_chunk_multiexon=0;
  else flag_chunk_multiexon=1;
    if transcript_count=1 then flag_chunk_multitranscript=0;
  else flag_chunk_multitranscript=1;
    if gene_count=1 then flag_chunk_multigene=0;
  else flag_chunk_multigene=1;
  
  if in1 then output;
run;


*cat exons;

proc sort data=exons_per_chunk;
  by descending count;
  run;
  
proc print data=exons_per_chunk(obs=10);
run;


*cat transcripts;

proc sort data=transcripts_per_chunk;
  by descending count;
  run;
  
proc print data=transcripts_per_chunk(obs=10);
run;

*cat genes;
proc sort data=genes_per_chunk;
  by descending count;
  run;
  
proc print data=genes_per_chunk(obs=10);
run;

%macro cat_feat(feature,count);

data cat_&feature.;
  array &feature.[&count.] $ 10.;
  retain &feature.1-&feature.&count.;
  set &feature.2chunk;
  by exonchunk_id;
  if first.exonchunk_id then do;
     call missing(of &feature.1-&feature.&count.);
     records = 0;
  end;
  records + 1;
  &feature.[records]=&feature._id;
  if last.exonchunk_id then output;
run;

  *clean up the output file;
data cat_&feature.2;
  set cat_&feature.;
  length &feature._id2 $ 4000.;
  rename records= &feature._count;
         &feature._id2= catx("|", OF &feature.1-&feature.&count.);
  drop &feature.1-&feature.&count. &feature._id;
  rename &feature._id2=&feature._id;
  run;


proc sort data=cat_&feature.2;
   by exonchunk_id;
   run;
   
%mend;

%cat_feat(exon,19);
%cat_feat(transcript,36);
%cat_feat(gene,7);


data chunk_info;
   set conesa.exon_chunk2fusion;
   keep exonchunk_id chrom chunk_start chunk_end;
   run;

proc sort data=chunk_info nodup;
   by exonchunk_id;
run;

data exon_chunk_w_info_uniq;
  merge chunk_info (in=in1) cat_exon2 cat_transcript2 cat_gene2;
  by exonchunk_id;
  if exon_count=1 then flag_chunk_multiexon=0;
  else flag_chunk_multiexon=1;
    if transcript_count=1 then flag_chunk_multitranscript=0;
  else flag_chunk_multitranscript=1;
    if gene_count=1 then flag_chunk_multigene=0;
  else flag_chunk_multigene=1;
    if in1 then output;
  run;
  

/* Make permenant */


data conesa.exon_chunk_w_info;
   set exon_chunk_w_info;
   run;
   
data conesa.exon_chunk_w_info_uniq;
 set exon_chunk_w_info_uniq;
run;

