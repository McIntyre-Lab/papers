/home/jrbnewman/McLab/ethanol/mel_pilot/analysis_output/dmel_fb611_fusion_info.csv


gene_id
gene_name
fusion_id
exon_id
flag_multigene
chrom
start
end

libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname splicing '!MCLAB/conesa_pacbio/sas_data/splicing';

/* Format and export fusion info */

/* need columns IN THIS ORDER:
   gene_id gene_name fusion_id exon_id flag_multigene chrom start end */

data fusion_info;
   retain gene_id gene_name fusion_id exon_id flag_multigene chrom fusion_start fusion_end;
   length exon_id $12.;
   length gene_id $15.;
   set conesa.fusions_si_info;
   if chrom='chr1' and gene_id="PB.1450" then do;
      gene_id=tranwrd(gene_id,"PB.1450","PB.1450_chr1");
      exon_id=tranwrd(exon_id,"PB.1450","PB.1450_chr1");
   end;
   if chrom='chrM' and gene_id="PB.33" then do;
      gene_id=tranwrd(gene_id,"PB.33","PB.33_chrM");
      exon_id=tranwrd(exon_id,"PB.33","PB.33_chrM");
   end;
   if chrom='chr2' and gene_id="PB.3783" then do;
      gene_id=tranwrd(gene_id,"PB.3783","PB.3783_chr2");
      exon_id=tranwrd(exon_id,"PB.3783","PB.3783_chr2");
   end;
   if chrom='chr13' and gene_id="PB.4412" then do;
      gene_id=tranwrd(gene_id,"PB.4412","PB.4412_chr13");
      exon_id=tranwrd(exon_id,"PB.4412","PB.4412_chr13");
   end;
   if chrom='chr1' and gene_id="PB.4642" then do;
      gene_id=tranwrd(gene_id,"PB.4642","PB.4642_chr1");
      exon_id=tranwrd(exon_id,"PB.4642","PB.4642_chr1");
   end;
   if chrom='chr18' and gene_id="PB.5261" then do;
      gene_id=tranwrd(gene_id,"PB.5261","PB.5261_chr18");
      exon_id=tranwrd(exon_id,"PB.5261","PB.5261_chr18");
   end;
   gene_name=gene_id;
   keep gene_id gene_name fusion_id exon_id flag_multigene chrom fusion_start fusion_end;
   rename fusion_start=start fusion_end=end;
run;


data check;
   set fusion_info;
   where gene_id ? "PB.5261";
run;


/* Update splicing catalogue */

data splicing_info;
   set splicing.splicing_events_annotations;
   length gene_id2 $15.;
   format gene_id2 $15.;
   informat gene_id2 $15.;
   length transcript_id2 $460.;
   format transcript_id2 $460.;
   informat transcript_id2 $460.;
   length event_id2 $2565.;
   format event_id2 $2565.;
   informat event_id2 $2565.;
   length feature1_id2 $260.;
   length feature2_id2 $260.;
   format feature1_id2 $260.;
   format feature2_id2 $260.;
   informat feature1_id2 $260.;
   informat feature2_id2 $260.;

   if chr='chr1' and gene_id="PB.1450" then do;
      gene_id2=tranwrd(gene_id,"PB.1450","PB.1450_chr1");
      transcript_id2=tranwrd(transcript_id,"PB.1450","PB.1450_chr1");
      event_id2=tranwrd(event_id,"PB.1450","PB.1450_chr1");
      feature1_id2=tranwrd(feature1_id,"PB.1450","PB.1450_chr1");
      feature2_id2=tranwrd(feature2_id,"PB.1450","PB.1450_chr1");
   end;
   else if chr='chrM' and gene_id="PB.33" then do;
      gene_id2=tranwrd(gene_id,"PB.33","PB.33_chrM");
      transcript_id2=tranwrd(transcript_id,"PB.33","PB.33_chrM");
      event_id2=tranwrd(event_id,"PB.3783","PB.3783_chr2");
      feature1_id2=tranwrd(feature1_id,"PB.3783","PB.3783_chr2");
      feature2_id2=tranwrd(feature2_id,"PB.3783","PB.3783_chr2");
   end;
   else if chr='chr2' and gene_id="PB.3783" then do;
      gene_id2=tranwrd(gene_id,"PB.3783","PB.3783_chr2");
      transcript_id2=tranwrd(transcript_id,"PB.3783","PB.3783_chr2");
      event_id2=tranwrd(event_id,"PB.3783","PB.3783_chr2");
      feature1_id2=tranwrd(feature1_id,"PB.3783","PB.3783_chr2");
      feature2_id2=tranwrd(feature2_id,"PB.3783","PB.3783_chr2");
   end;
   else if chr='chr13' and gene_id="PB.4412" then do;
      gene_id2=tranwrd(gene_id,"PB.4412","PB.4412_chr13");
      transcript_id2=tranwrd(transcript_id,"PB.4412","PB.4412_chr13");
      event_id2=tranwrd(event_id,"PB.4412","PB.4412_chr13");
      feature1_id2=tranwrd(feature1_id,"PB.4412","PB.4412_chr13");
      feature2_id2=tranwrd(feature2_id,"PB.4412","PB.4412_chr13");
   end;
   else if chr='chr1' and gene_id="PB.4642" then do;
      gene_id2=tranwrd(gene_id,"PB.4642","PB.4642_chr1");
      transcript_id2=tranwrd(transcript_id,"PB.4642","PB.4642_chr1");
      event_id2=tranwrd(event_id,"PB.4642","PB.4642_chr1");
      feature1_id2=tranwrd(feature1_id,"PB.4642","PB.4642_chr1");
      feature2_id2=tranwrd(feature2_id,"PB.4642","PB.4642_chr1");
   end;
   else if chr='chr18' and gene_id="PB.5261" then do;
      gene_id2=tranwrd(gene_id,"PB.5261","PB.5261_chr18");
      transcript_id2=tranwrd(transcript_id,"PB.5261","PB.5261_chr18");
      event_id2=tranwrd(event_id,"PB.5261","PB.5261_chr18");
      feature1_id2=tranwrd(feature1_id,"PB.5261","PB.5261_chr18");
      feature2_id2=tranwrd(feature2_id,"PB.5261","PB.5261_chr18");
   end;
   else do;
      gene_id2=gene_id;
      transcript_id2=transcript_id;
      event_id2=event_id;
      feature1_id2=feature1_id;
      feature2_id2=feature2_id; end;
   drop gene_id transcript_id event_id feature1_id feature2_id;
   rename gene_id2=gene_id transcript_id2=transcript_id
      event_id2=event_id feature1_id2=feature1_id
      feature2_id2=feature2_id;
run;

/* Make splicing info permenant, export fusion info */

data conesa.splicing_annotations_fixed;
   retain event_id event_type num_events gene_id chr strand event_size flag_event_short
   num_transcripts transcript_id feature1_id feature1_type feature2_id feature2_type;
   set splicing_info;
run;

proc export data=fusion_info outfile="!MCLAB/conesa_pacbio/analysis_output/conesa_pacbio_mm10_fusion_info.csv" dbms=csv replace;
run;


