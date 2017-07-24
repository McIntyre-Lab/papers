

/* Import fusion info */
proc import datafile='!MCLAB/conesa_pacbio/created_files/splicing_total/conesa_refseq_fusions_si.bed'
            out=fusions_bed
            dbms=tab replace; guessingrows=275177; getnames=no;
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


proc import datafile='!MCLAB/conesa_pacbio/created_files/splicing_total/conesa_refseq_fusions_si.tsv'
            out=fusions_info
            dbms=tab replace; guessingrows=355270;
run;

proc sort data=fusions_info;
   by fusion_id;
proc sort data=fusions_bed2;
   by fusion_id;
run;

data fusions_info_w_coord;
   merge fusions_info (in=in1) fusions_bed2 (in=in2);
   by fusion_id;
   if in1 and in2;
   gene_id=scan(exon_id,1,":");
run;


/* need columns IN THIS ORDER:
   gene_id gene_name fusion_id exon_id flag_multigene chrom start end */

data fusion_info2;
   retain gene_id gene_name fusion_id exon_id flag_multigene chrom fusion_start fusion_end;
   set fusions_info_w_coord;
   gene_name=gene_id;
   keep gene_id gene_name fusion_id exon_id flag_multigene chrom fusion_start fusion_end;
   rename fusion_start=start fusion_end=end;
run;



proc export data=fusion_info2 outfile="!MCLAB/conesa_pacbio/analysis_output/conesa_refseq_mm10_fusion_info.csv" dbms=csv replace;
run;

