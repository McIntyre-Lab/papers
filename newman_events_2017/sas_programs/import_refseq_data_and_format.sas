/* Import, format and export refseq counts */

proc import datafile='!MCLAB/conesa_pacbio/alignment_output/refseq_fusion_counts.csv'
            out=fusion_counts dbms=csv replace;
            guessingrows=550355;
run;


     data WORK.SPLICING_COUNTS    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile '!MCLAB/conesa_pacbio/alignment_output/refseq_splicing_counts.csv' delimiter = ','
 MISSOVER DSD lrecl=32767 firstobs=2 ;
        informat sample_id $4. ;
        informat fusion_id $335. ;
        informat mapped_reads best32. ;
        informat read_length best32. ;
        informat region_length best32. ;
        informat region_depth best32. ;
        informat reads_in_region best32. ;
        informat apn best32. ;
        informat rpkm best32. ;
        informat mean best32. ;
        informat std best32. ;
        informat cv best32. ;
        format sample_id $4. ;
        format fusion_id $335. ;
        format mapped_reads best12. ;
        format read_length best12. ;
        format region_length best12. ;
        format region_depth best12. ;
        format reads_in_region best12. ;
        format apn best12. ;
        format rpkm best12. ;
        format mean best12. ;
        format std best12. ;
        format cv best12. ;
     input
                 sample_id $
                 fusion_id $
                 mapped_reads
                 read_length
                 region_length
                 region_depth
                 reads_in_region
                 apn
                 rpkm
                 mean
                 std
                 cv
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;

data fusion_counts2;
   set fusion_counts;
   keep sample_id fusion_id apn;
run;

data splicing_counts2;
   set splicing_counts;
   keep sample_id fusion_id apn;
   rename fusion_id=event_id;
run;

proc export data=fusion_counts2
            outfile='!MCLAB/conesa_pacbio/analysis_output/refseq_fusion_counts.csv'
            dbms=csv replace;
run;

proc export data=splicing_counts2
            outfile='!MCLAB/conesa_pacbio/analysis_output/refseq_splicing_counts.csv'
            dbms=csv replace;
run;

/* Import and format fusion info */



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

