libname fru '!HOME/mclab/Fru_network/sasdata';

proc import out=work.all_coverage_spikes
            datafile="!HOME/mclab/Fru_network/pipeline_output/coverage_spikes/all_coverage_spikes.csv"
            dbms=csv replace;
            getnames=yes;
            datarow=2;
            run;

 data WORK.ALL_COVERAGE_SPIKES    ;
       %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
       infile '!HOME/mclab/Fru_network/pipeline_output/coverage_spikes/all_coverage_spikes.csv' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
          informat spike_id $4. ;
          informat sample_id $20. ;
          informat mapped_reads best32. ;
          informat reads_in_exon best32. ;
          informat coverage_in_exon best32. ;
          informat exon_length best32. ;
          informat apn best32. ;
          informat rpkm best32. ;
          format spike_id $4. ;
          format sample_id $20. ;
          format mapped_reads best12. ;
          format reads_in_exon best12. ;
          format coverage_in_exon best12. ;
          format exon_length best12. ;
          format apn best12. ;
          format rpkm best12. ;
       input
                   spike_id $
                   sample_id
                   mapped_reads
                   reads_in_exon
                   coverage_in_exon
                   exon_length
                   apn
                   rpkm
       ;
       if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
       run;

data fru.all_coverage_spikes;
    set all_coverage_spikes;
    logRPKM = log(RPKM+100);
    run;
