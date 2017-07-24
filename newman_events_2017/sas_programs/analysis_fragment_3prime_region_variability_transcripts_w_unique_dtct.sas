ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* 3' variability check: for transcripts that have detected unique pieces, how many have
   possible 3' variability? */

data xs_frag_flags;
   set event.xscripts_3prime_fragments;
run;

data xs_w_uniq;
   set event.xscripts_w_unique_by_bin;
   where num_unique_features_dtct > 0;
   keep transcript_id;
run;

proc sort data=xs_frag_flags;
   by transcript_id;
proc sort data=xs_w_uniq;
   by transcript_id;
run;

data uniq_xs_frag_flags;
   merge xs_w_uniq (in=in1) xs_frag_flags (in=in2);
   by transcript_id;
   if in1 and in2;
run;
*21108 of the 63807 xscripts with unique pieces retained;
*remainder either have their 3prime end  not detected at all;

proc freq data=uniq_xs_frag_flags;
   tables flag_possible_3prime_var
          flag_last_frag_in_xscript_on
          flag_possible_3prime_var*flag_last_frag_in_xscript_on  ;
run;


/*

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
           0 |      0 |  20487 |  20487
             |   0.00 |  97.06 |  97.06
             |   0.00 | 100.00 |
             |   0.00 |  98.42 |
    ---------+--------+--------+
           1 |    293 |    328 |    621
             |   1.39 |   1.55 |   2.94
             |  47.18 |  52.82 |
             | 100.00 |   1.58 |
    ---------+--------+--------+
    Total         293    20815    21108
                 1.39    98.61   100.00


*/
