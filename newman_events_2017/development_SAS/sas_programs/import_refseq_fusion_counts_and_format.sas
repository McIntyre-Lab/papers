/* Import, format and export refseq counts */

proc import datafile='!MCLAB/conesa_pacbio/alignment_output/refseq_fusion_counts.csv'
            out=fusion_counts dbms=csv replace;
            guessingrows=550355;
run;

data fusion_counts_nsc fusion_counts_old;
   set fusion_counts;
   if sample_id = "NSC1" or sample_id="NSC2" then output fusion_counts_nsc;
   if sample_id = "OLD1" or sample_id="OLD1" then output fusion_counts_old;
   keep sample_id fusion_id apn;
run;


proc export data=fusion_counts_nsc
            outfile='!MCLAB/conesa_pacbio/analysis_output/refseq_fusion_counts_nsc.csv'
            dbms=csv replace;
run;


proc export data=fusion_counts_old
            outfile='!MCLAB/conesa_pacbio/analysis_output/refseq_fusion_counts_old.csv'
            dbms=csv replace;
run;
