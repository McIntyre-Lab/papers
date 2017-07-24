
/* Import fusion info */
proc import datafile='!MCLAB/conesa_pacbio/alignment_output/mm10_refseq_fragment_counts.csv'
            out=fragment_counts
            dbms=csv replace; guessingrows=1440945;
run;


data fragment_counts_nsc fragment_counts_old;
   set fragment_counts;
   if sample_id="NSC1" or sample_id="NSC2" then output fragment_counts_nsc;
   if sample_id="OLD1" or sample_id="OLD2" then output fragment_counts_old;
   keep sample_id fusion_id apn;
   rename fusion_id=fragment_id;
run;


proc export data=fragment_counts_nsc
            outfile='!MCLAB/conesa_pacbio/analysis_output/refseq_fragment_counts_nsc.csv'
            dbms=csv replace;
run;

proc export data=fragment_counts_old
            outfile='!MCLAB/conesa_pacbio/analysis_output/refseq_fragment_counts_old.csv'
            dbms=csv replace;
run;
