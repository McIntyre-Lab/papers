ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Import PacBio read counts from Lorena
I need to figure out how best to match these since my IDs are different!

Note: I looked at the estimates I got from running all PB transcripts in RSEM and Lorena's counts.
They are pretty close,  so for the purposes of this I will assume the IDs I used are okay for now
Update this when new lists arrive from Manu (I can reuse the same code!)
*/

%macro importLR(sample,name);
proc import datafile="!MCLAB/event_analysis/text_data/LR_SR_cor/&sample..5merge.abundance.txt"
     out=lr_count dbms=tab replace; guessingrows=20000;
run;

data &sample._lr;
  set lr_count;
  keep pbid count_fl count_nfl;
  rename pbid=pacbio_id count_fl=reads_FL_&name. count_nfl=reads_nonFL_&name. ;
run;

proc sort data=&sample._lr;
  by pacbio_id;
run;

%mend;

%importLR(NSC1,NPC1);
%importLR(NSC2,NPC2);
%importLR(OLD1,OPC1);
%importLR(OLD2,OPC2);


data event.pacbio_LR_read_counts;
   merge NSC1_LR NSC2_LR OLD1_LR OLD2_LR;
   by pacbio_id;
run;

