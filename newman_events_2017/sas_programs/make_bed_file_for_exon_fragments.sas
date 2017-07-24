ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/*
Create and export a BED file of exon fragments, so that I can extract fragment sequences and BLAST them against
RefSeq and PacBio transcripts

Can use to flag ambiguous sequences. I am only interested in fragments above a length of 12
*/

data fragments;
  retain chr fragment_start fragment_end fragment_id;
  set mm10.mm10_exon_fragment_flagged;
  keep fragment_id fragment_start fragment_end chr;
run;

/* Export BED */

proc export data=fragments outfile="!MCLAB/event_analysis/references/mm10_refseq_exon_fragments.bed"
     dbms=tab replace; putnames=no;
run;

