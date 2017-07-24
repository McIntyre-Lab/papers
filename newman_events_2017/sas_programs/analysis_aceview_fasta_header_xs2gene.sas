ods listing; ods html close;
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';

/* Need to make an index for Aceview transcripts, their FASTA ID and their gene ID for RSEM */

proc import datafile="!PATCON/useful_human_data/aceview_hg19/transcripts/aceview_xscript_list.txt"
   out=xs_header dbms=tab replace; guessingrows=600000; getnames=no;
run;

/* Split FASTA headers to make the index */

data xs2gene_index;
   length gene_id $39.;
   length transcript_id $60.;
   length geneID_num $15.;
   set xs_header;
   gene_id=strip(scan(VAR1,3,"|"));
   transcript_id=strip(tranwrd(scan(VAR1,1,"|"), "MRNA:",""));
   geneID_num=strip(scan(VAR1,4,"|"));
   rename VAR1=transcript_fasta_id;
run;

/* Make permanant */

data event.hg19_aceview_xs2gene_fasta_index;
    retain gene_id transcript_fasta_id transcript_id geneID_num;
    set xs2gene_index;
run;

/* Export gene-to-xscript for all transcripts */

data xs2gene;
   set event.hg19_aceview_xs2gene_fasta_index;
   keep gene_id transcript_fasta_id;
run;

proc export data=xs2gene outfile="!PATCON/useful_human_data/aceview_hg19/transcripts/AceView.ncbi_37.all_mrnas_dna_gene2xs.txt" dbms=tab replace; putnames=no;
run;

