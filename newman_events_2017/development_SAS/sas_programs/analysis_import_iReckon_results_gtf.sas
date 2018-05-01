ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/splicing/sas_data';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';

/* Import iReckon results and process 

   iReckon output is a GTF file consisting of transcripts and exons
   Each transcript has the following attributes:
         gene_id	iReckon assigns this, so the actual "geneID" is not useful
         transcript_id	Gleaned from annotation, else assigned by iReckon
         RPKM		RPKM for transcript
         frac		Not sure? Fraction of gene expression?
         conf_lo	?
         conf_hi	?
         frac		?
         cov		?	

   Each exon has the following additional attribute:
         exon_number	number of exon in isoform
*/

%macro importGTF(sample);
proc import datafile="/mnt/store/event_sandbox/ireckon/&sample./result.gtf"
     out=&sample._ireckon dbms=tab replace;
     getnames=no; guessingrows=170000;
run;

/* Parse results GTF */

data &sample._ireckon2;
  set &sample._ireckon;
  length gene_id $100.;
  length transcript_id $100.;
  format exon_number best32. ;
  format rpkm best32. ;
  format frac1 best32. ;
  format conf_lo best32. ;
  format conf_hi best32. ;
  format frac2 best32. ;
  format cov best32. ;
  if VAR3="transcript" then do;
    gene_id=compress(scan(scan(VAR9,1,";"),2,'"'));
    transcript_id=compress(scan(scan(VAR9,2,";"),2,'"'));
    rpkm=compress(scan(scan(VAR9,3,";"),2,'"'))+0;
    frac1=compress(scan(scan(VAR9,4,";"),2,'"'))+0;
    conf_lo=compress(scan(scan(VAR9,5,";"),2,'"'))+0;
    conf_hi=compress(scan(scan(VAR9,6,";"),2,'"'))+0;
    frac2=compress(scan(scan(VAR9,7,";"),2,'"'))+0;
    cov=compress(scan(scan(VAR9,8,";"),2,'"'))+0;
  end;

  else if VAR3="exon" then do;
    gene_id=compress(scan(scan(VAR9,1,";"),2,'"'));
    transcript_id=compress(scan(scan(VAR9,2,";"),2,'"'));
    exon_number=compress(scan(scan(VAR9,3,";"),2,'"'))+0;
    rpkm=compress(scan(scan(VAR9,4,";"),2,'"'))+0;
    frac1=compress(scan(scan(VAR9,5,";"),2,'"'))+0;
    conf_lo=compress(scan(scan(VAR9,6,";"),2,'"'))+0;
    conf_hi=compress(scan(scan(VAR9,7,";"),2,'"'))+0;
    frac2=compress(scan(scan(VAR9,8,";"),2,'"'))+0;
    cov=compress(scan(scan(VAR9,9,";"),2,'"'))+0;
  end;
  else delete;
  drop VAR9;
  rename VAR1=chr VAR2=source VAR3=feature_type VAR4=start VAR5=stop VAR6=score
         VAR7=strand VAR8=frame;
  run;

data event.iReckon_&sample._output;
  set &sample._ireckon2;
run;

%mend;

%importGTF(NSC1);
%importGTF(NSC2);

