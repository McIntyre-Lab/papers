/* iReckon GTF file does not give unique identifiers for iReckon transcripts (if RefSeq/Ensembl/etc are not available)
   I am importing the iReckon GTF file into SAS so I can append the gene ID to the transcript ID to ensure that each transcript
   has a unique ID. Then I will output the GTF file to process for extracting transcript sequences */


%macro formatGTF(sample);

/* Import iReckon GTF file */ 

proc import datafile="!MCLAB/event_analysis/alignment_output/ireckon/&sample./result.gtf"
            out=gtf dbms=tab replace; getnames=no; guessingrows=max;
run;

/* Get IDs for transcripts and genes -- I need to only really keep gene, transcript and exon_number */

data gtf2;
  set gtf;
  length gene_id $100.;
  length transcript_id $100.;
  format exon_number best32. ;
  if VAR3="transcript" then do;
    gene_id=compress(scan(scan(VAR9,1,";"),2,'"'));
    transcript_id=compress(scan(scan(VAR9,2,";"),2,'"'));
  end;
  else if VAR3="exon" then do;
    gene_id=compress(scan(scan(VAR9,1,";"),2,'"'));
    transcript_id=compress(scan(scan(VAR9,2,";"),2,'"'));
    exon_number=compress(scan(scan(VAR9,3,";"),2,'"'))+0;
  end;
  rename VAR1=chr VAR2=source VAR3=feature_type VAR4=start VAR5=stop VAR6=score
         VAR7=strand VAR8=frame;
  run;

/* Join gene and transcript ID */

data gtf3;
   set gtf2;
   format transcript_id2 $200.;
   format attributes $512.;
   transcript_id2=catx("_", gene_id, transcript_id);
   if feature_type="transcript" then attributes=catt('gene_id "',gene_id,'"; transcript_id "',transcript_id2,'";');
   else if feature_type="exon" then attributes=catt('gene_id "',gene_id,'"; transcript_id "',transcript_id2,
                                                    '"; exon_number "',exon_number,'";');
   keep chr source feature_type start stop score strand frame attributes;
run;

/* output reformatted GTF */

     data _null_;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     %let _EFIREC_ = 0;     /* clear export record count macro variable */
     file "!MCLAB/event_analysis/alignment_output/ireckon/&sample./result.formatted.gtf"
 delimiter='09'x DSD DROPOVER lrecl=32767;
    set  GTF3   end=EFIEOD;
        format chr $5. ;
        format source $11. ;
        format feature_type $10. ;
        format start best12. ;
        format stop best12. ;
        format score best12. ;
        format strand $1. ;
        format frame best12. ;
        format attributes char512. ;
      do;
        EFIOUT + 1;
        put chr $ @;
        put source $ @;
        put feature_type $ @;
        put start @;
        put stop @;
        put score @;
        put strand $ @;
        put frame @;
        put attributes char512. ;
        ;
      end;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     if EFIEOD then call symputx('_EFIREC_',EFIOUT);
     run;


%mend;

%formatGTF(NSC1);
%formatGTF(NSC2);

