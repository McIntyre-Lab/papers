ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Import MISO gene annotations and relate to Refseq ID */

*import MISO SE-to-gene;

proc import datafile="/mnt/store/miso_sandbox/mm10/mm10/SE.mm10.gff3_to_ensGene.txt"
     out=se2gene dbms=tab replace;
     guessingrows=7936;
run;

*import Ensembl2Refseq;

    data WORK.ENS2REFSEQ    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile
'!MCLAB/useful_mouse_data/mm10/downloaded_files/mm10_ensembl2refseq_ucsc_browser.txt'
delimiter='09'x MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat _mm10_knownToRefSeq_name $10. ;
       informat mm10_knownToRefSeq_value $12. ;
       informat mm10_ensGene_name $3. ;
       informat mm10_ensGene_name2 $3. ;
       informat mm10_knownToEnsembl_name $10. ;
       informat mm10_knownToEnsembl_value $18. ;
       informat mm10_ensGtp_gene $18. ;
       informat mm10_ensGtp_transcript $18. ;
       informat mm10_ensemblToGeneName_name $18. ;
       informat mm10_ensemblToGeneName_value $14. ;
       informat mm10_knownGene_name $10. ;
       informat mm10_knownToLocusLink_name $10. ;
       informat mm10_knownToLocusLink_value $18. ;
       informat mm10_refGene_name $13. ;
       informat mm10_refGene_name2 $19. ;
       format _mm10_knownToRefSeq_name $10. ;
       format mm10_knownToRefSeq_value $12. ;
       format mm10_ensGene_name $3. ;
       format mm10_ensGene_name2 $3. ;
       format mm10_knownToEnsembl_name $10. ;
       format mm10_knownToEnsembl_value $18. ;
       format mm10_ensGtp_gene $18. ;
       format mm10_ensGtp_transcript $18. ;
       format mm10_ensemblToGeneName_name $18. ;
       format mm10_ensemblToGeneName_value $14. ;
       format mm10_knownGene_name $10. ;
       format mm10_knownToLocusLink_name $10. ;
       format mm10_knownToLocusLink_value $18. ;
       format mm10_refGene_name $13. ;
       format mm10_refGene_name2 $19. ;
    input
               _mm10_knownToRefSeq_name $
               mm10_knownToRefSeq_value $
               mm10_ensGene_name $
               mm10_ensGene_name2 $
               mm10_knownToEnsembl_name $
               mm10_knownToEnsembl_value $
               mm10_ensGtp_gene $
               mm10_ensGtp_transcript $
               mm10_ensemblToGeneName_name $
               mm10_ensemblToGeneName_value $
               mm10_knownGene_name $
               mm10_knownToLocusLink_name $
               mm10_knownToLocusLink_value $
               mm10_refGene_name $
               mm10_refGene_name2 $
   ;
   if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
   run;


data se2gene2;
   length gene_id2 $20.;
   set se2gene;
   do i=1 by 1 while(scan(gene_id,i,",") ^= "");
       gene_id2=scan(gene_id,i,",");
       output; end;
   keep event_id gene_id2;
   rename event_id=event_name gene_id2=ens_gene_id;
run;

data ens2refseq2;
  set ens2refseq;
  keep mm10_ensGtp_gene mm10_knownToLocusLink_value;
  rename mm10_ensGtp_gene=ens_gene_id mm10_knownToLocusLink_value=gene_id;
run;

proc sort data=se2gene2 nodup;
   by ens_gene_id event_name;
proc sort data=ens2refseq2 nodup;
   by ens_gene_id gene_id;
run;

/* Make permenant for now -- merge later */

data event.miso_se2gene;
  set se2gene2;
run;

data event.ensembl2refseq_gene_id;
  set ens2refseq2;
run;



