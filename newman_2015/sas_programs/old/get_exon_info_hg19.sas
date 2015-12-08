/* set libraries */

libname splice '/home/jrbnewman/McLab/junction_annotations/sas_data/';
libname splice2 '/media/jrbnewman/SAS_WRK1/';

/* import exons */

filename exons '/home/jrbnewman/McLab/junction_annotations/generated_files/hg19_exon_list.csv';


/*proc import datafile=exons
    out=hg19_exon_list
    dbms=csv
    replace;
    getnames=yes;
    guessingrows=1400000;
run;*/

   data WORK.HG19_EXON_LIST    ;
   %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
   infile EXONS delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
      informat chr $4. ;
      informat start best32. ;
      informat stop best32. ;
      informat strand $1. ;
      informat exon_name $49. ;
      informat exon_id $61. ;
      informat xscripts_cat $53. ;
      informat gene_cat $36. ;
      format chr $4. ;
      format start best12. ;
      format stop best12. ;
      format strand $1. ;
      format exon_name $49. ;
      format exon_id $61. ;
      format xscripts_cat $53. ;
      format gene_cat $36. ;
   input
               chr $
               start
               stop
               strand $
               exon_name $
               exon_id $
               xscripts_cat $
               gene_cat $
   ;
   if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
   run;


/* drop superfluous variables for now and drop dups */


data hg19_exon_list2;
    set hg19_exon_list;
    drop exon_name xscripts_cat;
run;


proc sort data=hg19_exon_list2 nodups;
    by chr start stop strand exon_id gene_cat;
run;

*dropped 721496 duplicated exons;


/* Group exons by a reference exon */

/* Sort by gene, chr, start, stop */

proc sort data=hg19_exon_list2;
   by gene_cat chr start stop;
run;

/* add in reference exon values - this will help with defining alt donors and acceptors */
/* if exons overlap then take the first exon as a reference exon */

data exon_references;
   set hg19_exon_list2;
   by gene_cat;
   length ref_exon $50.;
   retain ref_exon_stop;
   retain ref_exon;
   retain ref_exon_start;
   retain ref_exon_strand;
   retain exon_number;
   if first.gene_cat then do; *first exon is by default a reference exon;
      ref_exon=exon_id;
      ref_exon_start=start;
      ref_exon_stop=stop;
      ref_exon_strand=strand;
      exon_number=1;
      end;
   else do; *exon is not the first of gene, but does it overlap with previous reference?;
      if start lt ref_exon_stop then do; *exon overlaps with previous reference exon;
         ref_exon=ref_exon;
         ref_exon_start=ref_exon_start;
         ref_exon_stop=ref_exon_stop;
         ref_exon_strand=ref_exon_strand;
         exon_number=exon_number;
         end;
      else do; *exon does not overlap with previous reference, therefore new exon;
         ref_exon=exon_id;
         ref_exon_start=start;
         ref_exon_stop=stop;
         ref_exon_strand=strand;
         exon_number=exon_number+1;
         end;
      end;
run;

/* Flag first and alternative first exons */

data flag_first_exon;
    set exon_references;
    by gene_cat;
    if first.gene_cat then do;
        flag_firstexon=1;
        flag_altfirstexon=0;
        end;
    else if exon_number=1 then do;
        flag_firstexon=0;
        flag_altfirstexon=1;
        end;
    else do;
        flag_firstexon=0;
        flag_altfirstexon=0;
        end;
run;

/* Flag last and alternative last exons */

proc sort data=flag_first_exon;
    by gene_cat chr descending stop descending start;
run;

data flag_last_exon;
    set flag_first_exon;
    retain max_exon_num;
    final_exon_num=max_exon_num;
    by gene_cat;
    if first.gene_cat then do;
        flag_lastexon=1;
        flag_altlastexon=0;
        max_exon_num=exon_number;
        end;
    else if exon_number=final_exon_num then do;
        flag_lastexon=0;
        flag_altlastexon=1;
        end;
    else do;
        flag_lastexon=0;
        flag_altlastexon=0;
        end;
run;



/* Make list of all exons with donor info */
data donor_exons;
   set flag_last_exon;
   keep exon_id start stop ref_exon ref_exon_start ref_exon_stop flag_firstexon flag_altfirstexon gene_cat;
   rename exon_id=exonA;
   rename start=exonA_start;
   rename stop=exonA_stop;
   rename ref_exon=exonA_ref;
   rename ref_exon_start=exonA_ref_start;
   rename ref_exon_stop=exonA_ref_stop;
   rename gene_cat=geneA_id;
run;

/* Make list of all exons with acceptor info */
data acceptor_exons;
   set flag_last_exon;
   keep exon_id start stop ref_exon ref_exon_start ref_exon_stop flag_lastexon flag_altlastexon gene_cat; 
   rename exon_id=exonB;
   rename start=exonB_start;
   rename stop=exonB_stop;
   rename ref_exon=exonB_ref;
   rename ref_exon_start=exonB_ref_start;
   rename ref_exon_stop=exonB_ref_stop;
   rename gene_cat=geneB_id;
run;

/* Make datasets permenant */
data splice.donor_exons_hg19;
   set donor_exons;
run;

data splice.acceptor_exons_hg19;
   set acceptor_exons;
run;

data splice2.exon_info_hg19;
    set exon_references;
run;

/* Clean up */
    proc datasets nolist;
        delete exon_references
        hg19_exon_list flag_first_exon flag_last_exon
        hg19_exon_list2 donor_exons acceptor_exons
        ;
        run;
        quit;











