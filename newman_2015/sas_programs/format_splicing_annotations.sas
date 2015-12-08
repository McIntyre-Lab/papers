/***** Format splicing event annotations *******/

libname splice '/media/jrbnewman/ac89a883-cbf2-4ed0-8e2f-19c25fead575/splice';

data splicing_events_formatting;
   set splice.collapsed_splicing_events;
   length feature1_type $13.;
   length feature2_type $13.;
   rename new_event_id=event_id;
   if strand='+' then do;
       if event_type='exon_junction' then do;
            feature1_type='exon_donor';
            feature2_type='exon_acceptor';
            end;
       else if event_type='intron_retention' then do;
            feature1_type='exon_donor';
            feature2_type='intron';
            end;
       end;
   else if strand='-' then do; 
       if event_type='exon_junction' then do;
            feature1_type='exon_acceptor';
            feature2_type='exon_donor';
            end;
       else if event_type='intron_retention' then do;
            feature1_type='intron';
            feature2_type='exon_donor';
            end;
       end;
run;

/* Rearrange columns */


data splicing_events_reordered;
   retain event_id event_type num_events gene_id chr strand event_size flag_tiny_event num_transcripts transcript_id feature1_id feature1_type
feature2_id feature2_type donor_start donor_stop num_donor_exons donor_exons donor_size flag_short_donor acceptor_start acceptor_stop num_acceptor_exons
acceptor_exons acceptor_size flag_short_acceptor flag_junction_annotated flag_intron_retention num_skipped_exons flag_exonskip flag_alt_donor
flag_alt_acceptor;
   set splicing_events_formatting;
rename flag_tiny_event=flag_event_short donor_start=feature1_start donor_stop=feature1_stop num_donor_exons=num_feature1 donor_exons=feature1_list donor_size=feature1_size flag_short_donor=flag_feature1_short acceptor_start=feature2_start acceptor_stop=feature2_stop num_acceptor_exons=num_feature2 acceptor_exons=feature2_list acceptor_size=feature2_size flag_short_acceptor=flag_feature2_short;
drop acceptor_group donor_group;
run;


/* Make permenant */

data splice.splicing_events_annotations;
    set splicing_events_reordered;
run;


/* Delete temp datasets */

proc datasets noprint;
  delete splicing_events_formatting splicing_events_reordered;
run;quit;

/* Export catalogue CSV */

proc export data=splice.splicing_events_annotations
    outfile='/home/jrbnewman/McLab/junction_annotations/pipeline_output/aceview_hg19/hg19_splicing_catalogue_74bp.csv'
    dbms=csv
    replace;
run;

