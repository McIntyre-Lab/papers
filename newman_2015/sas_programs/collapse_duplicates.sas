/***** Collapse duplicate events *******/

libname splice '/media/jrbnewman/ac89a883-cbf2-4ed0-8e2f-19c25fead575/splice';
/* first count the number of "tiny" events - these are events that are smaller than the read size */
/* For the diabetes data this is 50bp */

data splicing_event_size;
   set splice.splicing_events_w_exon_info;
   event_size=donor_size+acceptor_size;
   if event_size lt 50 then flag_tiny_event=1;
   else flag_tiny_event=0;
   keep event_id event_size flag_tiny_event donor_size acceptor_size;
run;

proc sort data=splicing_event_size;
   by event_size donor_size acceptor_size;
run;

proc freq data=splicing_event_size noprint;
    tables flag_tiny_event / out=tiny_event_count;
run;

/* 12740030 events ge 50bp, 193971 events lt 50bp */

 data small_events;
   set splicing_event_size;
   if event_size lt 50 then output;
run;

proc export data=small_events
   outfile='/home/jrbnewman/McLab/junction_annotations/pipeline_output/aceview_hg19/small_exons.csv'
   dbms=csv
   replace;
run;

data splicing_events;
  set splice.splicing_events_w_exon_info;
  event_coords=cats(gene_id,"|",chr,":",donor_start,"-",donor_stop,":",acceptor_start,"-",acceptor_stop,":",strand);
  event_size=donor_size+acceptor_size;
  keep event_id donor_exon acceptor_exon event_coords transcript_id  event_size num_transcripts;
run;

proc sort data=splicing_events;
  by event_coords event_size event_id;
run;

data first_event;
   set splicing_events;
   by event_coords;
   if first.event_coords then output;
   drop donor_exon acceptor_exon;
   rename event_id = new_event_id;
run;

/* use new_event_id to cat event_id, donor_exon, acceptor_exon */
proc sort data=first_event;
  by event_coords;
run;

proc sort data=splicing_events;
  by event_coords;
run;

data splicing_events_W_new_id oops1 oops2;
   merge first_event (in=in1) splicing_events (in=in2);
   by event_coords;
   if in1 and in2 then output splicing_events_w_new_id;
   else if in1 then output oops1;
   else output oops2;
run;

/* Add short exon info back in */
proc sort data=splicing_events_W_new_id;
  by new_event_id;
run;

proc sort data=splicing_event_size(rename=(event_id=new_event_id));
   by new_event_id;
run;

data collapsed_events_w_size event_not_ref no_size_oops;
   merge splicing_event_size (in=in1) splicing_events_W_new_id (in=in2);
   by new_event_id;
   if in1 and in2 then output collapsed_events_w_size;
   else if in1 then output event_not_ref;
   else output no_size_oops;
run;

data collapsed_events_for_cat;
   set collapsed_events_w_size;
   length donor_acceptor_num $13.;
   length donor_num $6.;
   length acceptor_num $6.;
   length gene_id $36.;
   if donor_exon='intron' then donor_num='intron';
   else donor_num=scan(donor_exon,2,':');
   if acceptor_exon='intron' then acceptor_num='intron';
   else acceptor_num=scan(acceptor_exon,2,':');
   donor_acceptor_num=catx('_', donor_num, acceptor_num);
   gene_id=scan(donor_exon,1,':');
   keep new_event_id donor_acceptor_num gene_id transcript_id donor_num acceptor_num;
run;

/* Clean up here */
proc datasets noprint;
delete splicing_event_size tiny_event_count small_events splicing_events first_event splicing_events_W_new_id oops1 oops2
collapsed_events_w_size event_not_ref no_size_oops;
run; quit;

/* Cat event_ids, donor_exons, acceptor_exons, transcript_ids; sum num_transcripts */
/*** cat event_ids ***/
data uniq_events;
   set collapsed_events_for_cat;
   keep new_event_id donor_acceptor_num;
   run;

proc sort data=uniq_events nodup;
   by new_event_id donor_acceptor_num;
run;

/* get counts first */
proc freq noprint data=uniq_events;
   tables new_event_id / out=event_count;
run;

proc sort data=event_count;
  by descending count;
run;
*max=377 events per ref_event;

data cat_event; 
  array event[377] $ 13.;
  retain event1-event377;
  set uniq_events;
  by new_event_id;
  if first.new_event_id then do;
     call missing(of event1-event377);
     records = 0;
  end;
  records + 1;
  event[records]=donor_acceptor_num;
  if last.new_event_id then output;
run;

  *clean up the output file;
data cat_event2;
  set cat_event;
  length cat_event_id $ 2500.;
  rename records= num_events;
         cat_event_id= catx("|", OF event1-event377);
  drop event1-event377 donor_acceptor_num;
  run;

/* Cat together donor exons */
data uniq_donor_exons;
   set collapsed_events_for_cat;
   keep new_event_id donor_num;
   run;

/* Cat together donors by event */
proc sort data=uniq_donor_exons nodup; *drop duplicated donor exons!;
    by new_event_id donor_num;
run;

/* get counts first */
proc freq noprint data=uniq_donor_exons;
   tables new_event_id / out=donor_count;
run;

proc sort data=donor_count;
  by descending count;
run;
*max=24 donor exons per event;

data cat_donors; 
  array donors[24] $ 6;
  retain donors1-donors24;
  set uniq_donor_exons;
  by new_event_id;
  if first.new_event_id then do;
     call missing(of donors1-donors24);
     records = 0;
  end;
  records + 1;
  donors[records]=donor_num;
  if last.new_event_id then output;
run;

  *clean up the output file;
data cat_donors2;
  set cat_donors;
  length donor_exons $ 170;
  rename records= num_donor_exons;
         donor_exons= catx("|", OF donors1-donors24);
  drop donors1-donors24 donor_num;
  run;

/* Cat together acceptor exons */
data uniq_acceptor_exons;
   set collapsed_events_for_cat;
   keep new_event_id acceptor_num;
   run;

/* Cat together acceptors by event */
proc sort data=uniq_acceptor_exons nodup; *drop duplicated acceptor exons!;
    by new_event_id acceptor_num;
run;

/* get counts first */
proc freq noprint data=uniq_acceptor_exons;
   tables new_event_id / out=acceptor_count;
run;

proc sort data=acceptor_count;
  by descending count;
run;
*max=29 acceptor exons per event;


data cat_acceptors; 
  array acceptors[29] $ 6;
  retain acceptors1-acceptors29;
  set uniq_acceptor_exons;
  by new_event_id;
  if first.new_event_id then do;
     call missing(of acceptors1-acceptors29);
     records = 0;
  end;
  records + 1;
  acceptors[records]=acceptor_num;
  if last.new_event_id then output;
run;

  *clean up the output file;
data cat_acceptors2;
  set cat_acceptors;
  length acceptor_exons $ 200;
  rename records= num_acceptor_exons;
         acceptor_exons= catx("|", OF acceptors1-acceptors29);
  drop acceptors1-acceptors29 acceptor_num;
  run;

/* Cat together transcripts */
data uniq_xscripts;
   set collapsed_events_for_cat;
   keep new_event_id transcript_id;
   run;

/* Cat together acceptors by event */
proc sort data=uniq_xscripts nodup; *drop duplicated transcripts!;
    by new_event_id transcript_id;
run;

/* get counts first */
proc freq noprint data=uniq_xscripts;
   tables new_event_id / out=xscripts_count;
run;

proc sort data=xscripts_count;
  by descending count;
run;
*max=29 transcripts per event (actually more, but this is 29 "concatenated transcripts" per event);

data cat_xscripts; 
  array xscripts[29] $ 39;
  retain xscripts1-xscripts29;
  set uniq_xscripts;
  by new_event_id;
  if first.new_event_id then do;
     call missing(of xscripts1-xscripts29);
     records = 0;
  end;
  records + 1;
  xscripts[records]=transcript_id;
  if last.new_event_id then output;
run;

  *clean up the output file;
data cat_xscripts2;
  set cat_xscripts;
  length transcripts $ 763;
         transcripts= catx("|", OF xscripts1-xscripts29);
  drop xscripts1-xscripts29 transcript_id records;
  rename transcripts=transcript_id;
  run;

/* Count number of unique transcripts per collapsed event */

data uniq_xscript_counts;
   set uniq_xscripts;
   if transcript_id='' then xscript_count=0;
   else xscript_count=count(transcript_id,'|')+1;
   keep new_event_id xscript_count;
run;


/* Sum counts */
proc sort data=uniq_xscript_counts;
   by new_event_id;
run;

proc means data=uniq_xscript_counts noprint;
   by new_event_id;
   var xscript_count;
   output out=num_xscripts sum=total_transcripts;
run;

/* Merge together */
proc sort data=cat_event2;
by new_event_id;
run;

proc sort data=cat_xscripts2;
by new_event_id;
run;

proc sort data=cat_donors2;
by new_event_id;
run;

proc sort data=cat_acceptors2;
by new_event_id;
run;

proc sort data=num_xscripts;
by new_event_id;
run;

data collapsed_event_info;
   merge cat_event2 num_xscripts cat_xscripts2 cat_donors2 cat_acceptors2;
   by new_event_id;
   drop _TYPE_ _FREQ_;
   rename new_event_id=event_id;
run;

/* Collapse full event information */
data splice_events_all;
   set splice.splicing_events_w_exon_info;
   drop donor_exon acceptor_exon num_transcripts transcript_id;
run;

proc sort data=collapsed_event_info;
   by event_id;
run;

proc sort data=splice_events_all;
   by event_id;
run;

data collapsed_events_final not_ref no_event_oops;
   merge collapsed_event_info (in=in1) splice_events_all (in=in2);
   by event_id;
   if in1 and in2 then output collapsed_events_final;
   else if in1 then output no_event_oops;
   else output not_ref;
run;

/* Make permenant */
data splice.collapsed_splicing_events;
    length new_event_id $2550.;
    length feature1_id $250.;
    length feature2_id $250.;
    set collapsed_events_final;
    new_event_id=catx(':', gene_id, cat_event_id);
    feature1_id=catx(':', gene_id, donor_exons);
    feature2_id=catx(':', gene_id, acceptor_exons);
    rename total_transcripts=num_transcripts;
    event_size=donor_size+acceptor_size;
    if event_size lt 50 then flag_tiny_event=1;
    else flag_tiny_event=0;
    drop event_id cat_event_id;
run;

/* Delete temp datasets */
proc datasets noprint;
delete splicing_event_size tiny_event_count small_events
splicing_events first_event splicing_events_W_new_id oops1 oops2
collapsed_events_w_size event_not_ref no_size_oops
collapsed_events_for_cat uniq_events event_count
cat_event cat_event2 uniq_donor_exons
donor_count cat_donors cat_donors2 uniq_acceptor_exons
acceptor_count cat_acceptors cat_acceptors2 uniq_xscripts
xscripts_count cat_xscripts cat_xscripts2 uniq_xscript_counts
num_xscripts collapsed_event_info splice_events_all
collapsed_events_final not_ref no_event_oops;
run; quit;

