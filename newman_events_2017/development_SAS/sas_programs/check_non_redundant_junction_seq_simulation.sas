/* Create a non-redundant set of junction sequences (maximum coordinates) */

data juncs;
  set eventloc.splicing_events_annot_simulation;
  junction_id=catx(":",chr,feature1_stop,feature2_start,strand);
  keep event_id gene_id chr strand feature1_start feature1_stop feature2_start feature2_stop
  event_size junction_id event_type;
run;

proc sort data=juncs;
   by junction_id;
proc means data=juncs noprint;
   by junction_id;
   var feature1_start feature1_stop feature2_start feature2_stop;
   output out=max_coord_junc min(feature1_start)=donor_start_min
                             max(feature1_stop)=donor_stop
			     min(feature2_start)=acceptor_start
			     max(feature2_stop)=acceptor_stop_max;
run;

proc sort data=juncs;
   by junction_id;
proc sort data=max_coord_junc;
   by junction_id;
run;

data juncs_w_max;
  merge juncs (in=in1) max_coord_junc (in=in2);
  by junction_id;
  if in1 and in2;
run;

/* Make perm */

data eventloc.simul_junctions_distinct_coord;
   set juncs_w_max;
run;

data juncs_w_max2;
  set juncs_w_max;
  event_size=(donor_stop-donor_start_min)+(acceptor_stop_max-acceptor_start);
  keep junction_id event_size donor_start_min donor_stop acceptor_start acceptor_stop_max chr strand event_type;
run;

proc sort data=juncs_w_max2 nodup;
  by junction_id;
run;

data junc border;
  set juncs_w_max2;
  if event_type="exon_junction" then output junc;
  else if event_type="intron_retention" then output border;
  drop event_type;
run;

proc sort data=junc;
   by junction_id;
proc sort data=border;
   by junction_id;
run;

data juncs_w_max3;
  merge junc (in=in1) border (in=in2);
  by junction_id;
  length event_type $20.;
  if in1 and in2 then event_type="exon_junction";
  else if in1 then event_type="exon_junction";
  else event_type="intron_retention";
run;

data info_for_bed;
   length score $1.;
   length color $7.;
   set juncs_w_max3;
   score='.';
   color='255,0,0';
   block1_start=0;

   if event_type='intron_retention' then do;
      num_blocks=1;
      if strand='+' then do;
         block1_length=event_size+1;
         totalstart1=donor_start_min;
         totalstop1=acceptor_stop_max+1;
         totalstart2=donor_start_min;
         totalstop2=acceptor_stop_max+1;
         end;
      if strand='-' then do;
         block1_length=event_size+1;
         totalstart1=donor_start_min-1;
         totalstop1=acceptor_stop_max;
         totalstart2=donor_start_min-1;
         totalstop2=acceptor_stop_max;
         end;
      end;
   else do;
      num_blocks=2;
          totalstart1=donor_start_min;
          totalstop1=acceptor_stop_max;
          totalstart2=donor_start_min;
          totalstop2=acceptor_stop_max;
      block1_length=donor_stop-donor_start_min;
      block2_length=acceptor_stop_max-acceptor_start;
      block2_start=acceptor_start-donor_start_min;
      end;
   keep junction_id chr totalstart1 totalstop1 totalstart2 totalstop2 strand score color num_blocks block1_length block2_length block1_start block2_start;
run;

data assemble_bed;
    retain chr;
    retain totalstart1;
    retain totalstop1;
    retain junction_id;
    retain score;
    retain strand;
    retain totalstart2;
    retain totalstop2;
    retain color;
    retain num_blocks;
    length lengths_cat $10.;
    length starts_cat $10.;
    set info_for_bed;
    if num_blocks=2 then do;
        lengths_cat=catx(',', block1_length, block2_length);
        starts_cat=catx(',', block1_start, block2_start);
        end;
    else if num_blocks=1 then do;
        lengths_cat=put(block1_length, 3.);
        starts_cat=put(block1_start, 1.);
        end;
    drop block1_length block2_length block1_start block2_start;
run;

proc sort data=assemble_bed;
   by chr totalstart1 totalstop1 strand;
run;



proc export data=assemble_bed
        outfile='!MCLAB/event_analysis/simulation_distinct_junc_coord.bed'
        dbms=tab replace;
        putnames=no;
        run;


