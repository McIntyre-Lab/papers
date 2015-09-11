/***** Add exon information *******/

libname splice '!MCLAB/junction_annotations/sas_data';


/* Get donor/acceptor information to add to splicing events */
/* We want exon group, gene, flag_short_exon for donors and acceptors */

/* Get info from donors */

data donor_exon_info;
   set splice.donor_exons;
   drop donor_exon_start donor_exon_stop;
run;

/* Get info from acceptors */

data acceptor_exon_info;
   set splice.acceptor_exons;
   drop acceptor_exon_start acceptor_exon_stop;
run;


/* Sort donor and acceptor exon info */

proc sort data=donor_exon_info;
   by donor_exon;
run;

proc sort data=acceptor_exon_info;
   by acceptor_exon;
run;

/* Sort AS events by donor first */

proc sort data=splice.junction_and_ir_events;
   by donor_exon;
run;


/* Merge donor exon info - for now no_donor to check that all events without a donor exon are in fact IR events! */

data splicing_events_w_donors no_donor no_event;
   merge splice.junction_and_ir_events (in=in1) donor_exon_info (in=in2);
   by donor_exon;
   if in1 and in2 then output splicing_events_W_donors;
   else if in1 then output no_donor; * remove later if non-donors are ONLY introns!;
   else output no_event; *This can be more than zero, as not all exons will be used as donors! (ie, last exon per gene, single-exon genes);
run;

/* Check that all in no_donor are introns only */

proc freq data=no_donor noprint;
   tables donor_exon /out=no_donor_check;
run;

* all are introns!;


/* Sort AS events by acceptor first */

proc sort data=splice.junction_and_ir_events;
   by acceptor_exon;
run;


/* Merge donor exon info - for now no_donor to check that all events without a donor exon are in fact IR events! */

data splicing_events_w_acceptor no_acceptor no_event;
   merge splice.junction_and_ir_events (in=in1) acceptor_exon_info (in=in2);
   by acceptor_exon;
   if in1 and in2 then output splicing_events_W_acceptor;
   else if in1 then output no_acceptor; * remove later if non-acceptor are ONLY introns!;
   else output no_event; *This can be more than zero, as not all exons will be used as acceptor! (ie, last exon per gene, single-exon genes);
run;

/* Check that all in no_donor are introns only */

proc freq data=no_acceptor noprint;
   tables acceptor_exon /out=no_acceptor_check;
run;

* all are introns!;


/* Merge donor exon info - for now no_donor to check that all events without a donor exon are in fact IR events! */


proc sort data=splice.junction_and_ir_events;
   by donor_exon;
run;


data splicing_events_w_donors no_event;
   merge splice.junction_and_ir_events (in=in1) donor_exon_info (in=in2);
   by donor_exon;
   if in1 and in2 then output splicing_events_w_donors;
   else if in1 then do;
       donor_group=.;
       donor_gene='';
       flag_short_donor=.;
       flag_alt_donor=.;
       output splicing_events_w_donors;
       end;
   else output no_event;
run;


proc sort data=splicing_events_w_donors;
   by acceptor_exon;
run;


data splicing_events_w_acceptors no_event;
   merge splicing_events_w_donors (in=in1) acceptor_exon_info (in=in2);
   by acceptor_exon;
   if in1 and in2 then output splicing_events_w_acceptors;
   else if in1 then do;
       acceptor_group=.;
       acceptor_gene='';
       flag_short_acceptor=.;
       flag_alt_acceptor=.;
       output splicing_events_w_acceptors;
       end;
   else output no_event;
run;

/* Check genes */

data as_event_gene_check;
   set splicing_events_w_acceptors;
   if donor_gene = acceptor_gene then flag_bad_gene=0;
   else do;
       if donor_gene='' and acceptor_gene ne '' then flag_bad_gene=0;
       else if donor_gene ne '' and acceptor_gene = '' then flag_bad_gene=0;
       else flag_bad_gene=1;
       end;
run;


proc freq data=as_event_gene_check;
   tables flag_bad_gene;
run;
* all good!, collapse donor_gene and acceptor_gene into one variable: gene_id;
     

/* Make permenant */

data splice.splicing_events_w_exon_info;
   set as_event_gene_check;
   length gene_id $36.;
   if donor_exon='intron' then gene_id=acceptor_gene;
   else if acceptor_exon='intron' then gene_id=donor_gene;
   else gene_id=donor_gene;
   if flag_alt_donor=. then flag_alt_donor=0;
   if flag_alt_acceptor=. then flag_alt_acceptor=0;
   if flag_short_donor=. then flag_short_donor=0;
   if flag_short_acceptor=. then flag_short_acceptor=0;
   drop acceptor_gene donor_gene flag_bad_gene;
run;

proc datasets noprint;
  delete donor_exon_info acceptor_exon_info splicing_events_w_donors as_event_gene_check
   no_donor no_event splicing_events_w_acceptor splicing_events_w_acceptors no_acceptor;
run;
quit;

