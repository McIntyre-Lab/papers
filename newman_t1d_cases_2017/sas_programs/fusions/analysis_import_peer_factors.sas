/* Import PEER factors */

libname con '/home/jrbnewman/concannon/sas_data';

/* Import PEER factors */
    data WORK.PEER_FACTORS   ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile '/home/jrbnewman/concannon/peer/peer_factors.csv' delimiter = ',' MISSOVER DSD
lrecl=32767 firstobs=1 ;
       informat p1 best32. ;
       informat p2 best32. ;
       informat p3 best32. ;
       informat p4 best32. ;
       informat p5 best32. ;
       informat p6 best32. ;
       informat p7 best32. ;
       informat p8 best32. ;
       informat p9 best32. ;
       informat p10 best32. ;
       format p1 best12. ;
       format p2 best12. ;
       format p3 best12. ;
       format p4 best12. ;
       format p5 best12. ;
       format p6 best12. ;
       format p7 best12. ;
       format p8 best12. ;
       format p9 best12. ;
       format p10 best12. ;
    input
                p1
                p2
                p3
                p4
                p5
                p6
                p7
                p8
                p9
                p10
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;

/* Add sample names back -- the order of samples for factors is the same as the order of samples from input data */

data sample_ids;
  set con.design_file;
  if name =  '2009-PC-0221' then delete; *sample 75 cd8;
  if name =  '2009-PC-0144' then delete; *sample 48 cd4;
  if name =  '2009-PC-0236' then delete; *sample 80;
  if name =  '2009-PC-0237' then delete; *sample 80;
  if name =  '2009-PC-0235' then delete; *sample 80;
  keep Name;
run;

proc sort data=sample_ids nodup;
   by Name;
run;

data peer_factors_w_names;
  merge sample_ids peer_factors;
run;

/* Make permenant */

data con.peer_factors;
   set peer_factors_w_names;
run;




