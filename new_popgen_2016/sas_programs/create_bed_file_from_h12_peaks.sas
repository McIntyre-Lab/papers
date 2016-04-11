/* Create a bed file of the coordinates of the H12 peaks*/



 data WORK.COORD    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile '/home/fnew/dsim/h12/chr4_h12_window_coords.txt' delimiter='09'x MISSOVER DSD
 lrecl=32767 firstobs=2 ;
        informat center best32. ;
        informat start best32. ;
        informat end best32. ;
        format center best12. ;
        format start best12. ;
        format end best12. ;
     input
                 center
                 start
                 end
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;


data coord2;
  set coord;
  chr="4";
   run;

data coord3;
  retain chr start end center;
  set coord2;
  peak="peak";
  run;

data coord4;
  set coord3;
  center1=catx('_',vvalue(peak),vvalue(center));
  run;

data coord5;
  set coord4;
  drop center peak;
  rename center1=center;
  run;


proc export data=coord5
    outfile="/home/fnew/dsim/h12/chr4_h12_peak.bed"
    dbms=TAB REPLACE;
    run;
