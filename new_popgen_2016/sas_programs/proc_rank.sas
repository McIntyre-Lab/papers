/* proc rank and take mean of LD in bins of equal size by distance */

libname dsim "!MCLAB/ethanol/Sim_Pop_Gen/sas_data";




%macro ranking(chr);
data ld&chr;
  set dsim.ld&chr;
  distance =dist;
  run;
proc sort data=ld&chr;
  by dist;
  run;

proc rank data=ld&chr out=rank_&chr groups=1000 ;
  var dist;
  run;

proc means data=rank_&chr noprint;
  by dist;
  var R2;
  output out=rank_means_&chr mean=mean;
  run;

proc means data=rank_&chr noprint;
  by dist;
  var distance;
  output out=rank_median_&chr median=median;
  run;

data rank_mean_med_&chr;
  merge rank_means_&chr rank_median_&chr;
  by dist;
  run;


%mend;

%ranking(2l);
%ranking(2r);
%ranking(3l);
%ranking(3r);
%ranking(4);
%ranking(x);

