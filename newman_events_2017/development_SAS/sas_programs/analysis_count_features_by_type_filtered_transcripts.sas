ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Count the number of exons, fragments and junctions that are detected and/or unique for each transcript (detection threshold >0%). I also want to determine the number of transcripts that have unique exons, fragments and junctions, and what the overlap is between these sets (for making Venn diagrams)

  i.e. what transcripts have a unique singleton fragment, what has a unique fragment, what have a unique event, and what is overlap */

/* Set threshold here. I will make this into a command line variable later */
%let percUniqOn=0;
*%let percUniqOn=0.25;
*%let percUniqOn=0.5;
*%let percUniqOn=0.75;
*%let percUniqOn=1;

/* Calculate suffix */
%let percSuffIn=;
%let percSuffOut=;

data _null_;
   if &percUniqOn = 0 then call symput('percSuffIn',25);
   else call symput('percSuffIn', %SYSEVALF(&percUniqOn * 100));
   call symput('percSuffOut', %SYSEVALF(&percUniqOn * 100));
run;

%put &percSuffIn.;
%put &percSuffOut.;
%put &percUniqOn.;


%let type=color;

data _null_;
  if "&type" = 'color' then call symput('varlist','red blue');
  if "&type" = 'pattern' then call symput('varlist','solid stripe');
run;



data data1 (keep= &varlist);
  set datalib.data1;
run;

data _null_;
   if &percUniqOn. = 0 %then %let percSuff = 25;
   else %let percSuff = %SYSEVALF(&percUniqOn * 100);
run;

%put &percSuff.;


data xscript_list;
  set event.xscripts_w_uniq_dtct_gt_25;


