ods listing; ods html close;

libname event '!MCLAB/event_analysis/sas_data';

/* Import lists of PacBio transcripts that are only in NPCs or OLDs, or are in both.
   I am going to use these lists to calculate the "true" rate of AS, then compare
   this to the differential exon usage results */



%macro importPB(datain,dataout);

    data work.&dataout. ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile "!MCLAB/event_analysis/text_data/&datain..txt" delimiter='09'x MISSOVER
    DSD lrecl=32767 firstobs=2 ;
       informat row_num best32. ;
       informat pacbio_id $10. ;
       informat RAW_OLD1 best32. ;
       informat POST_OLD1 best32. ;
       informat RAW_OLD2 best32. ;
       informat POST_OLD2 best32. ;
       format row_num best12. ;
       format pacbio_id $10. ;
       format RAW_OLD1 best12. ;
       format POST_OLD1 best12. ;
       format RAW_OLD2 best12. ;
       format POST_OLD2 best12. ;
    input
                row_num
                pacbio_id $
                RAW_OLD1
                POST_OLD1
                RAW_OLD2
                POST_OLD2
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;

    data &dataout._2;
       set &dataout.;
       keep pacbio_id;
    run;

    proc sort data=&dataout._2;
        by pacbio_id;
    run;

%mend;

%importPB(TOTAL_NPC_exclusive,npc_only); * NPC-only list;
%importPB(TOTAL_OPC_exclusive,old_only); * OLD-only list;
%importPB(TOTAL_Intersection,npc_and_old); * NPC and OLD;

/* Merge together and flag */

data pacbio_list;
   merge npc_and_old_2 (in=in1) npc_only_2 (in=in2) old_only_2 (in=in3);
   by pacbio_id;
   if in1 then flag_in_both=1; else flag_in_both=0;
   if in2 then flag_in_npc_only=1; else flag_in_npc_only=0;
   if in3 then flag_in_old_only=1; else flag_in_old_only=0;
run;

*check to make sure that transcripts are unique to each list;

proc freq data=pacbio_list noprint;
   tables flag_in_both*flag_in_npc_only*flag_in_old_only / out=xs_check;
proc print data=xs_check;
run;

/* flag_in_    flag_in_    flag_in_
  both      npc_only    old_only    COUNT    PERCENT

    0           0           1         784     5.1766
    0           1           0         499     3.2948
    1           0           0       13862    91.5286

Good: all are exclusive to one list */

data flag_pacbio;
   set pacbio_list;
   if flag_in_both=1 then do;
      flag_npc_expressed=1;
      flag_old_expressed=1;
      end;
   else if flag_in_npc_only=1 then do;
      flag_npc_expressed=1;
      flag_old_expressed=0;
      end;
   else if flag_in_old_only=1 then do;
      flag_npc_expressed=0;
      flag_old_expressed=1;
      end;
   keep pacbio_id flag_npc_expressed flag_old_expressed;
run;

proc freq data=flag_pacbio;
  tables flag_npc_expressed*flag_old_expressed;
run;

/*
 flag_npc_expressed
           flag_old_expressed

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |      0 |    784 |    784
          |   0.00 |   5.18 |   5.18
          |   0.00 | 100.00 |
          |   0.00 |   5.35 |
 ---------+--------+--------+
        1 |    499 |  13862 |  14361
          |   3.29 |  91.53 |  94.82
          |   3.47 |  96.53 |
          | 100.00 |  94.65 |
 ---------+--------+--------+
 Total         499    14646    15145
              3.29    96.71   100.00

Total 15145 transcripts
92% common to both cell types
5% OLD exclusive
3% NPC exckusive
or 8% that are differentially detected
*/

/* Make permenant */

data event.pacbio_xscripts_flag_cell_exp;
   set flag_pacbio;
run;

