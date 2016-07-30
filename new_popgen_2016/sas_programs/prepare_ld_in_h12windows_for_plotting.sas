/* Merge the LD results with the H12 output to get coordinates together for
 * plotting H12 vs LD */


 libname dsim "!MCLAB/ethanol/Sim_Pop_Gen/sas_data";

%macro ld_h12(chr);

 proc sort data=dsim.h12_output_chr&chr;
   by center;
   run;

proc sort data=dsim.chr&chr._ld_h12windows;
  by center;
  run;

data dsim.chr&chr._h12_and_LD;
  merge dsim.h12_output_chr&chr (in=in1) dsim.chr&chr._ld_h12windows (in=in2);
  by center;
  if in1 and in2;
  run;
%mend;

%ld_h12(2L);
%ld_h12(2R);
%ld_h12(3L);
%ld_h12(3R);
%ld_h12(4);
%ld_h12(X);

/* Add flag for quantiles */
data chr2L;
  set dsim.chr2L_h12_and_ld;
  if H12 > 0.099100 then flag_top_25 = 1; else flag_top_25 =0;
  if H12 < 0.036060 then flag_bot_25 = 1; else flag_bot_25 =0;
  if H12 < 0.099100 and H12 > 0.036060 then flag_mid_50 = 1; else flag_mid_50 =0;
  if flag_top_25=1 then flag_bin='TOP'; else if flag_mid_50=1 then flag_bin='MIDDLE'; else if flag_bot_25=1 then flag_bin='BOTTOM';
  run;
data chr2R;
  set dsim.chr2R_h12_and_ld;
  if H12 > 0.110200 then flag_top_25 = 1; else flag_top_25 =0;
  if H12 < 0.033490 then flag_bot_25 = 1; else flag_bot_25 =0;
  if H12 < 0.110200 and H12 > 0.033490 then flag_mid_50 = 1; else flag_mid_50 =0;
    if flag_top_25=1 then flag_bin='TOP'; else if flag_mid_50=1 then flag_bin='MIDDLE'; else if flag_bot_25=1 then flag_bin='BOTTOM';

  run;
data chr3L;
  set dsim.chr3L_h12_and_ld;
  if H12 > 0.094810 then flag_top_25 = 1; else flag_top_25 =0;
  if H12 < 0.032730 then flag_bot_25 = 1; else flag_bot_25 =0;
  if H12 < 0.094810 and H12 > 0.032730 then flag_mid_50 = 1; else flag_mid_50 =0;
    if flag_top_25=1 then flag_bin='TOP'; else if flag_mid_50=1 then flag_bin='MIDDLE'; else if flag_bot_25=1 then flag_bin='BOTTOM';

  run;
data chr3r;
  set dsim.chr3r_h12_and_ld;
  if H12 > 0.095570 then flag_top_25 = 1; else flag_top_25 =0;
  if H12 < 0.034330 then flag_bot_25 = 1; else flag_bot_25 =0;
  if H12 < 0.095570 and H12 > 0.034330 then flag_mid_50 = 1; else flag_mid_50 =0;
    if flag_top_25=1 then flag_bin='TOP'; else if flag_mid_50=1 then flag_bin='MIDDLE'; else if flag_bot_25=1 then flag_bin='BOTTOM';

  run;
data chr4;
  set dsim.chr4_h12_and_ld;
  if H12 > 0.4224 then flag_top_25 = 1; else flag_top_25 =0;
  if H12 < 0.2519 then flag_bot_25 = 1; else flag_bot_25 =0;
  if H12 < 0.4224 and H12 > 0.2519 then flag_mid_50 = 1; else flag_mid_50 =0;
    if flag_top_25=1 then flag_bin='TOP'; else if flag_mid_50=1 then flag_bin='MIDDLE'; else if flag_bot_25=1 then flag_bin='BOTTOM';

  run;
data chrx;
  set dsim.chrx_h12_and_ld;
  if H12 > 0.164600 then flag_top_25 = 1; else flag_top_25 =0;
  if H12 < 0.034760 then flag_bot_25 = 1; else flag_bot_25 =0;
  if H12 < 0.164600 and H12 > 0.034760 then flag_mid_50 = 1; else flag_mid_50 =0;
    if flag_top_25=1 then flag_bin='TOP'; else if flag_mid_50=1 then flag_bin='MIDDLE'; else if flag_bot_25=1 then flag_bin='BOTTOM';

  run;



%macro out(chr);
  proc export data=chr&chr
    outfile="!MCLAB/ethanol/Sim_Pop_Gen/output/LD_vs_h12/chr&chr._h12_and_LD.txt"
    DBMS=TAB REPLACE;
    run;
%mend;

%out(2L);
%out(2R);
%out(3L);
%out(3R);
%out(4);
%out(X);

