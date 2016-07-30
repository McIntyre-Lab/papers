/* H12 peak output for all chromosomes are in SAS now. I want to stack the
 * datasets together and get the top 50 peaks overall in a new dataset. I will
 * also do the top 50 for each chromosome and make a separate dataset */

 libname dsim "!MCLAB/ethanol/Sim_Pop_Gen/sas_data";

 * Combine the chromosome datasets into one big dataset, sort, and then take the
 * top 50 rows. ;

 proc sort data=dsim.H12_peaks_chr2L; by H12; run;
 proc sort data=dsim.H12_peaks_chr2R; by H12; run;
 proc sort data=dsim.H12_peaks_chr3L; by H12; run;
 proc sort data=dsim.H12_peaks_chr3R; by H12; run;
 proc sort data=dsim.H12_peaks_chrX; by H12; run;
 proc sort data=dsim.H12_peaks_chr4; by H12; run;


data H12_peaks_all_chrom;
  set dsim.H12_peaks_chr2L 
  dsim.H12_peaks_chr2R 
  dsim.H12_peaks_chr3L 
  dsim.H12_peaks_chr3R 
  dsim.H12_peaks_chrX 
  dsim.H12_peaks_chr4 ;
  run;

proc sort data=H12_peaks_all_chrom;
  by H12;
  run;


data H12_top_50_peaks_overall;
  set H12_peaks_all_chrom (obs=50);
  run;

proc export data= H12_top_50_peaks_overall
    outfile="!MCLAB/ethanol/Sim_Pop_Gen/output/H12/H12_peak_intervals/H12_top_50_intervals_overall.txt"
    DBMS=TAB REPLACE;
    run;

data dsim.H12_top_50_peaks_overall; 
  set H12_top_50_peaks_overall; 
  run;


/* Now I want the top 50 peaks from each chromosome to be in one dataset */
%macro top50(chr);

data dsim.H12_top_50_peaks_chr&chr;
  set dsim.H12_peaks_chr&chr (obs=50);
  run;

%mend;

%top50(2L);
%top50(2R);
%top50(3L);
%top50(3R);
%top50(X);
%top50(4); *Only 14 obs;


data H12_top_50_peaks_from_each_chr;
  set dsim.H12_top_50_peaks_chr2L 
  dsim.H12_top_50_peaks_chr2R 
  dsim.H12_top_50_peaks_chr3L 
  dsim.H12_top_50_peaks_chr3R 
  dsim.H12_top_50_peaks_chrX 
  dsim.H12_top_50_peaks_chr4 ;
  run;

proc sort data=H12_top_50_peaks_from_each_chr;
  by H12;
  run;

proc export data=H12_top_50_peaks_from_each_chr
    outfile="!MCLAB/ethanol/Sim_Pop_Gen/output/H12/H12_peak_intervals/H12_top50_intervals_from_each_chrom.csv"
    DBMS=CSV REPLACE;
    run;

data dsim.H12_top_50_peaks_from_each_chr;
  set H12_top_50_peaks_from_each_chr; 
  run;
