/* Import the h12 peaks output for each chromosome -- These are ALL the peaks
 * for each chromosome right now */
/* I have removed the columns 5,6 from the output files prior to import here.
 * These columns are comma-delim lists of haplotype numbers, not useful */


libname dsim "!MCLAB/ethanol/Sim_Pop_Gen/sas_data";

%macro importh12(chr);

  data WORK.H12_&chr    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile
 "/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/output/H12/H12_peak_intervals/chr&chr._H12_peaks_10col.txt" delimiter='09'x MISSOVER DSD lrecl=32767 ;
        informat center best32. ;
        informat left best32. ;
        informat right best32. ;
        informat num_uniq_hap best32. ;
        informat H1 best32. ;
        informat H2 best32. ;
        informat H12 best32. ;
        informat H2H1 best32. ;
        informat smallest_edge_peak best32. ;
        informat largest_edge_peak best32. ;
        format center best12. ;
        format left best12. ;
        format right best12. ;
        format num_uniq_hap best12. ;
        format H1 best12. ;
        format H2 best12. ;
        format H12 best12. ;
        format H2H1 best12. ;
        format smallest_edge_peak best12. ;
        format largest_edge_peak best12. ;
     input
                 center
                 left
                 right
                 num_uniq_hap
                 H1
                 H2
                 H12
                 H2H1
                 smallest_edge_peak
                 largest_edge_peak
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;

%mend;

%importh12(2L);
%importh12(2R);
%importh12(3L);
%importh12(3R);
%importh12(4);
%importh12(X);





* Save the datasets;

data dsim.H12_peaks_chr2L;
  set h12_2L;
  chr="2L";
  run;

data dsim.H12_peaks_chr2R;
  set h12_2R;
  chr="2R";
  run;

data dsim.H12_peaks_chr3L;
  set h12_3L;
  chr="3L";
  run;

data dsim.H12_peaks_chr3R;
  set h12_3r;
  chr="3R";
  run;

data dsim.H12_peaks_chr4;
  set h12_4;
  chr="4";
  run;

data dsim.H12_peaks_chrX;
  set h12_X;
  chr="X";
  run;

