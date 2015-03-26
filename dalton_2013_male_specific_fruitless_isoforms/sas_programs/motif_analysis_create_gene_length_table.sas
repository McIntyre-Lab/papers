/********************************************************************************
* Calculate the region length that genes were searched for motifs
********************************************************************************/

*libname FRU '!MCLAB/arbeitman/arbeitman_fru_network/sasdata';
*libname DMEL '!MCLAB/useful_dmel_data/flybase530/sasdata';

/* Create dataset with region legnth and motif flags */
proc sort data=DMEL.fbgn2coord;
    by primary_fbgn;
    run;

data FRU.motif_search_regions;
    set DMEL.fbgn2coord;
    region_start = start - 2000;
    region_end = end + 2000;
    region_length = region_end - region_start;
    keep primary_fbgn region_start region_end region_length;
    run;
