/* First import the gene lists. Then merge the lists into the full dataset, add flags               */
/* Lists: 1253 genes union Fru ABC dalton/fear. These are the genes that are regulated by Fru.      */
/*        3121 genes male IP vs female IP (fdr20), these are the genes that are male biased.        */
/*        575 these are the genes that overlap between the male biased and the fru regulated genes. */


libname ribo '!MCLAB/arbeitman/arbeitman_ribotag/sas_data';
libname dmel '!MCLAB/useful_dmel_data/flybase551/sasdata';

/* Import Motif Data */
%macro import_motif(indata,outdata);

        proc import datafile=&indata out=&outdata dbms=CSV replace;
            getnames=yes;
            guessingrows=1000; 
            run;

        proc sort data=&outdata;
            by primary_fbgn;
            run;

%mend import_motif;

    %import_motif('!MCLAB/arbeitman/arbeitman_ribotag/motif_analysis/fru_a_results_up_and_down.csv',fru_a);
    %import_motif('!MCLAB/arbeitman/arbeitman_ribotag/motif_analysis/fru_b_results_up_and_down.csv',fru_b);
    %import_motif('!MCLAB/arbeitman/arbeitman_ribotag/motif_analysis/fru_c_results_up_and_down.csv',fru_c);

/*Combine into full dataset*/
%macro create_full_dataset(letter);

        data fru_&letter.2;
            set fru_&letter;
            rename motif_count = fru_&letter._motif_cnt;
            rename motif_positions = fru_&letter._motif_pos;
            rename motif_eval = fru_&letter._motif_eval;
            if motif_count = 0  then flag_fru_&letter._motif = 0; else flag_fru_&letter._motif = 1;
            if motif_count = 0  then flag_fru_&letter._multi_motif = 0;
            if motif_count > 0  and  motif_count <= 5  then flag_fru_&letter._multi_motif = 1;
            if motif_count > 5  and  motif_count <= 10 then flag_fru_&letter._multi_motif = 2;
            if motif_count > 10 then flag_multi_motif = 3;
            keep primary_fbgn motif_count motif_positions motif_eval flag_fru_&letter._motif flag_fru_&letter._multi_motif;
            run;

        proc sort data=fru_&letter.2;
            by primary_fbgn;
            run;

%mend create_full_dataset;

    %create_full_dataset(a); * 16379 obs;
    %create_full_dataset(b); * 16379 obs;
    %create_full_dataset(c); * 16379 obs;

data ribo.motif_flags_and_cnts;
    merge fru_a2 (in=in1) fru_b2 (in=in2) fru_c2 (in=in3);
    by primary_fbgn;
    run;

/* Clean up */
proc datasets nolist;
    delete fru_a fru_a2 fru_b fru_b2 fru_c fru_c2;
    run; quit;




/* Import gene lists from Michelle */
 data WORK.GENE_LIST575    ;
      %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
      infile '/home/fnew/mclab/arbeitman/arbeitman_ribotag/from_michelle/575genesintersectIPandDalton.csv' delimiter =
 ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
         informat Gene_Secondary_Identifier $7. ;
         informat symbol $8. ;
         informat primary_fbgn $11. ;
         informat Organism $23. ;
         format Gene_Secondary_Identifier $7. ;
         format symbol $8. ;
         format primary_fbgn $11. ;
         format Organism $23. ;
      input
                  Gene_Secondary_Identifier $
                  symbol $
                  primary_fbgn $
                  Organism $
      ;
      if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
      run;


 data WORK.GENE_LIST1253    ;
      %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
      infile '/home/fnew/mclab/arbeitman/arbeitman_ribotag/from_michelle/1253genesUnionFruABCDaltonandFear.csv' delimiter =
 ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
         informat Gene_Secondary_Identifier $7. ;
         informat symbol $8. ;
         informat primary_fbgn $11. ;
         informat Organism $23. ;
         format Gene_Secondary_Identifier $7. ;
         format symbol $8. ;
         format primary_fbgn $11. ;
         format Organism $23. ;
      input
                  Gene_Secondary_Identifier $
                  symbol $
                  primary_fbgn $
                  Organism $
      ;
      if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
      run;



 data WORK.GENE_LIST3121    ;
      %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
      infile '/home/fnew/mclab/arbeitman/arbeitman_ribotag/from_michelle/3121genesmaleIPvsFemaleIPFDR20.csv' delimiter =
 ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
         informat Gene_Secondary_Identifier $7. ;
         informat symbol $8. ;
         informat primary_fbgn $11. ;
         informat Organism $23. ;
         format Gene_Secondary_Identifier $7. ;
         format symbol $8. ;
         format primary_fbgn $11. ;
         format Organism $23. ;
      input
                  Gene_Secondary_Identifier $
                  symbol $
                  primary_fbgn $
                  Organism $
      ;
      if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
      run;


* Sort the gene lists;

proc sort data=gene_list1253;
  by primary_fbgn;
  run;
proc sort data=gene_list575;
  by primary_fbgn;
  run;
proc sort data=gene_list3121;
  by primary_fbgn;
  run;



/*******************************************************************************/
/* Now merge the gene lists with the motif flags and counts */


* 575 - genes intersect IP and Dalton list ;

data list575_fru_abc ;
  merge gene_list575 (in=in1) ribo.motif_flags_and_cnts (in=in2);
  by primary_fbgn;
  if in1;
  if in1 and in2 then flag_motif=1; else flag_motif=0;
  run;

data list575_fru_abc;
  set list575_fru_abc;
  array a(*) _numeric_ ;
  do i=1 to dim(a);
  if a(i) = . then a(i) =0;
  end;
  drop i;
  run;


* 3121 - genes male IP vs female IP at FDR .20 ;

data list3121_fru_abc ;
  merge gene_list3121 (in=in1) ribo.motif_flags_and_cnts (in=in2);
  by primary_fbgn;
  if in1;
  if in1 and in2 then flag_motif=1; else flag_motif=0;
  run;

data list3121_fru_abc;
  set list3121_fru_abc;
  array a(*) _numeric_ ;
  do i=1 to dim(a);
  if a(i) = . then a(i) =0;
  end;
  drop i;
  run;



* 1253 - genes union FRU ABC Dalton and Fear ;
data list1253_fru_abc only1253;
  merge gene_list1253 (in=in1) ribo.motif_flags_and_cnts (in=in2);
  by primary_fbgn;
  if in1 and in2 then output list1253_fru_abc;
  if in1 and not in2 then output only1253;
  run;


data list1253_fru_abc;
  merge gene_list1253 (in=in1) ribo.motif_flags_and_cnts (in=in2);
  by primary_fbgn;
  if in1;
  if in1 and in2 then flag_motif=1; else flag_motif=0;
  run;

data list1253_fru_abc;
  set list1253_fru_abc;
  array a(*) _numeric_ ;
  do i=1 to dim(a);
  if a(i) = . then a(i) =0;
  end;
  drop i;
  run;

* Make datasets permanent ;

data ribo.intersectIPandDalton_575;
  set list575_fru_abc;
  run;

data ribo.unionFruABC_DaltonandFear_1253;
  set list1253_fru_abc;
  run;

data ribo.maleIPvsFemaleIPFDR20_3121;
  set list3121_fru_abc;
  run;
