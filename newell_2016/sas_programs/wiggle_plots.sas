/* Prepare data for making wiggles. Macro at end makes wiggle plots */

libname ribo '!MCLAB/arbeitman/arbeitman_ribotag/sas_data';
libname dmel '!MCLAB/useful_dmel_data/flybase551/sasdata';


* Manipulate exon annotation ;

* Order needs to be gene_symbol chrom start end transcript_id;
data ribo.exon_anno;
  retain FBgn;
  set dmel.exon2symbol;
  drop exon_id exon_name strand FBtrs_per_exon symbol FBpp; 
  run;

* Export exon annotation to project folder;
proc export data=ribo.exon_anno
	outfile='!MCLAB/arbeitman/arbeitman_ribotag/wiggle_plots/exon_annotation.csv'
	dbms=csv replace;
	run;




* Import mpileup files. Need to average counts across bio reps. Only want chrom pos and count. ;
%macro import (ID) ;
	data WORK.pileups_&ID. ;

	%let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile "/home/fnew/ribo/&ID..txt"
 delimiter='09'x MISSOVER DSD lrecl=32767 ;
    
    informat chrom $9. ;
         informat pos best32. ;
         informat count_&ID. best32. ;
         format chrom $9. ;
         format pos best12. ;
         format count_&ID. best12. ;
      input
                  chrom $
                  pos
                  count_&ID.
      ;
      if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
      run;

proc sort data=pileups_&ID.;
  by chrom pos;
  run;


%mend;

* Import all of the mpileup files;
%import(Input_Female1);
%import(Input_Female2);
%import(Input_Female3);
%import(Input_Female4);
%import(Input_Female5);
%import(Input_Male1);
%import(Input_Male2);
%import(Input_Male3);
%import(Input_Male4);
%import(Input_Male5);
%import(Input_Male6);
%import(IP_Female1);
%import(IP_Female2);
%import(IP_Female3);
%import(IP_Female4);
%import(IP_Female5);
%import(IP_Male1);
%import(IP_Male2);
%import(IP_Male3);
%import(IP_Male4);
%import(IP_Male5);
%import(IP_Male6);


* Average counts across bio reps;
 
* Input_Female;
data input_female_all;
  merge pileups_input_female1 pileups_input_female2 pileups_input_female3 pileups_input_female4 pileups_input_female5;
  by chrom pos;
  if count_input_female1 eq " " then count_input_female1=0;
  if count_input_female2 eq " " then count_input_female2=0;
  if count_input_female3 eq " " then count_input_female3=0;
  if count_input_female4 eq " " then count_input_female4=0;
  if count_input_female5 eq " " then count_input_female5=0;
  run;
	
data input_female_all_for_wiggle;
  set input_female_all;
  avg_count_input_female = ((count_input_female1 + count_input_female2 + count_input_female3 + count_input_female4 + count_input_female5) / 5) / 0.6260817 ;
  keep chrom pos avg_count_input_female;
  run;


* Input Male;
data input_male_all;
  merge pileups_input_male1 pileups_input_male2 pileups_input_male3 pileups_input_male4 pileups_input_male5 pileups_input_male6;
  by chrom pos;
  if count_input_male1 eq " " then count_input_male1=0;
  if count_input_male2 eq " " then count_input_male2=0;
  if count_input_male3 eq " " then count_input_male3=0;
  if count_input_male4 eq " " then count_input_male4=0;
  if count_input_male5 eq " " then count_input_male5=0;
  if count_input_male6 eq " " then count_input_male6=0;
 run;

data input_male_all_for_wiggle;
  set input_male_all;
  avg_count_input_male = ((count_input_male1 + count_input_male2 + count_input_male3 + count_input_male4 + count_input_male5 + count_input_male6) / 6) / .7058248 ;
  keep chrom pos avg_count_input_male;
  run;


* IP Female;
data ip_female_all;
  merge pileups_ip_female1 pileups_ip_female2 pileups_ip_female3 pileups_ip_female4 pileups_ip_female5;
  by chrom pos;
  if count_ip_female1 eq " " then count_ip_female1=0;
  if count_ip_female2 eq " " then count_ip_female2=0;
  if count_ip_female3 eq " " then count_ip_female3=0;
  if count_ip_female4 eq " " then count_ip_female4=0;
  if count_ip_female5 eq " " then count_ip_female5=0;
  run;
	
data ip_female_all_for_wiggle;
  set ip_female_all;
  avg_count_ip_female = ((count_ip_female1 + count_ip_female2 + count_ip_female3 + count_ip_female4 + count_ip_female5) / 5) / 0.0732265 ;
  keep chrom pos avg_count_ip_female;
  run;


* IP Male;
  data ip_male_all;
  merge pileups_ip_male1 pileups_ip_male2 pileups_ip_male3 pileups_ip_male4 pileups_ip_male5 pileups_ip_male6;
  by chrom pos;
  if count_ip_male1 eq " " then count_ip_male1=0;
  if count_ip_male2 eq " " then count_ip_male2=0;
  if count_ip_male3 eq " " then count_ip_male3=0;
  if count_ip_male4 eq " " then count_ip_male4=0;
  if count_ip_male5 eq " " then count_ip_male5=0;
  if count_ip_male6 eq " " then count_ip_male6=0;
 run;

data ip_male_all_for_wiggle;
  set ip_male_all;
  avg_count_ip_male = ((count_ip_male1 + count_ip_male2 + count_ip_male3 + count_ip_male4 + count_ip_male5 + count_ip_male6) / 6) / .11127968 ;
  keep chrom pos avg_count_ip_male;
  run;



* Sort by chrom pos ;
proc sort data=input_female_all_for_wiggle;
  by chrom pos;
proc sort data=input_male_all_for_wiggle;
  by chrom pos;
proc sort data=ip_female_all_for_wiggle;
  by chrom pos;
proc sort data=ip_male_all_for_wiggle;
  by chrom pos;
  run;



* Merge counts together and replace missing values with 0 ;
data ribo_all_for_wiggle;
  merge input_female_all_for_wiggle input_male_all_for_wiggle ip_female_all_for_wiggle ip_male_all_for_wiggle ;
  by chrom pos;
  if avg_count_input_female eq " " then avg_count_input_female=0;
  if avg_count_input_male eq " " then avg_count_input_male=0;
  if avg_count_ip_female eq " " then avg_count_ip_female=0;
  if avg_count_ip_male eq " " then avg_count_ip_male=0;
  run;


* Merge only the male data ;
data ribo_male_counts;
  merge input_male_all_for_wiggle ip_male_all_for_wiggle;
  by chrom pos;
	if avg_count_ip_male
        eq " " then avg_count_ip_male=0;
	if avg_count_input_male eq " " then avg_count_input_male=0;
  run;

* Merge IP male and IP female ;
data ribo_ip_male_female;
  merge ip_female_all_for_wiggle ip_male_all_for_wiggle;
  by chrom pos;
	if avg_count_ip_male eq " " then avg_count_ip_male=0;
	if avg_count_ip_female eq " " then avg_count_ip_female=0;
  run;

* Merge Input male and Input female;	
data ribo_input_male_female;
  merge input_female_all_for_wiggle input_male_all_for_wiggle;
  by chrom pos;
	if avg_count_input_male eq " " then avg_count_input_male=0;
	if avg_count_input_female eq " " then avg_count_input_female=0;
  run;

	


* Import BED file of genes ;
 data DMEL.GENES_bed_file    ;
      %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
      infile '!MCLAB/useful_dmel_data/flybase551/output/fb551_genes.bed'
 delimiter='09'x MISSOVER DSD lrecl=32767 ;
         informat chrom $25. ;
         informat gstart best32. ;
         informat gend best32. ;
         informat gene_symbol $11. ;
         format chrom $25. ;
         format gstart best12. ;
         format gend best12. ;
         format gene_symbol $11. ;
      input
                  chrom $
                  gstart
                  gend
                  gene_symbol $
      ;
      if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
      run;


%let anno=dmel.genes_bed_file;
%let pile=ribo_ip_male_female;
%let gene='FBgn0038498';



/* Make wiggle plots */

* Macro for making wiggle plots ;
%macro merge_sym_pile(anno, pile, gene);
* Pull chrom and coordinates;
    data _null_;

set &anno(where=(gene_symbol eq &gene));
        call symput('chrom', chrom);
        call symput('gstart', gstart);
        call symput('gend', gend);
        run;

    data gene_subset;
	set &pile;
    if chrom eq "&chrom" and pos ge &gstart and pos le &gend;
        run;

    proc export data=gene_subset 
	outfile="/home/fnew/mclab/arbeitman/arbeitman_ribotag/wiggle_plots/gene_counts_file.csv"
	dbms=csv replace;
	run;


data _null_;

call system("Rscript !MCLAB/scripts/R/wiggleplots_example.R !MCLAB/arbeitman/arbeitman_ribotag/wiggle_plots/exon_annotation.csv !MCLAB/arbeitman/arbeitman_ribotag/wiggle_plots/gene_counts_file.csv &gene !MCLAB/arbeitman/arbeitman_ribotag/wiggle_plots/&gene._input_male_female.png");
    run;

%mend;

*Input here the genes of interest;
* Compare all four treatments;
%merge_sym_pile(dmel.genes_bed_file, ribo_all_for_wiggle, 'FBgn0004652');
%merge_sym_pile(dmel.genes_bed_file, ribo_all_for_wiggle, 'FBgn0015381');
%merge_sym_pile(dmel.genes_bed_file, ribo_all_for_wiggle, 'FBgn0005631');
%merge_sym_pile(dmel.genes_bed_file, ribo_all_for_wiggle, 'FBgn0051632');
%merge_sym_pile(dmel.genes_bed_file, ribo_all_for_wiggle, 'FBgn0038953');
%merge_sym_pile(dmel.genes_bed_file, ribo_all_for_wiggle, 'FBgn0000504');

* Compare the input and ip male;
%merge_sym_pile(dmel.genes_bed_file, ribo_male_counts, 'FBgn0004652');
%merge_sym_pile(dmel.genes_bed_file, ribo_male_counts, 'FBgn0015381');
%merge_sym_pile(dmel.genes_bed_file, ribo_male_counts, 'FBgn0005631');
%merge_sym_pile(dmel.genes_bed_file, ribo_male_counts, 'FBgn0051632');
%merge_sym_pile(dmel.genes_bed_file, ribo_male_counts, 'FBgn0038953');
%merge_sym_pile(dmel.genes_bed_file, ribo_male_counts, 'FBgn0000504');

* Compare IP male and IP female;
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0004652');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0015381');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0005631');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0051632');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0038953');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0000504');


* Compare Input male and Input female;
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0004652');
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0015381');
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0005631');
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0051632');
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0038953');
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0000504');
*Input here the genes of interest;
* Compare all four treatments;
%merge_sym_pile(dmel.genes_bed_file, ribo_all_for_wiggle, 'FBgn0004652');
%merge_sym_pile(dmel.genes_bed_file, ribo_all_for_wiggle, 'FBgn0015381');
%merge_sym_pile(dmel.genes_bed_file, ribo_all_for_wiggle, 'FBgn0005631');
%merge_sym_pile(dmel.genes_bed_file, ribo_all_for_wiggle, 'FBgn0051632');
%merge_sym_pile(dmel.genes_bed_file, ribo_all_for_wiggle, 'FBgn0038953');
%merge_sym_pile(dmel.genes_bed_file, ribo_all_for_wiggle, 'FBgn0000504');

* Compare the input and ip male;
%merge_sym_pile(dmel.genes_bed_file, ribo_male_counts, 'FBgn0004652');
%merge_sym_pile(dmel.genes_bed_file, ribo_male_counts, 'FBgn0015381');
%merge_sym_pile(dmel.genes_bed_file, ribo_male_counts, 'FBgn0005631');
%merge_sym_pile(dmel.genes_bed_file, ribo_male_counts, 'FBgn0051632');
%merge_sym_pile(dmel.genes_bed_file, ribo_male_counts, 'FBgn0038953');
%merge_sym_pile(dmel.genes_bed_file, ribo_male_counts, 'FBgn0000504');

* Compare IP male and IP female;
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0004652');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0015381');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0005631');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0051632');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0038953');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0000504');


/*
Genes: (from Michelle)

1)fruitless, all exons, since some are enriched in female and some in male. fru CG14307 FBgn0004652

2) dsf (F3320_SI) CG9019  FBgn0015381

3) robo (S46170_SI)	CG13521  FBgn0005631

4) sens-2 (S3938_SI)	CG31632  FBgn0051632

5) CG18596 (S57971_SI)	FBgn0038953

6) DSX - FBgn0000504 CG11094 (adding this in!)

*/


/* 11/30/15, new genes of interest */

* Macro for making wiggle plots ;
%macro merge_sym_pile(anno, pile, gene);
* Pull chrom and coordinates;
    data _null_;

set &anno(where=(gene_symbol eq &gene));
        call symput('chrom', chrom);
        call symput('gstart', gstart);
        call symput('gend', gend);
        run;

    data gene_subset;
	set &pile;
    if chrom eq "&chrom" and pos ge &gstart and pos le &gend;
        run;

    proc export data=gene_subset 
	outfile="/home/fnew/ufgi_share/SHARE/McIntyre_Lab/arbeitman/arbeitman_ribotag/wiggle_plots/gene_counts_file.csv"
	dbms=csv replace;
	run;


data _null_;

call system("Rscript /home/fnew/ufgi_share/SHARE/McIntyre_Lab/scripts/R/wiggleplots_example.R /home/fnew/ufgi_share/SHARE/McIntyre_Lab/arbeitman/arbeitman_ribotag/wiggle_plots/exon_annotation.csv /home/fnew/ufgi_share/SHARE/McIntyre_Lab/arbeitman/arbeitman_ribotag/wiggle_plots/gene_counts_file.csv &gene /home/fnew/ufgi_share/SHARE/McIntyre_Lab/arbeitman/arbeitman_ribotag/wiggle_plots/&gene._Input_male_female.png");
    run;

%mend;



* Compare IP male and IP female;
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0038498');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0004828');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0034012');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0026262');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0000479');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0003091');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0003475');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0004652');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0014020');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0020653');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0023023');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0029687');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0034970');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0037643');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0039924');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0261113');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0261238');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0261570');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0263397');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0264357');


* Compare Input male and Input female;
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0004652');
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0015381');
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0005631');
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0051632');
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0038953');
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0000504');



* Additional genes for wiggles from Nikki;
* 3/1/16;

* Compare IP Male and IP female;
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0043005');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0000592');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0010240');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0259240');
%merge_sym_pile(dmel.genes_bed_file, ribo_ip_male_female, 'FBgn0035649');


*Compare Input male and Input female;
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0265991');
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0031752');
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0015269');
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0262593');
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0261794');
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0035464');
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0036480');
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0263998');
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0001139');
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0036888');
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0259734');
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0030648');
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0001215');
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0003178');
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0034971');
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0000479');
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0033240');
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0263097');
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0001612');
%merge_sym_pile(dmel.genes_bed_file, ribo_input_male_female, 'FBgn0052683');

