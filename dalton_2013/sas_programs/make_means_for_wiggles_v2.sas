/********************************************************************************
* This script will import mpileups and then calculate means for the different
* biological replicates
********************************************************************************/

* libname fru '!MCLAB/arbeitman_fru_network/sasdata';
libname sandbox '/home/jfear/sandbox/arbeitman_fru_network/sasdata';


%let piledir='/home/jfear/mclab/arbeitman_fru_network/pipeline_output/mpileup_fb530_genome';
filename mpile pipe "ls -1 &piledir | grep mpileup";

data myfiles;
    length fname $65;
    infile mpile truncover;
    input fname $65.;
    run;

data myvars;
    set myfiles;
    date=scan(fname,1,'_');
    lane=scan(fname,2,'_');
    rep=scan(fname,-3,'_');
    sample=prxchange("s!\d+-\d+-\d+_\d_(.*)_\d_fb530_genome.mpileup!$1!",99,fname);
    call symput ('num_file',_n_);
    run;

proc sort data=myvars;
    by sample rep;
    run;

%macro file_read;
    %let count=1;
    %do j=1 %to &num_file;

        data _null_;
            set myvars;
            if _n_ = &j;
            call symput ('filein',fname);
            call symput ('sample',sample);
            call symput ('rep',rep);
            run;

        data tmp;
            infile "/home/jfear/mclab/arbeitman_fru_network/pipeline_output/mpileup_fb530_genome/&filein." delimiter='09'x MISSOVER DSD lrecl=32767;
            informat chrom $32. ;
            informat pos best32. ;
            informat base $1. ;
            informat rep_&rep. best32. ;
            informat qual1 $9. ;
            informat qual2 $7. ;
            format chrom $32. ;
            format pos best32. ;
            format base $1. ;
            format rep_&rep. best32. ;
            format qual1 $9. ;
            format qual2 $7. ;
            input chrom $ pos base $ rep_&rep. qual1 $ qual2 $ ;
            drop base qual1 qual2;
            run;

        proc sort data=tmp;
            by chrom pos;
            run;

        %if %sysfunc(exist(in_&sample)) %then %do;

            %let count=&count + 1;

            data tmp2;
                set in_&sample;
                run;

            proc sort data = tmp2;
                by chrom pos;
                run;

            data in_&sample;
                retain chrom pos num_reps;
                merge tmp2 (in=in1) tmp (in=in2);
                by chrom pos;
                num_reps=&count;
                run;
        %end;
        %else %do;
            data in_&sample;
                set tmp;
                run;

            %let count=1;
        %end;
    %end;
    proc datasets;
        delete tmp;
        delete tmp2;
        run;
%mend file_read;
%file_read;

%macro take_mean(sample);
    data _null_;
        set in_&sample;
        if _n_ = 1;
        call symput ('num_reps',num_reps);
        run;

    data tmp_in_&sample;
        set in_&sample;
        run;

    %do j=1 %to &num_reps;
        data tmp_in_&sample;
            set tmp_in_&sample;
            if rep_&j = '' then rep_&j = 0;
            if &j = 1 then total_count = rep_&j;
            else total_count = total_count + rep_&j;
            run;
    %end;

    data sandbox.in_&sample._means;
        set tmp_in_&sample;
        mean = total_count / &num_reps;
        run;

    data to_export;
        set sandbox.in_&sample._means;
        keep chrom pos mean;
        rename mean=count;
        run;

    proc export data=to_export
                outfile="/home/jfear/mclab/arbeitman_fru_network/data/for_wiggles/&sample..csv"
                dbms=CSV replace;
                putnames=yes;
                run;

    proc datasets;
        delete tmp_in_&sample;
        run;
%mend take_mean;
%take_mean(48APF_CS);
%take_mean(48APF_FRUP14_440);
%take_mean(48APF_FRUW12_CHAM5);
%take_mean(AB_CS);
%take_mean(AB_FRUP14_440);
%take_mean(AB_FRUW12_CHAM5);
%take_mean(AH_BERF);
%take_mean(AH_BERM);
%take_mean(AH_CS);
%take_mean(AH_CSFEMALE);
%take_mean(AH_DSXD);
%take_mean(AH_DSXNULLF);
%take_mean(AH_DSXNULLM);
%take_mean(AH_FEMALE_FRUM_A);
%take_mean(AH_FEMALE_FRUM_B);
%take_mean(AH_FEMALE_FRUM_C);
%take_mean(AH_FRUP14_440);
%take_mean(AH_FRUW12_CHAM5);
%take_mean(AH_MALE_FRUM_A);
%take_mean(AH_MALE_FRUM_B);
%take_mean(AH_MALE_FRUM_C);
%take_mean(WPP_CS);
%take_mean(WPP_FRUP14_440);
%take_mean(WPP_FRUW12_CHAM5);
