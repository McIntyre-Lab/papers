ods listing gpath="!MCLAB/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/length_figs";

/* Create a Dataset with Flags for Induced and Repressed along with Motif information */
    proc sort data = FRU.motif_search_regions;
        by primary_fbgn;
        run;

    proc sort data = FRU.Flag_ind_rep;
        by primary_fbgn;
        run;

    proc sort data = FRU.motif_flags_and_cnts;
        by primary_fbgn;
        run;

    data ind_rep_length;
        merge FRU.motif_search_regions FRU.Flag_ind_rep FRU.motif_flags_and_cnts;
        by primary_fbgn;
        run;

/* Create Histogram of the distribution of region lengths for each treatment */
    %macro create_figures(name);
        ods graphics on /reset imagename="&name";
            proc sgplot data = ind_rep_length (where=(flag_&name = 1));
                histogram region_length;
                run;
        ods graphics off;
    %mend;

    %create_figures(male_a_ind);
    %create_figures(male_b_ind);
    %create_figures(male_c_ind);
    %create_figures(male_a_rep);
    %create_figures(male_b_rep);
    %create_figures(male_c_rep);
    %create_figures(female_a_ind);
    %create_figures(female_b_ind);
    %create_figures(female_c_ind);
    %create_figures(female_a_rep);
    %create_figures(female_b_rep);
    %create_figures(female_c_rep);
    %create_figures(null_ind);
    %create_figures(null_rep);

/* Create Hisogram of the distribution of region lengths for each treatment * if it had a motif identified */
    %macro create_figures(name,letter);
        ods graphics on /reset imagename="&name._with_&letter._motif";
            proc sgplot data = ind_rep_length (where=(flag_&name = 1 and flag_fru_&letter._motif = 1));
                histogram region_length;
                run;
        ods graphics off;
    %mend;
    %create_figures(male_a_ind,a);
    %create_figures(male_b_ind,b);
    %create_figures(male_c_ind,c);
    %create_figures(male_a_rep,a);
    %create_figures(male_b_rep,b);
    %create_figures(male_c_rep,c);
    %create_figures(female_a_ind,a);
    %create_figures(female_b_ind,b);
    %create_figures(female_c_ind,c);
    %create_figures(female_a_rep,a);
    %create_figures(female_b_rep,b);
    %create_figures(female_c_rep,c);

/* Calculate the Mean and SD for Region length for each of the treatment groups */
    %macro create_mean_sd(name);
        proc means data=ind_rep_length noprint;
            class flag_&name;
            output out=tmp mean(region_length)= std(region_length)= /autoname;
            run;

        data tmp2;
            length type $20;
            retain type;
            set tmp;
            where _type_ = 1;
            type = "flag_&name";
            rename flag_&name = flag;
            drop _TYPE_;
            run;

        proc append base=mean_out data=tmp2;
        run;

        proc datasets nolist;
            delete tmp tmp2;
            run;quit;
    %mend;

    %create_mean_sd(male_a_ind);
    %create_mean_sd(male_b_ind);
    %create_mean_sd(male_c_ind);
    %create_mean_sd(male_a_rep);
    %create_mean_sd(male_b_rep);
    %create_mean_sd(male_c_rep);
    %create_mean_sd(female_a_ind);
    %create_mean_sd(female_b_ind);
    %create_mean_sd(female_c_ind);
    %create_mean_sd(female_a_rep);
    %create_mean_sd(female_b_rep);
    %create_mean_sd(female_c_rep);
    %create_mean_sd(null_ind);
    %create_mean_sd(null_rep);

    data FRU.motif_region_length_stats;
        set mean_out;
        run;

    proc datasets nolist;
        delete mean_out;
        run;quit;

/* Create a stacked dataset */

    %macro make_stack(name,num);
        data temp;
            set ind_rep_length;
            keep primary_fbgn flag_&name;
            run;

        data tmp&num;
            set temp;
            if flag_&name = 1 then trt = "flag_&name";
            else delete;
            keep primary_fbgn trt;
            run;
    %mend;
    %make_stack(male_a_ind,1);
    %make_stack(male_b_ind,2);
    %make_stack(male_c_ind,3);
    %make_stack(male_a_rep,4);
    %make_stack(male_b_rep,5);
    %make_stack(male_c_rep,6);
    %make_stack(female_a_ind,7);
    %make_stack(female_b_ind,8);
    %make_stack(female_c_ind,9);
    %make_stack(female_a_rep,10);
    %make_stack(female_b_rep,11);
    %make_stack(female_c_rep,12);
    %make_stack(null_ind,13);
    %make_stack(null_rep,14);

    data stacked;
        set tmp:;
        run;

    proc sort data=stacked;
        by primary_fbgn;
        run;

    data info;
        set ind_rep_length;
        drop flag_male_a_ind flag_male_b_ind flag_male_c_ind flag_male_a_rep
        flag_male_b_rep flag_male_c_rep flag_male_ind flag_male_rep
        flag_female_a_ind flag_female_b_ind flag_female_c_ind flag_female_a_rep
        flag_female_b_rep flag_female_c_rep flag_female_ind flag_female_rep
        flag_null_ind flag_null_rep;
        run;

    proc sort data=info;
        by primary_fbgn;
        run;

    data for_boxplot;
        merge info (in=in1) stacked (in=in2);
        by primary_fbgn;
        if in1 and not in2 then trt = 'nodiff';
        run;
    
    proc datasets nolist;
        delete temp tmp: info stacked;
        run; quit;

/* Make Lots of Box Plots */

    proc sort data=for_boxplot;
        by trt;
        run;

    /* Make Box plot of length vs treatment */
        ods graphics on / reset imagename="boxplot_length_by_trt";
        title "Region Length By Treatment";
        proc sgplot data=for_boxplot;
            vbox region_length /category=trt;
            yaxis max=150000;
            run;
        ods graphics off;

    /* Make Box plot of length vs treatment for genes that had a motif {a,b,c} identified in them */
        ods graphics on / reset imagename="boxplot_length_by_trt_with_motif";
        title "Region Length By Treatment for Genes that had a Motif";
        proc sgplot data=for_boxplot(where=(flag_fru_b_motif = 1 or flag_fru_b_motif = 1 or flag_fru_c_motif = 1));
            vbox region_length /category=trt;
            yaxis max=150000;
            run;
        ods graphics off;

    /* Make Box plot of length vs Flag Multi motif {a,b,c} to see if identifing multiple motifs is effected by length */
        ods graphics on / reset imagename="boxplot_fru_a_multi";
        title "Region Length By FruA Multi Motif Flag";
        proc sgplot data=for_boxplot;
            vbox region_length /category=flag_fru_a_multi_motif;
            yaxis max=150000;
            run;
        ods graphics off;

        ods graphics on / reset imagename="boxplot_fru_b_multi";
        title "Region Length By FruB Multi Motif Flag";
        proc sgplot data=for_boxplot;
            vbox region_length /category=flag_fru_b_multi_motif;
            yaxis max=150000;
            run;
        ods graphics off;

        ods graphics on / reset imagename="boxplot_fru_c_multi";
        title "Region Length By FruC Multi Motif Flag";
        proc sgplot data=for_boxplot;
            vbox region_length /category=flag_fru_c_multi_motif;
            yaxis max=150000;
            run;
        ods graphics off;


    /* Make Box plot of length vs Flag motif {a,b,c} to see if identifing a motifs is effected by length */
        ods graphics on / reset imagename="boxplot_fru_a_motif";
        title "Region Length By FruA Motif Flag";
        proc sgplot data=for_boxplot;
            vbox region_length /category=flag_fru_a_motif;
            yaxis max=150000;
            run;
        ods graphics off;

        ods graphics on / reset imagename="boxplot_fru_b_multi";
        title "Region Length By FruB Motif Flag";
        proc sgplot data=for_boxplot;
            vbox region_length /category=flag_fru_b_motif;
            yaxis max=150000;
            run;
        ods graphics off;

        ods graphics on / reset imagename="boxplot_fru_c_multi";
        title "Region Length By FruC Motif Flag";
        proc sgplot data=for_boxplot;
            vbox region_length /category=flag_fru_c_motif;
            yaxis max=150000;
            run;
        ods graphics off;
