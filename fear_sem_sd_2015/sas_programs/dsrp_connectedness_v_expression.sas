/********************************************************************************
* Scrip to look at the relationship of averge gene expression and the level of
* connectedness that a GGM identifies. 
********************************************************************************/

/* Data munge cgnumbers */
    * In GGMs I replaced sex det gene's cgnumber with their symbol. I need to do
    * the same here so I can merge on the ggm results.
    ;

    data stack;
        set SEM.dsrp_stack;
        num = index(cgnumber, '_');
        if num > 0 then suffix = substr(cgnumber,num+1);
        else suffix = '';
        drop FBgn FBtr;
        run;

    data stack_w_sex;
        set stack;
            if 
            symbol eq 'B52' or 
            symbol eq 'dsx' or 
            symbol eq 'fl(2)d' or 
            symbol eq 'fru' or 
            symbol eq 'her' or 
            symbol eq 'ix' or 
            symbol eq 'msl-2' or 
            symbol eq 'mub' or 
            symbol eq 'ps'  or
            symbol eq 'Psi' or 
            symbol eq 'Rbp1' or 
            symbol eq 'Rm62' or 
            symbol eq 'snf' or 
            symbol eq 'Spf45' or 
            symbol eq 'sqd' or 
            symbol eq 'Sxl' or 
            symbol eq 'tra' or 
            symbol eq 'tra2' or 
            symbol eq 'vir' or 
            symbol eq 'Yp1' or 
            symbol eq 'Yp2' or 
            symbol eq 'Yp3' then do;
                if suffix ne '' and symbol eq 'fl(2)d' then sym = 'fl_2_d_' || trim(suffix);
                else if suffix ne '' then sym = trim(symbol) || '_' || trim(suffix);
                else sym = symbol;
            end;
            else sym = cgnumber;
            drop num suffix symbol cgnumber;
            run;

/* Calculate Averge expression on a per gene level */
    proc sort data=stack_w_sex;
        by sym;
        run;

    proc means data=stack_w_sex noprint;
        by sym;
        output mean(col1)=mean_exp std(col1)=std_exp sum(col1)=sum_exp out=means;
        run;

    data means2;
        set means;
        rename sym = FG;
        drop _type_ _freq_;
        run;

/* Import GGM Neighborhood */
    data WORK.NEIGH    ;
    infile '!MCLAB/cegs_sem_sd_paper/analysis_output/ggm/dsrp_neighborhood_analysis_v3.csv' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
      informat FG $54. ;
      informat FG_sex best32. ;
      informat FG_splice best32. ;
      informat num_primary best32. ;
      informat num_primary_sex best32. ;
      informat num_secondary best32. ;
      informat num_secondary_sex best32. ;
      informat total_neighbor best32. ;
      informat total_neighbor_sex best32. ;
      format FG $54. ;
      format FG_sex best12. ;
      format FG_splice best12. ;
      format num_primary best12. ;
      format num_primary_sex best12. ;
      format num_secondary best12. ;
      format num_secondary_sex best12. ;
      format total_neighbor best12. ;
      format total_neighbor_sex best12. ;
    input
               FG $
               FG_sex
               FG_splice
               num_primary
               num_primary_sex
               num_secondary
               num_secondary_sex
               total_neighbor
               total_neighbor_sex
    ;
    run;

/* Merge Neighborhood and Expression information */
    proc sort data=means2;
        by FG;
        run;

    proc sort data=Neigh;
        by FG;
        run;

    data merged oops;
        merge means2 (in=in1) Neigh (in=in2);
        by FG;
        if in1 and in2 then output merged;
        else output oops;
        run;
    
/* Export to external dataset for import to R */
    proc export data=merged outfile='/tmp/neigh_means.csv' dbms=csv replace;
        putnames=yes;
        run;

/* Use R to Create Figures */
    * Paste the following code in R;
    /*
        mclab <- Sys.getenv("MCLAB")

        library(ggplot2)
        source(paste0(mclab,'/scripts/R/multiplot.R'))
        mydata <- read.csv('/tmp/neigh_means.csv')


        p1 <- qplot(data=mydata, x=total_neighbor, y=mean_exp, geom='density2d', main="Total in Neighborhood")
        p2 <- qplot(data=mydata, x=total_neighbor_sex, y=mean_exp, geom='density2d', main="Total SEX in Neighborhood")

        p3 <- qplot(data=mydata, x=num_primary, y=mean_exp, geom='density2d', main="Number in Primary")
        p4 <- qplot(data=mydata, x=num_primary_sex, y=mean_exp, main="Number SEX in Primary")

        p5 <- qplot(data=mydata, x=num_secondary, y=mean_exp, geom='density2d', main="Number in Secondary")
        p6 <- qplot(data=mydata, x=num_secondary_sex, y=mean_exp, main="Number SEX in Secondary")

        png(paste0(mclab,'/cegs_sem_sd_paper/analysis_output/ggm/expression_v_neighbohood.png'), width=1200, height=800)
            multiplot(p1, p3, p5, p2, p4, p6, cols=2)
        dev.off()


        png(paste0(mclab,'/cegs_sem_sd_paper/analysis_output/ggm/expression_v_std.png'))
            qplot(data=mydata, x=mean_exp, y=std_exp, geom='density2d', main="Variation and Avg Expression")
        dev.off()

    */
