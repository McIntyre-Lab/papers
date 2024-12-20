/********************************************************************************
* This script imports coverage counts against fusions, from the CSV files and
* creates a permanent stacked dataset.
********************************************************************************/

/* Import Complete Coverage Counts */

    %macro import_cc(line,mv,rep);
        %let name="!MCLAB/cegs_sergey/pipeline_output/OE_coverage_counts/coverage_count_fb551_non-redundant_fusions_on_fusions_nodup_platek/&line._&mv.&rep..csv";

        data coverage_counts    ;
            infile  &name delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
            informat fusion_id $9. ;
            informat mapped_reads best32. ;
            informat read_length best32. ;
            informat region_length best32. ;
            informat region_depth best32. ;
            informat reads_in_region best32. ;
            informat apn best32. ;
            informat rpkm best32. ;
            format fusion_id $9. ;
            format mapped_reads best12. ;
            format read_length best12. ;
            format region_length best12. ;
            format region_depth best12. ;
            format reads_in_region best12. ;
            format apn best12. ;
            format rpkm best12. ;
            input
                 fusion_id $
                 mapped_reads
                 read_length
                 region_length
                 region_depth
                 reads_in_region
                 apn
                 rpkm
            ;
            run;

        data &line._&mv.&rep._cc;
            retain line mating_status rep;
            set coverage_counts;
            length line $4.;
            length mating_status $1.;
            length mating_status $1.;
            length rep $2.;
            line = "&line";
            mating_status = "&mv";
            rep = "&rep.k";
            run;

        proc append base=all_cc data=&line._&mv.&rep._cc;
            run;

        proc datasets nolist;
            delete coverage_counts ;
            delete &line._&mv.&rep._cc;
            run;
            quit;

    %mend import_cc;
    %iterdataset(dataset=CEGS.platek_design,function=%nrstr(%import_cc(&line,&mating_status,&rep);));

/* Create Permanant Dataset */
    data CEGS.platek_ccfus_stack;
        set all_cc;
        sample = trim(line) || '_' || trim(mating_status) || '_' || trim(rep);
        run;

/* Create a Local Copy for Faster Access */
    data CEGLOCAL.platek_ccfus_stack;
        set CEGS.platek_ccfus_stack;
        run;

/* Clean Up */
    proc datasets nolist;
        delete all_cc;
        run; quit;
