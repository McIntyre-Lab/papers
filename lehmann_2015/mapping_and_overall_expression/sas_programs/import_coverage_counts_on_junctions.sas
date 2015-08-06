/********************************************************************************
* This script imports coverage counts against fusions, from the CSV files and
* creates a permanent stacked dataset.
********************************************************************************/

/* Running all of the datasets takes a really long time, so if you are testing
 * things then change this temp design file to have only a few datasets */

/* Import the number of mapped reads from cannonical junctions Complete */
    data tmp_design;
        set CEGS.complete_design_by_rep;
        run;

    %macro import_cc2(line,mv,rep);
        %let name = "!MCLAB/cegs_sergey/pipeline_output/OE_coverage_counts/coverage_count_fb551_canonical_junctions_nodup/&line._&mv.&rep..csv";

        data coverage_counts    ;
            infile  &name delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
            informat fusion_id $90. ;
            informat junc_mapped_reads best32. ;
            informat junc_read_length best32. ;
            informat junc_region_length best32. ;
            informat junc_region_depth best32. ;
            informat junc_reads_in_region best32. ;
            informat junc_apn best32. ;
            informat junc_rpkm best32. ;
            format fusion_id $90. ;
            format junc_mapped_reads best12. ;
            format junc_read_length best12. ;
            format junc_region_length best12. ;
            format junc_region_depth best12. ;
            format junc_reads_in_region best12. ;
            format junc_apn best12. ;
            format junc_rpkm best12. ;
            input
                 fusion_id $
                 junc_mapped_reads
                 junc_read_length
                 junc_region_length
                 junc_region_depth
                 junc_reads_in_region
                 junc_apn
                 junc_rpkm
            ;
            run;

        data &line._&mv.&rep._junc;
            retain line mating_status rep;
            set coverage_counts;
            where junc_apn gt 0;
            length line $4.;
            length mating_status $1.;
            length mating_status $1.;
            length rep $1.;
            line = "&line";
            mating_status = "&mv";
            rep = "&rep";
            keep line mating_status rep fusion_id junc_mapped_reads junc_reads_in_region;
            run;

        proc append base=all_junc data=&line._&mv.&rep._junc;
            run;

        proc datasets nolist;
            delete coverage_counts ;
            delete &line._&mv.&rep._junc;
            run;
            quit;

    %mend import_cc2;
    %iterdataset(dataset=tmp_design,function=%nrstr(%import_cc2(&line,&mating_status,&rep);));

/* Import the number of mapped reads from cannonical junctions Incomplete */
    data tmp_design;
        set CEGS.incomplete_design_by_rep;
        run;

    %macro import_cc2(line,mv,rep);
        %let name = "!MCLAB/cegs_sergey/pipeline_output/OE_coverage_counts/coverage_count_fb551_canonical_junctions_incomplete_nodup/&line._&mv.&rep..csv";

        data coverage_counts    ;
            infile  &name delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
            informat fusion_id $90. ;
            informat junc_mapped_reads best32. ;
            informat junc_read_length best32. ;
            informat junc_region_length best32. ;
            informat junc_region_depth best32. ;
            informat junc_reads_in_region best32. ;
            informat junc_apn best32. ;
            informat junc_rpkm best32. ;
            format fusion_id $90. ;
            format junc_mapped_reads best12. ;
            format junc_read_length best12. ;
            format junc_region_length best12. ;
            format junc_region_depth best12. ;
            format junc_reads_in_region best12. ;
            format junc_apn best12. ;
            format junc_rpkm best12. ;
            input
                 fusion_id $
                 junc_mapped_reads
                 junc_read_length
                 junc_region_length
                 junc_region_depth
                 junc_reads_in_region
                 junc_apn
                 junc_rpkm
            ;
            run;

        data &line._&mv.&rep._junc;
            retain line mating_status rep;
            set coverage_counts;
            where junc_apn gt 0;
            length line $4.;
            length mating_status $1.;
            length mating_status $1.;
            length rep $1.;
            line = "&line";
            mating_status = "&mv";
            rep = "&rep";
            keep line mating_status rep fusion_id junc_mapped_reads junc_reads_in_region;
            run;

        proc append base=all_junc data=&line._&mv.&rep._junc;
            run;

        proc datasets nolist;
            delete coverage_counts ;
            delete &line._&mv.&rep._junc;
            run;
            quit;

    %mend import_cc2;
    %iterdataset(dataset=tmp_design,function=%nrstr(%import_cc2(&line,&mating_status,&rep);));

/* Create Permanant Dataset */
    data CEGS.junc_cnts;
        set all_junc;
        run;

/* Create Local Dataset for Faster Access */
    data CEGLOCAL.junc_cnts;
        set all_junc;
        run;

/* Clean Up */
    proc datasets nolist;
        delete all_junc;
        delete tmp_design;
        run; quit;
