/* Testing my iterative merge macro */

libname eqtl '/mnt/data/eqtls/sas_data';

/* Create a temporary empty dataset with proper formats */
data temp_fus;
   set eqtl.fusion_counts_for_eqtls;
   if _n_ = 1;
   keep Name cell_type subject_id fusion_id log_q3_q3_apn gene_id;
   rename fusion_id=feature_id log_q3_q3_apn=log_measurement;
run;

data temp_index;
   set eqtl.snp2gene_index;
   if _n_ = 1;
run;

data temp_geno;
   set eqtl.genotype_data;
   if _n_ = 1;
run;

proc sort data=temp_index;
   by gene_id snp_id;
proc sort data=temp_geno;
   by gene_id snp_id;
run;

data temp_geno_index;
   merge temp_index (in=in1) temp_geno(in=in2);
   by gene_id snp_id;
   if in1 and in2;
run;
 
proc sort data=temp_geno_index;
   by gene_id subject_id;
proc sort data=temp_fus;
   by gene_id subject_id;
run;


data eqtl_data_stacked temp;
   merge temp_geno_index (in=in1) temp_fus (in=in2);
   by gene_id subject_id;
   if in1 and in2 then output temp;
run;

proc datasets noprint;
  delete temp;
run;quit;

/* Make a test sets for testing the macro */

data test_index;
   set eqtl.snp2gene_index;
   if _n_ le 8;
run;

data test_genotypes;
   set eqtl.genotype_data;
run;

data test_exp;
   set eqtl.fusion_counts_for_eqtls;
   keep Name cell_type subject_id fusion_id log_q3_q3_apn gene_id;
   rename fusion_id=feature_id log_q3_q3_apn=log_measurement;
run;


/* This is outside the macro since I shouldn't need to sort it each time it gets subsetted (time saver)! */        
        proc sort data=test_exp;
             by gene_id subject_id;
             run;


/* Iterate through SNP-Gene index to merge SNP data and expression data */

%macro loopOverDatasets();
    /* declare macro variables */
    %local iter inSNP inGene;

    /*initiate loop*/
    %let iter=1;
    %do %while (&iter.<= 8);

       /* Sort snp-gene index of iteration and all genotypes, then merge */

        /*get libref and dataset name for dataset you will work on during this iteration*/
        data _NULL_;
            set test_index (firstobs=&iter. obs=&iter.); *only read 1 record;
            *write the Name and Event to the macro variables;
            call symput("inSNP",strip(snp_id));
            call symput("inGene",strip(gene_id));
            *NOTE: strip function is to remove any trailing blanks. This should not be a problem, but it is here just in case;
        run;

        /* get SNP-gene pair */
        data temp_snp; *assuming you want to apply the changes to the dataset itself;
            set test_genotypes;
            if snp_id = "&inSNP." and gene_id="&inGene." then output;
        run;

	/* Get expression info for gene */

        data temp_exp;
          set test_exp;
          if gene_id="&inGene.";
        run;

        proc sort data=temp_snp;
             by gene_id subject_id;
             run;
      

        /* Merge in all expression data for a gene for each individual.
         Since each individual only has ONE genotype, this should be a 1:Many merge! */

        data temp_eqtl_data;
           merge temp_snp (in=in1) temp_exp (in=in2);
           by gene_id subject_id ;
           if in1 and in2 then output;
        run;

        /* Stack dataset */
        
        data eqtl_data_stacked;
           set eqtl_data_stacked temp_eqtl_data;
	run;

        /*increment the iterator of the loop*/
        %let iter=%eval(&iter.+1);
    %end;
%mend;

/*call the macro*/
%loopOverDatasets();


/* Make permenant */

data eqtl.all_data;
   set eqtl_data_stacked;
run;


