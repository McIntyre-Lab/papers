
/* Iterative merging of SNP data and expression data */

libname eqtl '/scratch/lfs/patcon/jnewman/sas_analysis/eqtls';

/* Create a temporary empty dataset with proper formats */
data temp_fus;
   set eqtl.expression_data;
   if _n_ = 1;
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

/* Making some working copies of the expression and genotype data to avoid any potential conflicts */        


data expression_data_&sysparm.;
   set eqtl.expression_data;
run;

data genotype_data_&sysparm.;
   set eqtl.genotype_data;
run;

data snp2gene_index_&sysparm.;
   set eqtl.snp2gene_index;
   if _n_ gt (&sysparm. - 1)*500 and _n_ le (&sysparm. + 0)*500; *only takes in 20 observations!;
run;

        proc sort data=expression_data_&sysparm.;
             by gene_id subject_id;
             run;

/* Iterate through SNP-Gene index to merge SNP data and expression data */

%macro loopOverDatasets();
    /* declare macro variables */
    %local iter inSNP inGene NObs;


    /* set number of iterations */
    data _null_;
    set snp2gene_index_&sysparm.;
    call symput("NObs", _N_);
    run;
 
    %put NObs;

    /*initiate loop*/
    %let iter=1;
    %do %while (&iter.<= &NObs.);

       /* Sort snp-gene index of iteration and all genotypes, then merge */

        /*get libref and dataset name for dataset you will work on during this iteration*/
        data _NULL_;
            set snp2gene_index_&sysparm. (firstobs=&iter. obs=&iter.); *only read 1 record;
            *write the Name and Event to the macro variables;
            call symput("inSNP",strip(snp_id));
            call symput("inGene",strip(gene_id));
            *NOTE: strip function is to remove any trailing blanks. This should not be a problem, but it is here just in case;
        run;

        /* get SNP-gene pair */
        data temp_snp; *assuming you want to apply the changes to the dataset itself;
            set genotype_data_&sysparm.;
            if snp_id = "&inSNP." and gene_id="&inGene." then output;
        run;

	/* Get expression info for gene */

        data temp_exp;
          set expression_data_&sysparm.;
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

data eqtl.all_exp_data_w_snps_&sysparm.;
   set eqtl_data_stacked;
run;



