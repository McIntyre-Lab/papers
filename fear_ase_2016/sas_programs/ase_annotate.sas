
libname DMEL551 's:\McIntyre_lab\useful_dmel_data\flybase551\sasdata';
libname cegs_a "s:\McIntyre_lab\alison_g\cegs_ase_explore\sas_data";
libname species "S:\McIntyre_Lab\berlin-c167-hyb-solexa\Analysis RMG\SAS Data\Analysis_Step_6_ASE_model";
libname rita2 "S:\McIntyre_Lab\berlin-c167-hyb-solexa\Analysis RMG\SAS Data\Bayesian";
libname begun "S:\McIntyre_Lab\Gene lists for comparisons\Begun et al PLoS 2007";
data WORK.SEX_DET;
        infile 's:\McIntyre_lab\useful_dmel_data\gene_lists\pathways\sex_det.csv' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
        informat symbol $43. ;
        informat primary_fbgn $11. ;
        format symbol $43. ;
        format primary_fbgn $11. ;
        input
            symbol $
            primary_fbgn $
            ;
        run;

    proc sort data=sex_det;
        by primary_fbgn;
        run;

    proc sort data=DMEL551.FB551_si_fusion_2_gene_id;
        by FBgn;
        run;

    data fus_sxlist oops;
        merge  DMEL551.FB551_si_fusion_2_gene_id (in=in2 rename=(FBgn=primary_fbgn)) sex_det (in=in1);
        by primary_fbgn;
        if in1 and not in2 then output oops;
        else if in1 and in2 then output fus_sxlist;
        keep fusion_id;
        run; * oops had none!!;

    proc sort data=DMEL551.FB551_si_fusion_2_gene_id;
        by fusion_id;
        run;
		
    proc sort data=fus_sxlist;
        by fusion_id;
        run;


		data si_fus_plussx;
		merge DMEL551.FB551_si_fusion_2_gene_id(in=in1) fus_sxlist(in=in2);
		by fusion_id;
		if in2 then flag_sexdet=1;
			else flag_sexdet=0;
			run;



proc contents data=species.All_results_bygene2;
*7515;
run;

proc contents data=DMEL551.fusions2go;
run;

proc contents data= cegs_a.rita_results_w_flags;
run;
 *7515 obs;
proc contents data= cegs_a.pool_2015_gene_lists;
run;

proc contents data=rita2.collapsed_bayes_to_til_wmktests;
run;
