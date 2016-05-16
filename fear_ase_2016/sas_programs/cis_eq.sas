/*******************************************************************************
* Filename: cis_eq.sas
*
* Author: Justin M Fear | jfear@ufl.edu
*
* Description: Sergey and Lauren developed a set of equtions found here:
* Nuzhdin, S. V, Friesen, M. L., & McIntyre, L. M. (2012).
* Genotype-phenotype mapping in a post-GWAS world. Trends in Genetics : TIG,
* 28(9), 421–6. doi:10.1016/j.tig.2012.06.003
* 
* which potentially allow for the identification of cis- and trans-effects.
* Here I try using these questions and test if they give reasonable results.
* 
* Basics: For a given gene the expression level of Eii of allele i in F1
* genotype i.
* 
* Eii=μ+Ci+(Ti+Tt)/2
* 
* Eti=μ+Ct+(Ti+Tt)/2
* 
* For each allele the cis- and trans-effects are deviations from the
* population means, we expect that they will sum to zero:
* 
* ∑ni=1Ci=0
* 
* ∑ni=1Ti=0
* 
* Then the expected difference in expression between the Line and Tester
* allele over the entire population is:
* 
* ∑ni=1Eti−Eiin
* 
* Which can be re-written as
* 
* ∑ni=1Ct−Cin=Ct
* 
* The cis-effect of allele i can be estimated by:
* 
* ˆCi=ˆEii−ˆEti+ˆCt
*
*******************************************************************************/

/* Libraries
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
libname DMEL551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
libname GENELIST '!MCLAB/useful_dmel_data/gene_lists/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Run this data set if you want to use proportion adjusted values */
    * sum_line / (sum_line + sum_tester) * 1000
    * sum_tester / (sum_line + sum_tester) * 1000
    ;

    /*
    data maren;
        set CEGS.data_for_cis_eq;
        Eti = tester_prop_adj;
        Eii = line_prop_adj;
        diffEtiEii = Eti - Eii;
        diffEiiEti = Eii - Eti;
        sumEtiEii = Eti + Eii;
        meanEtiEii = sumEtiEii / 2;
        diff_n = diffEtiEii / geno_count;
        Eti_n = Eti / geno_count;
        run;
    */

/* Run this data set if you want to use mean center values */
    * mu = mean( (sum_line + sum_tester) / 2 )
    * sum_line - mu
    * sum_tester - mu
    ;

    data maren;
        set CEGS.data_for_cis_eq;
        Eti = tester_mean_center;
        Eii = line_mean_center;
        diffEtiEii = Eti - Eii;
        diffEiiEti = Eii - Eti;
        sumEtiEii = Eti + Eii;
        meanEtiEii = sumEtiEii / 2;
        diff_n = diffEtiEii / geno_count;
        Eti_n = Eti / geno_count;
        run;

/* Calculate Cis-tester */
    * Ct = sum( Eti - Eii ) / n
    ;

    proc sort data=maren;
        by mating_status fusion_id;
        run;

    proc means data= maren noprint;
        by mating_status fusion_id;
        output out=sums sum(diff_n)=cis_tester sum(Eti_n)=sum_Eti_n mean(meanEtiEii)=mu;
        run;

    data ct;
        merge maren (in=in1) sums (in=in2);
        by mating_status fusion_id;
        if mu < 0.000000001 then mu = 0; * if mu is really small, say from mean centering then set to 0;
        drop _type_ _freq_;
        run;

/* Calculate additional effects */
    * Cis-Line: Ci = Eii - Eti + Ct
    * Trans-Tester: Tt = 2( sum(Eti)/n - mu - Ct)
    * Trans-Line: Ti = 2(Eti - mu - Ct) - Tt
    ;
    data effects;
        set ct;
        cis_line = diffEiiEti + cis_tester;
        trans_tester = 2 * (sum_Eti_n - mu - cis_tester);
        trans_line = 2 * (Eti - mu - cis_tester) - trans_tester;
        run;

/* Create Permanent Data Set */
    data CEGS.cis_eq;
        set effects;
        keep mating_status fusion_id line q5_mean_theta flag_AI_combined mean_apn
        cis_line cis_tester trans_line trans_tester mu;
        run;

/* Export Data Set */
    proc export data=CEGS.cis_eq outfile='!MCLAB/cegs_ase_paper/pipeline_output/cis_effects/cis_eq_full_output_w_ai_calls.csv' dbms=csv replace;
        putnames=yes;
        run;
