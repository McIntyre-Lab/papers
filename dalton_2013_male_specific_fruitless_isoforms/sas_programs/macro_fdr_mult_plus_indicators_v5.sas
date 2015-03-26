%macro fdr_mult(data_in,pvalue,id,level1,level2,level3);

    data raw_p;
        set &data_in;
        rename  &pvalue=raw_p;
        *keep &id &pvalue;
        run;

    proc sort data=raw_p;
        by raw_p;
        run;

    proc multtest pdata=raw_P fdr out=pout1 noprint;
        run; quit;

    data fdr_&pvalue;
        set pout1;
        length fdr_level_&pvalue $8;
        rename fdr_p=fdr_&pvalue;
        if fdr_p =. and raw_p>.99 then do;
            fdr_level_&pvalue="0_tan";
            fdr_p=1;
            end;
        else if fdr_p =. then fdr_level_&pvalue="";
        else if fdr_p le &level3 then fdr_level_&pvalue="3_red";
        else if fdr_p le &level2 then fdr_level_&pvalue="2_orange";
        else if fdr_p le &level1 then fdr_level_&pvalue="1_yellow";
        else if fdr_p > &level1 then fdr_level_&pvalue="0_tan";

        if fdr_level_&pvalue="3_red" then fdr_&pvalue._red=1;
            else if fdr_level_&pvalue="" then fdr_&pvalue._red="";
            else fdr_&pvalue._red=0;

        if fdr_level_&pvalue="3_red" or fdr_level_&pvalue="2_orange"  then fdr_&pvalue._orange=1;
            else if fdr_level_&pvalue="" then fdr_&pvalue._orange="";
            else fdr_&pvalue._orange=0;
            drop raw_p;

        if fdr_level_&pvalue="3_red" or fdr_level_&pvalue="2_orange" or fdr_level_&pvalue="1_yellow" then fdr_&pvalue._yellow=1;
            else if fdr_level_&pvalue="" then fdr_&pvalue._yellow="";
            else fdr_&pvalue._yellow=0;
            drop raw_p;

    run;

    proc sort data=fdr_&pvalue;
        by &id;
        run;

%mend;
