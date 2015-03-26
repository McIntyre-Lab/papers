/********************************************************************************
********************************************************************************/

* libname fru '!MCLAB/arbeitman_fru_network/sasdata';
* libname dmel530 '!MCLAB/useful_dmel_data/flybase530/sasdata';

/* Multiple Motif Enrichment Analysis */

%macro enrichment_analysis(letter);

    %macro response_type(type);

        proc freq data=fru_&letter._motif ;
            tables flag_&letter._&type.*flag_multi_motif /agree cmh chisq ;
            run;

    %mend response_type;
    %response_type(ind);
    %response_type(rep);

%mend enrichment_analysis;
%enrichment_analysis(a);
%enrichment_analysis(b);
%enrichment_analysis(c);

