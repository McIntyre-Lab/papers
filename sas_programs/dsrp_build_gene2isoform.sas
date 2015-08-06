/********************************************************************************
* I need to create a dataset which is by gene, but then lists isoforms.
* I will use this gene list to iterate through the genome and feed my
* python script to generate models.
********************************************************************************/

/*
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

proc sort data=SEM.dsrp_sbs_combine_sym;
    by patRIL matRIL;
    run;

proc transpose data=SEM.dsrp_sbs_combine_sym out=flip;
    by patRIL matRIL;
    var :;
    run;

data flip2;
    set flip;
    keep _name_;
    label _name_ = ' ';
    run;

data iso2gene;
    length symbol $12. ;
    set flip2;
    count = index(_name_, '_');
    if count > 0 then do;
        symbol = trim(substr(_name_, 1, count-1));
    end;
    else symbol = _name_;
    pos = prxmatch("/^fl_2_d/", _name_);
    if pos > 0 then symbol = 'fl_2_d';
    rename _name_ = isoform;
    drop count pos;
    run;

proc sort data=iso2gene nodupkey;
    by symbol isoform ;
    run;

proc transpose data=iso2gene out=isoflip;
    by symbol;
    var isoform;
    run;
    
data combined;
    length iso_cat $300.;
    set isoflip;
    iso_cat = catx(" ", OF col1-col8);
    run;

data gene2isoform;
    retain symbol iso_cat;
    set combined;
    keep symbol iso_cat;
    run;

data SEM.dsrp_gene2isoform;
    set gene2isoform;
    run;

proc export data=SEM.dsrp_gene2isoform outfile='!MCLAB/cegs_sem_sd_paper/design_file/dsrp_gene2isoform.csv' dbms=csv replace;
    putnames=no;
    run;
