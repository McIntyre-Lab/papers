/********************************************************************************
* Create and export dataset for gaussian graphical models using the sex
* determination subset. Then run GGM R script. GGM allows you to select
* edges by an FDR cutoff or by specifiying a number of edges. Two
* different FDR cutoffs were used (0.2, 0.8) along with outputing the
* top 20 edges. GGM was run on the entire DSRP dataset.
********************************************************************************/

data stack;
    set SEM.dsrp_stack_rank;
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
            if suffix ne '' then sym = trim(symbol) || '_' || trim(suffix);
            else sym = symbol;
        end;
        else sym = cgnumber;
        drop num suffix symbol cgnumber;
        run;

data stack_w_sex2;
    set stack_w_sex;
    rename sym = symbol;
    run;

proc sort data=stack_w_sex2;
    by patRIL matRIL;
    run;

proc transpose data=stack_w_sex2 out=cgflip;
    by patRIL matRIL;
    var expRank;
    id symbol;
    run;

data w_sample;
    retain sample;
    format sample $ 20.;
    set cgflip;
    sample = strip(patRIL) || '_' || strip(matRIL);
    drop patRIL matRIL _name_ _label_;
    run;

* I could not get this to export all of the columns (11065). So I ended up
* exporting the file to a csv from JMP genomics. When saving under options select csv.;

/*
proc export data=w_sample outfile='/tmp/data_for_ggm.csv' dbms=csv replace;
    putnames=yes;
    run;

* Run GGM On Mated and Virgin;
data _null_;
    call system('Rscript $MCLAB/cegs_sem_sd_paper/r_programs/dsrp_ggm_gene.R');
    run;
*/

/* Clean up */
proc datasets nolist;
    delete cgflip;
    delete stack;
    delete stack_w_sex;
    delete stack_w_sex2;
    delete w_sample;
    run; quit;
