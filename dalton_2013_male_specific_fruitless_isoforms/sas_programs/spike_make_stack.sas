libname fru '!HOME/mclab/Fru_network/sasdata';

proc sort data=fru.all_coverage_spikes;
    by spike_id;
    run;

proc sort data=fru.spike_design_file;
    by spike_id;
    run;

data spike_stack;
    merge fru.all_coverage_spikes fru.spike_design_file;
    by spike_id;
    run;

proc means data = spike_stack noprint;
    by spike_id;
    var logrpkm;
    output  out=spike_avg mean=logrpkm_avg;
    run;

data spike_avg2;
    set spike_avg;
    drop _type_ _freq_;
    run;

proc sort data=spike_avg2;
    by spike_id;
    run;

proc sort data=spike_stack;
    by spike_id;
    run;

data fru.spike_stack_avg;
    merge spike_stack spike_avg2;
    by spike_id;
    run;


