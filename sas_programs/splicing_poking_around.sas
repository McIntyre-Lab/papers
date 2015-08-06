
libname sem 'S:\McIntyre_Lab\cegs_sem_sd_paper\sas_data';

ods html;
title '607 vs 609';
proc glimmix data=sem.cegsV_splice_data plots=studentpanel;
where fusion_id="F57609_SI" or fusion_id="S57607_SI";
    by symbol_cat;
	where symbol_cat="InR";
    class fusion_id line rep;
    model uq_log_uq_center = line|fusion_id /htype=1;
    lsmeans line*fusion_id /slice=line slice=fusion_id;
    *output out=resid_by_symbol r=resid p=pred;
    *ods output Tests1=model_ts FitStatistics=model_fs;
    run;


	
ods html;
title '607 vs 609';
proc glimmix data=sem.cegsV_splice_data plots=studentpanel;
where fusion_id="F57609_SI" or fusion_id="S57607_SI";
    by symbol_cat;
	where symbol_cat="InR";
    class fusion_id line rep;
    model uq_log_uq_center = line|fusion_id /htype=1;
    lsmeans line*fusion_id /slice=line slice=fusion_id;
    *output out=resid_by_symbol r=resid p=pred;
    *ods output Tests1=model_ts FitStatistics=model_fs;
    run;
title 'by fusion_is';
proc sort data=sem.cegsV_splice_data;
by symbol_cat fusion_id;
	proc glm data=sem.cegsV_splice_data ;
	by symbol_cat fusion_id;
	class line;
	model uq_log_uq_center=line;
	ods output  modelanova=glmbyfusion;
	run;
quit;


ods html;
title '607 vs 601';
proc glimmix data=sem.cegsV_splice_data plots=studentpanel;
where fusion_id="F57601_SI" or fusion_id="S57607_SI";
    by symbol_cat;
	where symbol_cat="InR";
    class fusion_id line rep;
    model uq_log_uq_center = line|fusion_id /htype=1;
    lsmeans line*fusion_id /slice=line slice=fusion_id;
    *output out=resid_by_symbol r=resid p=pred;
    *ods output Tests1=model_ts FitStatistics=model_fs;
    run;


	data small_test;
	set sem.cegsV_splice_data;
	where line='r371' or line='r502';
	run;

	proc glimmix data=small_test plots=studentpanel;
where fusion_id="F57601_SI" or fusion_id="S57607_SI";
    by symbol_cat;
	where symbol_cat="InR";
    class fusion_id line rep;
    model uq_log_uq_center = line|fusion_id /htype=1;
    lsmeans line*fusion_id /slice=line slice=fusion_id;
    *output out=resid_by_symbol r=resid p=pred;
    *ods output Tests1=model_ts FitStatistics=model_fs;
    run;
