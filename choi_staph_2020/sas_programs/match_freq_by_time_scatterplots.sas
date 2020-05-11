libname relapse "!MCLAB/staph/seoung_ho_project/sasdata";
libname TB14 "/home/ammorse/TB14/staph_relapse/sasdata";

ods listing ;
ods _all_ close ;
ods select all ;

data snp_counts  ;
set relapse.relapse_snp_counts  ;
log_match_freq = log(match_freq) ;
run ;


proc univariate data = snp_counts ;
var match_freq interval_in_days;
run ;  /*  0.35 to 1.0, 15 to 1108 */

ods pdf file = "!MCLAB/staph/seoung_ho_project/output/match_freq_by_interval_scatter_plots.pdf" ;

title "match_freq by interval" ;
proc sgplot data = snp_counts ;
where is_a_pair = "YES";
scatter x = interval_in_days y = match_freq / markerattrs=(symbol=CircleFilled size=10px) ;
run;

title "no_match_freq by interval" ;
proc sgplot data = snp_counts ;
where is_a_pair = "YES";
scatter x = interval_in_days y = noMatch_freq / markerattrs=(symbol=CircleFilled size=10px) ;
run;

ods pdf close ;


/* playing with plotting below..... */
data snp_counts_pairs  ;
set relapse.relapse_snp_counts  ;
where is_a_pair = "YES";
run ;

proc sgplot data = snp_counts_pairs ;
vline interval_in_days/response = match_freq ;
vline interval_in_days/response = noMatch_freq ;
run ;

plots=matrix(histogram) ;
var match_freq interval_in_days ;
run;


plots=matrix(




data pair_cnts ;
set snp_counts ;
where is_a_pair = "YES" ;
run ;

proc univariate data = snp_counts noprint;
var match_freq;
class is_a_pair ;
histogram match_freq / endpoints=(0 to 1500 by 50);
label is_a_pair = "Where isolate pair ";
run ;




title "" ;
proc univariate data = snp_counts noprint;
var sum_no_match;
class is_a_pair ;
histogram sum_no_match / endpoints=(0 to 1500 by 50);
label is_a_pair = "Where isolate pair ";
run ;

proc univariate data = snp_counts noprint;
var noMatch_freq;
class is_a_pair ;
histogram noMatch_freq / endpoints=(0 to 1 by .1);
label is_a_pair = "Where isolate pair ";
run ;


proc sgplot data = snp_counts ;
scatter x=



