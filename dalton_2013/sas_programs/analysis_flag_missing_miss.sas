/* I re-ran the analysis without the dsx and did several models
 * For the fusion level model I need to identify fusions that have
 * missing values so they can be separated.
 */

data flag;
    set fru.resid_by_fusion_miss; * 1801390 obs, with 47405 uniq fusions;
    if logrpkm = '' then flag_missing = 1;
    else flag_missing = 0;
    run;

proc means data=flag noprint;
    by fusion_id;
    output out=flag_sum sum(flag_missing)=sum_flags;
    run;

data fru.flag_missing_miss;
    set flag_sum;
    if sum_flags > 0 then flag_missing = 1;
    else flag_missing = 0;
    keep fusion_id flag_missing;
    run; * 47405 obs;

