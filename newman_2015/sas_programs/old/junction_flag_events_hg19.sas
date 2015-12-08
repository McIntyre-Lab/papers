/* Add libraries */

libname splice '/home/jrbnewman/McLab/junction_annotations/sas_data/';
libname splice2 '/home/jrbnewman/Desktop/sas_temp/';

/* import junctions */

data junctions_annotated_all;
   set splice2.junctions_annotated_all_hg19;
run;


/* Merge in exon reference data */

proc sort data=junctions_annotated_all;
   by exonA;
run;

proc sort data=splice.donor_exons_hg19;
   by exonA;
run;

data junction_donor_ref no_junct oops;
   merge junctions_annotated_all (in=in1) splice.donor_exons_hg19 (in=in2);
   by exonA;
   if in1 and in2 then output junction_donor_ref;
   else if in2 then output no_junct; *489520 obs, this is ok;
   else output oops; *0 obs;
   rename geneA_id=gene_id;
run;


proc sort data=junction_donor_ref;
   by exonB;
run;

proc sort data=splice.acceptor_exons_hg19;
   by exonB;
run;

data junction_donor_accptr_refs no_junct oops;
   merge junction_donor_ref (in=in1) splice.acceptor_exons_hg19 (in=in2);
   by exonB;
   if in1 and in2 then output junction_donor_accptr_refs;
   else if in2 then output no_junct; *489262 obs, this is ok;
   else output oops; *0 obs;
   drop geneB_id; *don't need duplicated gene_ids;
run;



/* Add in entries for previous donor exon */
/* We want to check if previous donor exon is in the same grouping as the current one */
/* Sort by chr, gene, start, stop */

proc sort data=junction_donor_accptr_refs;
   by gene_id chr exonA_ref_stop exonA_ref_start exonA_ref; * we want to sort by the same donor site and keep exons together;
run;

data junction_list;
    set junction_donor_accptr_refs;
    exonA_prev=lag1(exonA); *get previous exonA_id;
    exonA_start_prev=lag1(exonA_start); *get previous exonA_id;
    exonA_stop_prev=lag1(exonA_stop); *get previous exonA_id;
    exonA_ref_prev=lag1(exonA_ref); *get previous exonA_id;
    exonA_ref_start_prev=lag1(exonA_ref_start); *get previous exonA_id;
    exonA_ref_stop_prev=lag1(exonA_ref_stop); *get previous exonA_id;
run;



/* Add in entries for previous acceptor exon */
/* We want to check if previous acceptor exon is in the same grouping as the current one */
/* Sort by chr, gene, start, stop */
proc sort data=junction_list;
   by chr gene_id exonB_ref_stop exonB_ref_start exonB_ref; * we want to sort by the same acceptor site and keep exons together;
run;

data junction_list2;
    set junction_list;
    exonB_prev=lag1(exonB); *get previous exonB_id;
    exonB_start_prev=lag1(exonB_start); *get previous exonB_id;
    exonB_stop_prev=lag1(exonB_stop); *get previous exonB_id;
    exonB_ref_prev=lag1(exonB_ref); *get previous exonB_id;
    exonB_ref_start_prev=lag1(exonB_ref_start); *get previous exonB_id;
    exonB_ref_stop_prev=lag1(exonB_ref_stop); *get previous exonB_id;
run;

/* Sort for genes */
proc sort data=junction_list2;
   by gene_id chr exonB_start exonB_stop exonA_start exonA_stop;
run;

/* Flagging exon skipping junctions */
data set_flag_exonBskip;
    set junction_list2;
    by gene_id;
    retain flag_exonBskip;
    prev_exonflag=lag(flag_exonBskip);
    if first.gene_id then flag_exonBskip=0; *if first junction for gene, then no skipped exon;
    else if exonA_ref = exonA_ref_prev then do; *check if exonA is the same;
        if exonB_ref=exonB_ref_prev then do; *check if same ref;
            if prev_exonflag=1 then flag_exonBskip=1;
            end;
        else if exonB_ref ne exonB_ref_prev then flag_exonBskip=1; *check if exonB is different;
        else flag_exonBskip=0;
        end;
    else do; *if exonA2 is not identical to exonA1 do this;
        if exonA_start ge exonA_ref_start then do; *check if exonA2 overlaps exonA1;
            if exonA_stop le exonA_ref_stop then do; *if this is not met, then is not an overlapping exon;
                if exonB_start gt exonB_stop_prev then flag_exonBskip=1;  *check if exonB is different;
                else flag_exonBskip=0;
                end;
            else flag_exonBskip=0;
            end;
        else flag_exonBskip=0; *otherwise exonB is the same;
        end;
    drop prev_exonflag;
run;

/* Check to see if junctions are identical */
data set_flag_identical;
    set set_flag_exonBskip;
    by gene_id;
    if first.gene_id then flag_identical=1; *if first junction for gene, then no skipped exon;
    else if exonA = exonA_prev and exonB = exonB_prev then flag_identical=1; *check if same exons;
    else if exonA_stop=exonA_stop_prev then do; *check if same donor;
        if exonB_start=exonB_start_prev then flag_identical=1; *check if same acceptor;
        else flag_identical=0; *otherwise not identical;
        end;
    else flag_identical=0; *otherwise not identical;
run;


/* Check to see if alternative acceptor */
/* Scenario 1 - first exon of gene - altdonor=0, altacceptor=0 */
/* Scenario 2 - donor and acceptor sites are the same as reference - altdonor=0, altacceptor=0 */
/* Scenario 2a - same donor, exon skip - altdonor=0, altacceptor=0, exonskip=1 */
/* Scenario 3 - same donor, alternative acceptor to ref  - altdonor=0, altacceptor=1 */
/* Scenario 3a - same donor, exon skip, alternative acceptor to ref  - altdonor=0, altacceptor=1, exonskip=1 */
/* Scenario 4 - alternative donor site, same acceptor to ref - altdonor=1, altacceptor=0 */
/* Scenario 4a - alt donor, exon skip - altdonor=1, altacceptor=0, exonskip=1 */
/* Scenario 5 - alt donor, alt acceptor to ref  - altdonor=1, altacceptor=1 */
/* Scenario 5a - alt donor, exon skip, alt acceptor - altdonor=1, altacceptor=1, exonskip=1 */ 

data set_flag_acceptor_donor;
    set set_flag_identical;
    by gene_id;
/* Scenario 1 - first exon of gene */
    if first.gene_id then flag_alt_donor=0; *first junction for gene is not alt donor;
    if first.gene_id then flag_alt_acceptor=0; *first junction for gene is not alt acceptor;
/* Scenario 2 and 2a - donor and acceptor sites same as ref */
    if exonA_stop = exonA_ref_stop then do; *only look at events with shared donors;
       flag_alt_donor=0;
       if exonB_start = exonB_ref_start then flag_alt_acceptor=0;
/* Scenario 3 and 3a - same donor, alt acceptor */
       else flag_alt_acceptor=1;
       end;
/* Scenario 4 and 4a - alt donor, same acceptor */
    else if exonA_stop ne exonA_ref_stop then do; *only look at events with different donors;
        flag_alt_donor=1;
        if exonB_start = exonB_ref_start then flag_alt_acceptor=0;
/* Scenario 5 and 5a - same donor, alt acceptor */
        else flag_alt_acceptor=1;
        end;
run;


/* add a "fusion junction id" and alt fusion-donor/acceptor

data junction_fusion_cat;
    length fusion_cat $20.;
    set set_flag_acceptor_donor;
    fusion_cat=catx("|",fusionA_id,fusionB_id);
    if flag_fusiondonor=1 then flag_alt_fusiondonor=0; *invert fusiondonor;
    else flag_alt_fusiondonor=1;
    if flag_fusionacceptor=1 then flag_alt_fusionacceptor=0; *invert fusionacceptor;
    else flag_alt_fusionacceptor=1;
run;  */
    
    
/* Make dataset permenant */

/* Only keeping flags donor/acceptor info */
data splice2.junctions_w_flags_hg19;
    set set_flag_acceptor_donor;
    keep
        chr
        gene_id
        xscripts_cat
        genes_cat
        junction_id
        exonA
        exonB
        donor_start
        donor_stop
        acceptor_start
        acceptor_stop
        strand
        num_genes
        num_xscripts
        flag_multigene
        flag_junction_annotated
        flag_firstexon
        flag_altfirstexon
        flag_lastexon
        flag_altlastexon
        flag_exonBskip
        flag_alt_donor
        flag_alt_acceptor
    ;
run;

/* Clean up */
    proc datasets nolist;
        delete junctions_annotated_all junction_donor_ref no_junct oops
        junction_donor_accptr_refs no_junct oops junction_list junction_list2
        set_flag_exonBskip set_flag_identical set_flag_acceptor_donor
        junction_fusion_cat
        ;
        run;
        quit;

