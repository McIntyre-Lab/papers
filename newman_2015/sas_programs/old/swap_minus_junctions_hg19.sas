/* swapping acceptor and donor information for minus-strand junctions */

libname splice '/home/jrbnewman/McLab/junction_annotations/sas_data/';
libname splice2 '/media/jrbnewman/SAS_WRK1/';

/* split on strand */

data junctions_w_flags_plus junctions_w_flags_minus oops;
   set splice2.junctions_w_flags_hg19; *2610258 junctions;
   if strand="+" then output junctions_w_flags_plus; *2573138 junctions;
   else if strand="-" then output junctions_w_flags_minus; *37120 junctions;
   else output oops; *0 obs, woo!;
run;

/* swapping flags on minus */
/* going to keep the exonA/B, donor/acceptor coords the same for now */
/* just in case we need to extract these to rebuild coordinates - will probably get messy otherwise! */

data junctions_minus_flagswap;
   set junctions_w_flags_minus;
   flag_firstexon_swap=flag_lastexon;
   flag_altfirstexon_swap=flag_altlastexon;
   flag_lastexon_swap=flag_firstexon;
   flag_altlastexon_swap=flag_altfirstexon;
   flag_alt_donor_swap=flag_alt_acceptor;
   flag_alt_acceptor_swap=flag_alt_donor;
run;

data junctions_minus_flagswap2;
   set junctions_minus_flagswap;
   flag_firstexon=flag_firstexon_swap;
   flag_altfirstexon=flag_altfirstexon_swap;
   flag_lastexon=flag_lastexon_swap;
   flag_altlastexon=flag_altlastexon_swap;
   flag_alt_donor=flag_alt_donor_swap;
   flag_alt_acceptor=flag_alt_acceptor_swap;

   drop
      flag_firstexon_swap
      flag_altfirstexon_swap
      flag_lastexon_swap
      flag_altlastexon_swap
      flag_alt_donor_swap
      flag_alt_acceptor_swap
     ;
run;

/* merge back with plus junctions */

data junctions_w_flags_fixed;
    set junctions_w_flags_plus junctions_minus_flagswap2;
run;

/* make permenant */

data splice2.junctions_w_flags_fixed_hg19;
   set junctions_w_flags_fixed;
run;

/* Clean up */
    proc datasets nolist;
        delete junctions_w_flags_plus junctions_w_flags_minus oops
        junctions_minus_flagswap junctions_minus_flagswap2
        junctions_w_flags_fixed
        ;
        run;
        quit;

