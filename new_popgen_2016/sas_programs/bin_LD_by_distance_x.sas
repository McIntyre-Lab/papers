* Written by R;
*  write.foreign(ldx_out, "/home/fnew/dsim/ld/ldx_dist.txt", "/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/sas_programs/bin_LD_by_distance_x.sas",  ;

DATA  rdata ;
INFILE  "/home/fnew/dsim/ld/ldx_dist.txt" 
     DSD 
     LRECL= 24 ;
INPUT
 CHR_A
 R2
 dist
;
RUN;
