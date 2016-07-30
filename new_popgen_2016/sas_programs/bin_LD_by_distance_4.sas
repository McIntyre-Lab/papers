* Written by R;
*  write.foreign(ld4_out, "/home/fnew/dsim/ld/ld4_dist.txt", "/home/fnew/ufgi_share/SHARE/McIntyre_Lab/ethanol/Sim_Pop_Gen/sas_programs/bin_LD_by_distance_4.sas",  ;

DATA  rdata ;
INFILE  "/home/fnew/dsim/ld/ld4_dist.txt" 
     DSD 
     LRECL= 22 ;
INPUT
 CHR_A
 R2
 dist
;
RUN;
