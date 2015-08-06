     data WORK.SIM    ;
     infile '/home/jfear/mclab/cegs_sem_sd_paper/simulation/test_cfa.csv' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
        informat sample $4. ;
        informat y1 best32. ;
        informat y2 best32. ;
        informat y3 best32. ;
        informat y4 best32. ;
        informat y5 best32. ;
        informat y6 best32. ;
        format sample $4. ;
        format y1 best12. ;
        format y2 best12. ;
        format y3 best12. ;
        format y4 best12. ;
        format y5 best12. ;
        format y6 best12. ;
     input
                 sample $
                 y1
                 y2
                 y3
                 y4
                 y5
                 y6
     ;
     run;


proc calis data=WORK.SIM method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000;
    path
        f1 -> y1 y2 y3,
        f2 -> y4 y5 y6
    ;
    run;
