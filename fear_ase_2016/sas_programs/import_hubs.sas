PROC IMPORT OUT= WORK.hubs 
            DATAFILE= "/home/agerken/McLab/alison_g/cegs_ase_explore/Poo
l_2015_gene_lists/Table7_gene_list.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;
