%macro unix_chk_dir(dir=) ; 
   %local rc fileref ; 
   %let rc = %sysfunc(filename(fileref,&dir)) ; 
   %if %sysfunc(fexist(&fileref))  %then 
      %put NOTE: The directory "&dir" exists ; 
   %else 
     %do ; 
         %sysexec(mkdir -p &dir) ; 
         %put %sysfunc(sysmsg()) The directory has been created. ; 
   %end ; 
   %let rc=%sysfunc(filename(fileref)) ; 
%mend chk_dir ; 
