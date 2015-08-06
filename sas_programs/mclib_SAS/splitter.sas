/* This macro splits a string and creates global variables with the contents of the string */

%macro splitter(string=,wordpfx=SPLITVAR,dlm=%str( ));

    %do cnt=1 %to %sysfunc(countw(&string,&dlm));
        %global &wordpfx&cnt;
        %let &wordpfx&cnt = %scan(&string,&cnt,%str(&dlm));
        %* echo macro var result to log window;
        %put &wordpfx&cnt=&&&wordpfx&cnt;
    %end;

%mend splitter;
