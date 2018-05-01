/********************************************************************************
* This macro is given a design file as a data set and another macro call using
* variable names from that dataset.
*
* To call macro use the following syntax:
*   %iterdataset(dataset=<DATASET>, function=%nrstr(<MACRONAME with vars>));
*
* Real world example
*   %iterdataset(dataset=SEM.design_by_bio_rep, function=%nrstr(%test(&vars)));
*
********************************************************************************/
%macro iterdataset(dataset=,function=);
    %local dsid now total rows cols rc;
    %let dsid = %sysfunc(open(&dataset));
    %let now = 0;
    %let rows = %sysfunc(attrn(&dsid, nobs));
    %let cols = %sysfunc(attrn(&dsid, nvars));

    %do %while(%sysfunc(fetch(&dsid)) = 0); %* outer loop across rows;
        %let now = %eval(&now + 1);

        %do i = 1 %to &cols; %* inner loop across coloumns;
            %local v t;
            %let v=%sysfunc(varname(&dsid,&i));
            %local &v;
            %let t = %sysfunc(vartype(&dsid,&i));
            %let &v = %sysfunc(getvar&t(&dsid,&i));
        %end;

        %unquote(&function);

    %end;
    %let rc = %sysfunc(close(&dsid));
%mend;
