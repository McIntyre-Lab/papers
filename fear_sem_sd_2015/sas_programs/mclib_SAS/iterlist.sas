*----------------------------| MACRO ITERLIST() |----------------------------*;
* Iterate items of a list through a fragment of SAS code                     *;
*----------------------------------------------------------------------------*;
* CODE = Fragment of SAS code that is to be repeated, using                  *;
* the question mark (?) as a token to be replaced by items in                *;
* the space-delimited list embodied by &LIST                                 *;
* LIST = List of items                                                       *;
*Copied from --                                                              *;
* http://www.wuss.org/proceedings08/08WUSS%20Proceedings/papers/cod/cod06.pdf*;
*-RMG 05/18/2012                                                             *;
*----------------------------------------------------------------------------*;
%macro iterlist
(
code =
,list =
)
;
%*** ASSIGN EACH ITEM IN THE LIST TO AN INDEXED MACRO VARIABLE &&ITEM&I ;
%let i = 1;
%do %while (%cmpres(%scan(&list., &i.)) ne );
%let item&i. = %cmpres(%scan(&list., &i.));
%let i = %eval((&i. + 1);
%end;
%*** STORE THE COUNT OF THE NUMBER OF ITEMS IN A MACRO VARIABLE: &CNTITEM;
%let cntitem = %eval((&i. - 1);
%*** EXPRESS CODE, REPLACING TOKENS WITH ELEMENTS OF THE LIST, IN SEQUENCE;
%do i = 1 %to &cntitem.;
%let codeprp = %qsysfunc(tranwrd(&code.,?,%nrstr(&&item&i..)));
%unquote(&codeprp.)
%end;
%mend iterlist;
