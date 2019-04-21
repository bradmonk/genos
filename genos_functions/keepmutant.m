function [KEEP] = keepmutant(ADMX,L)
%{

GENOMIX_filtereffect takes the following inputs

ADMX  : Variant table
L     : 1x10 logical array


This function will return a Vx1 logical array 
(where V is the number of rows in ADMX)
that indicates which rows to include and which
to remove based on the variant EFFECT tag in ADMX.

The order of these values is shown below, with an example
logical array...

    [SYN MIS STG STL SPA SPD INT NCE UTR CUK]
L = [ 0   1   1   1   1   1   0   1   1   0]


In this example, variants with the effect label: SYN INT CUK 
will be assigned zeros in the logical array. This array can
be used to then filter the dataset...

    ADMX    = ADMX( INCLUDE , : );
    SNPCASE = SNPCASE(INCLUDE);
    SNPCTRL = SNPCTRL(INCLUDE);
    ADMX.VID = ( 1:size(ADMX,1) )';

%}



% find rows with the following variant types:
SYN = (strcmp(cellstr(ADMX.EFFECT),'SYN')); % SYN: synonymous
MIS = (strcmp(cellstr(ADMX.EFFECT),'MIS')); % MIS: missense
STG = (strcmp(cellstr(ADMX.EFFECT),'STG')); % STG: stop-gained
STL = (strcmp(cellstr(ADMX.EFFECT),'STL')); % STL: stop-lost
SPA = (strcmp(cellstr(ADMX.EFFECT),'SPA')); % SPA: splice-acceptor
SPD = (strcmp(cellstr(ADMX.EFFECT),'SPD')); % SPD: splice-donor
INT = (strcmp(cellstr(ADMX.EFFECT),'INT')); % INT: intron
NCE = (strcmp(cellstr(ADMX.EFFECT),'NCE')); % NCE: non-coding-exon
UTR = (strcmp(cellstr(ADMX.EFFECT),'UTR')); % UTR: 3-/5-prime UTR
CUK = (strcmp(cellstr(ADMX.EFFECT),'CUK')); % CUK: coding-unknown


VSUM = sum(sum([MIS SYN INT NCE STG UTR SPD SPA STL CUK]));
fprintf('missense        (MIS): % 8.0f \n', sum(MIS) )
fprintf('synonymous      (SYN): % 8.0f \n', sum(SYN) )
fprintf('intron          (INT): % 8.0f \n', sum(INT) )
fprintf('non-coding-exon (NCE): % 8.0f \n', sum(NCE) )
fprintf('stop-gained     (STG): % 8.0f \n', sum(STG) )
fprintf('3-/5-prime UTR  (UTR): % 8.0f \n', sum(UTR) )
fprintf('splice-donor    (SPD): % 8.0f \n', sum(SPD) )
fprintf('splice-acceptor (SPA): % 8.0f \n', sum(SPA) )
fprintf('coding-unknown  (CUK): % 8.0f \n', sum(CUK) )
fprintf('stop-lost       (STL): % 8.0f \n', sum(STL) )
fprintf('-------------------------------- \n')
fprintf('total                  % 8.0f \n',   VSUM   )


A = [];

if L(1);   A=[A SYN];  end;
if L(2);   A=[A MIS];  end;
if L(3);   A=[A STG];  end;
if L(4);   A=[A STL];  end;
if L(5);   A=[A SPA];  end;
if L(6);   A=[A SPD];  end;
if L(7);   A=[A INT];  end;
if L(8);   A=[A NCE];  end;
if L(9);   A=[A UTR];  end;
if L(10);  A=[A CUK];  end;

KEEP = sum(A,2) > 0;


end

% IN _MAIN DATASET
% ------------------------------
% missense        (MIS):  811382
% synonymous      (SYN):  477751
% intron          (INT):  101808
% non-coding-exon (NCE):   45477
% stop-gained     (STG):   21778
% 3-/5-prime UTR  (UTR):   12885
% splice-donor    (SPD):    6422
% splice-acceptor (SPA):    4830
% stop-lost       (STL):     791
% coding-unknown  (CUK):    1847*removed
% ------------------------------
% total                1,483,123

% IN _EQUAL DATASET
% ------------------------------
% missense        (MIS):   740908 
% synonymous      (SYN):   438945 
% intron          (INT):        0 
% non-coding-exon (NCE):        0 
% stop-gained     (STG):    19690 
% 3-/5-prime UTR  (UTR):        0 
% splice-donor    (SPD):     5807 
% splice-acceptor (SPA):     4380 
% ------------------------------
% total


% IN _RAND DATASET
% ------------------------------
% missense        (MIS):  811382
% synonymous      (SYN):  477752
% intron          (INT):  101807
% non-coding-exon (NCE):   45476
% stop-gained     (STG):   21778
% 3-/5-prime UTR  (UTR):   12885
% splice-donor    (SPD):    6422
% splice-acceptor (SPA):    4830
% coding-unknown  (CUK):    1847
% stop-lost       (STL):     791
% ------------------------------
% total                1,483,123
