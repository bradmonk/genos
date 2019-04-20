function [TRX,TRL,TEX,TEL,TRi,VAi,TEi,W] = invalidmx(PVTRX, LPTR, PVTEX, LPTE, RATIO, NHOLD)
%--------------------------------------------------------------------------
% PVTRX(: , 1)  =  PHE.SRR;        % COL1: ID
% PVTRX(: , 2)  =  PHE.AD .*2 -1;  % COL2: AD
% PVTRX(: , 3)  =  PHE.COHORTNUM;  % COL3: COHORT
% PVTRX(: , 4)  =  PHE.APOE;       % COL4: APOE
% PVTRX(: , 5)  =  PHE.WEIGHTS;    % COL5: WEIGHTS
% PVTRX(: , 6)  =  PHE.AGE;        % COL6: AGE
% PVTRX(: , 7)  =  PHE.BRAAK;      % COL7: BRAAK
% PVTRX(: , 8)  =  PHE.BRAGEZ;     % COL8: BRAGEZ
% PVTRX(: , 9)  =  PHE.AGEZ;       % COL9: AGEZ
%--------------------------------------------------------------------------


% CONCATINATE PV MATRIX WITH QUADRADIC INTERACTION TERM
% TRAINX = [VPTR VPTR.^2]';
% TESTX = [VPTE VPTE.^2]';

TRAINX  =  PVTRX;
TESTX   =  PVTEX;

TRAINL  =  LPTR';
TESTL   =  LPTE';



%--- ADD 1/2 OF TEST SET TO TRAINING MATRIX FOR VALIDATION
nTe = size(TESTX,1);
r  = round(nTe .* RATIO);
VALIX  = TESTX(1:r,:);
TEVX   = TESTX(r+1:end,:);
VALIL  = TESTL(1:r,:);
TEVL   = TESTL(r+1:end,:);


szTR = size(TRAINX,1);
szVA = size(VALIX,1)-NHOLD;
szTE = NHOLD;
sTR = 1;
eTR = szTR;
sVA = eTR+1;
eVA = sVA + szVA - 1;
sTE = eVA+1;
eTE = sTE + szTE - 1;
iTR = sTR:eTR;
iVA = sVA:eVA;
iTE = sTE:eTE;

TRVX = [TRAINX; VALIX];

TRVL = [TRAINL; VALIL];

[TRi,VAi,TEi] = divideind(size(TRVX,1),iTR,iVA,iTE);


TRX  =  TRVX;
TRL  =  TRVL;

TEX  =  TEVX;
TEL  =  TEVL;


disp(' ');
fprintf('size(TRX) : %5.f x %0.f \n', size(TRX))
fprintf('size(TEX) : %5.f x %0.f \n', size(TEX))
fprintf('size(TRL) : %5.f x %0.f \n', size(TRL))
fprintf('size(TEL) : %5.f x %0.f \n', size(TEL))

disp(' ');
fprintf('TRAINING N : %6.f  \n', numel(iTR))
fprintf('VALIDATE N : %6.f  \n', numel(iVA))
fprintf('PRETEST  N : %6.f  \n', numel(iTE))
fprintf('HOLDOUT  N : %6.f  \n', length(TEL))


end