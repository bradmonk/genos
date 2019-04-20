function [PBTRX, PBTEX] = bragematrix(PVTRX, PVTEX)
%--------------------------------------------------------------------------
% PVTRX(: , 1)  =  PHE.SRR;        % COL1: ID
% PVTRX(: , 2)  =  PHE.AD .*2 -1;  % COL2: AD
% PVTRX(: , 3)  =  PHE.COHORTNUM;  % COL3: COHORT
% PVTRX(: , 4)  =  PHE.APOE;       % COL4: APOE
% PVTRX(: , 5)  =  PHE.WEIGHT;     % COL5: WEIGHT
% PVTRX(: , 6)  =  PHE.AGE;        % COL6: AGE
% PVTRX(: , 7)  =  PHE.BRAAK;      % COL7: BRAAK
% PVTRX(: , 8)  =  PHE.BRAGEZ;     % COL8: BRAGEZ
% PVTRX(: , 9)  =  PHE.AGEZ;       % COL9: AGEZ
%--------------------------------------------------------------------------
%  AGEZ & BRAGEZ EXPLAINED
%--------------------------------------------------------------------------
%{

When training a classifier to recognize contributors to Alzheimers,
we want the classifier to know that young CASES and old CONTROLS
should be weighted more heavily than old CASES and young CONTROLS.

Furthermore, we want the classifier to know that CASES with high BRAAK
scores and CONTROLS with low BRAAK scores should be weighted more heavily.

Since weights need to be a single value per person, we compute a weight
called 'BRAGE' that is based on both Age and BRAAK scores, such that...

    BRAGE = AGE - BRAK

...the worse their 'BRAK' value the more their age is reduced (the first
'BRAK' value below is 5.4):

    BRAGE(BRAAK==0)   : AGE + 5.4
    BRAGE(BRAAK==1)   : AGE + 3.4
    BRAGE(BRAAK==2)   : AGE + 1.4
    BRAGE(BRAAK==NaN) : AGE - 0.0
    BRAGE(BRAAK==3)   : AGE - 2.5
    BRAGE(BRAAK==4)   : AGE - 5.5
    BRAGE(BRAAK==5)   : AGE - 7.5
    BRAGE(BRAAK==6)   : AGE - 9.5

So for example if someone is 80 years old and has a BRAAK==0 their
BRAGE equals 85.4; if they had a BRAAK==6 their BRAGE equals 70.5

From there a Z-score is computed for both AGE & BRAGE, but for CASES
the Z +/- sign is flipped...

    CASE_AGEZ   = zscore(CASE_AGE)   * -1
    CASE_BRAGEZ = zscore(CASE_BRAGE) * -1

    CTRL_AGEZ   = zscore(CTRL_AGE)   
    CTRL_BRAGEZ = zscore(CTRL_BRAGE)

the sign for CASES is flipped because the youngest cases with the highest
BRAAK scores will have the lowest BRAGE, but we want those CASES to have
the greatest weight, so the sign is flipped.

%}
%--------------------------------------------------------------------------


% SCALE TRAINING GROUP
WEIGHT = rescale(PVTRX(:,5), -1,1);
V = PVTRX(:,10:end);
for i = 1:size(V,1)
    x = V(i,:);
    x(x>0) = x(x>0) + WEIGHT(i);
    V(i,:) = x;
end
PBTRX = [PVTRX(:,1:9) V];


% SCALE TEST GROUP
WEIGHT = rescale(PVTEX(:,5), -1,1);
V = PVTEX(:,10:end);
for i = 1:size(V,1)
    x = V(i,:);
    x(x>0) = x(x>0) + WEIGHT(i);
    V(i,:) = x;
end
PBTEX = [PVTEX(:,1:9) V];


% TRAX = [V VTRX(:,2)];
% TRAX = [VTRPXV(:,2:end) VTRL(:,1)];


end