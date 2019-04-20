function [PVMX, PVMXi, PVL] = featurematrix(ADMX,CASE,CTRL,USNP,PHE,varargin)



% USE CUSTOM OR DEFAULT NUMERICAL VARIANT LABELS
%------------------------------------------------
if nargin > 5
    REFUNKALT = varargin{1};
else
              %refref unkunk refalt altalt
    REFUNKALT = [-1     0     1     3];
end






% CREATE VARIANT FEATURE MATRIX
%------------------------------------------------
CACOID = PHE.SRR;

vMX  = zeros( size(CACOID,1) , size(ADMX,1) );


for nn = 1:size(CASE,1)

    CASES = CASE{nn};
    CTRLS = CTRL{nn};
    CACO = [CASES; CTRLS];
    UNCAS = USNP{nn};


    if any(UNCAS)
        [~,Uj] = ismember(CACOID , UNCAS );
        vMX(Uj>0,nn) = 1;
    end
    if any(CACO)
        CACOSRR = CACO(:,1);
        CACOHH  = CACO(:,2);

        [~,Aj] = ismember(CACOID , CACOSRR );
        [~,Uj] = ismember(CACOID , UNCAS );

        Af = Aj(Aj>0);

        vMX(Aj>0,nn) = (CACOHH(Af)+1); %HETALT:2  HOMALT:3
        vMX(Uj>0,nn) = 1;              %UNCALL:1  HOMREF:0

    end

end





% UPDATE NUMERICAL VARIANT LABELS
%------------------------------------------------
refref = vMX==0;
unkunk = vMX==1;
refalt = vMX==2;
altalt = vMX==3;

PVX = vMX;
PVX(refref) =  REFUNKALT(1);     % HOMREF: -1
PVX(unkunk) =  REFUNKALT(2);     % UNCALL:  0
PVX(refalt) =  REFUNKALT(3);     % HETALT:  2
PVX(altalt) =  REFUNKALT(4);     % HOMALT:  5



% ADD PHENOTYPE INFO TO VARIANT MATRIX
%------------------------------------------------
PVMX = padarray(PVX,[0 9],0,'pre');
PVMX(: , 1)  =  PHE.SRR;        % COL1: ID
PVMX(: , 2)  =  PHE.AD .*2 -1;  % COL2: AD
PVMX(: , 3)  =  PHE.COHORTNUM;  % COL3: COHORT
PVMX(: , 4)  =  PHE.APOE;       % COL5: APOE
PVMX(: , 5)  =  PHE.SEX;        % COL4: SEX
PVMX(: , 6)  =  PHE.AGE;        % COL5: AGE
PVMX(: , 7)  =  PHE.BRAAK;      % COL5: BRAAK
PVMX(: , 8)  =  PHE.BRAGEZ;     % COL6: BRAGEZ
PVMX(: , 9)  =  PHE.AGEZ;       % COL4: AGEZ



% CREATE LABEL MATRIX
%------------------------------------------------
PVL = zeros( 2 , numel(PHE.AD) );  % Make 2-D class label matrix
PVL(1,:) = PHE.AD==1;              % In col-1 set CASEs index 1's
PVL(2,:) = PHE.AD==0;              % In col-2 set CTRLs index 1's












% DISPLAY MATRIX RANK
%------------------------------------------------

fprintf('\n\n---------------- MATRIX RANK --------------------\n')
fprintf('MATRIX ROWS: %4.f  (PEOPLE)\n', size(PVX,1));
fprintf('MATRIX COLS: %4.f  (VARIANTS)\n', size(PVX,2));
fprintf('MATRIX RANK: %4.f  (FULL-RANK WHEN COLS==RANK)\n', rank(PVX));




% DETERMINE ALT ALLELE COVERAGE
%------------------------------------------------
PPL_ALTs = sum(sum( PVMX(:,10:end) >= REFUNKALT(3), 2)>0);
VAR_ALTs = sum( PVMX(:,10:end) >= REFUNKALT(3), 1);

fprintf('\n\n----------- MATRIX ALT-ALLELE STATS -------------\n')
fprintf('ALTs Per Locus  MIN: %6.f \n',  min(VAR_ALTs));
fprintf('ALTs Per Locus  MAX: %6.f \n',  max(VAR_ALTs));
fprintf('ALTs Per Locus MEAN: %6.1f \n', mean(VAR_ALTs));
fprintf('\nPEOPLE WITH ALTS: %4.f of %4.f \n\n\n', PPL_ALTs , size(PVMX,1));


% DISPLAY REFREF UNKUNK REFALT ALTALT HISTOGRAMS
%------------------------------------------------
% PP_REFREF = sum( PVMX(:,10:end) == REFUNKALT(1), 2);
% PP_UNKUNK = sum( PVMX(:,10:end) == REFUNKALT(2), 2);
% PP_REFALT = sum( PVMX(:,10:end) == REFUNKALT(3), 2);
% PP_ALTALT = sum( PVMX(:,10:end) == REFUNKALT(4), 2);

% fh01 = figure('Units','normalized','OuterPosition',[.01 .05 .95 .90],...
%               'Color','w','MenuBar','none');
% subplot(2,2,1); histogram(PP_REFREF); title('REF/REFs   Per Person')
% subplot(2,2,2); histogram(PP_UNKUNK); title('UNK/UNKs   Per Person')
% subplot(2,2,3); histogram(PP_REFALT); title('REF/ALTs   Per Person')
% subplot(2,2,4); histogram(PP_ALTALT); title('ALT/ALTs   Per Person')



% DISPLAY FEATURE MATRIX PREVIEW
%------------------------------------------------
% T = PVMX(1:9,1:13);
% T = array2table(T);
% T.Properties.VariableNames = {'ID','AD','COH','APOE','SEX','AGE',...
%                               'BRAAK','BRAGEZ','AGEZ',...
%                              ['V1_' char(ADMX.GENE(1))],...
%                              ['V2_' char(ADMX.GENE(2))],...
%                              ['V3_' char(ADMX.GENE(3))],...
%                              ['V4_' char(ADMX.GENE(4))]};
% T.ID = int64(T.ID); disp(T)



% ADD AN INTERCEPT COLUMN IF TRAINING NEURAL NETS
%------------------------------------------------
PVMXi = [PVMX(:,1:9) ones(size(PVMX,1),1) PVMX(:,10:end)];


end