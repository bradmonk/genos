%% GENOS: 
% 
%--------------------------------------------------------------------------
% 
% SUMMARY TABLE OF THE 24 COHORTS
% 
% COHID    CONSOR    STUDY        COHORT    CASES    CTRLS    TOTAL    %CASE    EQL  BRAAK ID GOOD
% 01       DGC       Adult_Chng    ACT        323      945     1268       25    323    1   01    1
% 02       ADGC      AD_Centers    ADC       2438      817     3255       75    817    0   02    1
% 03       CHARGE    Athrosclro    ASC         39       18       57       68     18    0   03    0
% 04       CHARGE    Aus_Stroke    SKE        121        5      126       96      5    0   04    0
% 05       ADGC      ChiT_Aging    CHA         27      204      231       12     27    0   05    0
% 06       CHARGE    CardioHlth    CHS        250      583      833       30    250    1   06    1
% 07       ADGC      Hispanic_S    HSP        160      171      331       48    160    0   07    0
% 08       CHARGE    Erasmus_Er    ERF         45        0       45      100      0    0   08    0
% 09       CHARGE    Framingham    FHS        157      424      581       27    157    1   09    1
% 10       ADGC      Gene_Diffs    GDF        111       96      207       54     96    1   10    1
% 11       ADGC      NIA_LOAD      LOD        367      109      476       77    109    1   11    1
% 12       ADGC      Aging_Proj    MAP        138      277      415       33    138    1   12    1
% 13       ADGC      Mayo_Clini    MAY        250       99      349       72     99    1   13    1
% 14       ADGC      Miami_Univ    MIA        186       14      200       93     14    1   14    0
% 15       ADGC      AD_Genetic    MIR        316       15      331       95     15    0   15    0
% 16       ADGC      Mayo_cl_PD    MPD          0       20       20        0      0    1   16    0
% 17       ADGC      NationC_AD    NCA        160        0      160      100      0    1   17    0
% 18       ADGC      Wash_Unive    RAS         46        0       46      100      0    1   18    0
% 19       ADGC      Relig_Ordr    ROS        154      197      351       44    154    1   19    1
% 20       CHARGE    RotterdamS    RDS        276      813     1089       25    276    0   20    1
% 21       ADGC      Texas_AD_S    TAR        132       12      144       92     12    0   21    0
% 22       ADGC      Un_Toronto    TOR          9        0        9      100      0    0   22    0
% 23       ADGC      Vanderbilt    VAN        210       26      236       89     26    1   23    1
% 24       ADGC      WashNY_Age    WCA         34      116      150       23     34    0   24    1
% 
% 
% GOODCOHORTS = [1 2         6 7   9 10 11 12 13               19 20     23 24]
% BRAKCOHORTS = [1           6     9 10 11 12 13 14   16 17 18 19        23   ]
%--------------------------------------------------------------------------
% 
% THE ADSP DATASET - WHAT'S BEING IMPORTED?
%
% 
% The dataset that will be loaded in STEP-1 below will import 5 variables
% and store them into a structural array named 'ADSP'. If you type ADSP
% into the command prompt you will see...
% 
% 
% >> ADSP
% 
% ADSP = 
%   struct with fields:
% 
%     PHEN: [10910×22 table]
%     LOCI: [94483×24 table]
%     CASE: {94483×1 cell}
%     CTRL: {94483×1 cell}
%     USNP: {94483×1 cell}
% 
% 
% ...(maybe not in this specific order) the 5 container variables.
% 
% 
% PHEN    a table containing the phenotype information for each
%         participant. If you type head(ADSP.PHEN) you can see what
%         data each column contains.
% 
% 
% LOCI    a table containing genotype info for each exome variant locus.
%         if you type head(ADSP.LOCI) you can see what data each column
%         contains.
% 
% 
% 
% CASE    the last three are cell arrays, containing a list of 
% CTRL    participant IDs & HET/HOM status. Each have 1 cell per row
% USNP    of LOCI (~94483 cells); they are (at least upon import) 
%         pre-sorted in corresponding order, which allows us to
%         iterate over each loci/cell and tally each HET (+1) or 
%         HOM (+2) that matches subsets of participant IDs. (there
%         is an optimized function specifically designed to
%         perform this task, as you will see below).
%         
% 
%--------------------------------------------------------------------------





%==========================================================================
%% STEP-1: LOAD THE DATASET
%==========================================================================
%{.
close all; clear; clc; rng('shuffle');
genosdir = fileparts(which('GENOS.m'));
matlabdir = fileparts(which('MATLABS.m'));
cd(genosdir);


subfuncpath = [genosdir filesep 'genosfunctions'];
datasetpath = [matlabdir filesep 'genosdata'];
gpath = [genosdir pathsep subfuncpath pathsep datasetpath];
addpath(gpath)


which('GENOSDATA.mat')
ADSP = load('GENOSDATA.mat');



disp('dataset loaded')
clearvars -except FILES MATDAT IJ ADSP

%}


%==========================================================================
%%   CARBON COPY MAIN VARIABLES FROM ADSP.STRUCT
%==========================================================================
%
% After evaluating this section, each variable will be copied from the
% ADSP structural array to their own base variable. This is done so
% that (1) you can access their data directly (e.g. LOCI.GENE(1:5) instead
% of ADSP.LOCI.GENE(1:5) ) and so that (2) you can always restart fresh
% here, by running this segment of code, rather than having to import the
% data from the .mat file in the section above.

LOCI = ADSP.LOCI;
CASE = ADSP.CASE;
CTRL = ADSP.CTRL;
USNP = ADSP.USNP;
PHEN = ADSP.PHEN;


clc; clearvars -except FILES MATDAT IJ ADSP PHEN LOCI CASE CTRL USNP
head(PHEN)
head(LOCI)



%==========================================================================
%==========================================================================
%==========================================================================
%%
%
% GET DATASET PATHS FOR MACHINE LEARNING
%
%==========================================================================
%==========================================================================
%==========================================================================
clearvars -except FILES MATDAT IJ ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL


MATDIR = 'F:\ML\genosdat\APOE3x\APOE33_SUBSAM50_FISHP';
% MATDIR = 'F:\ML\genosdat\APOE3x\APOE33e34_SUBSAM50_FISHP';
% MATDIR = 'F:\ML\genosdat\APOE2x4x\APOE2x4x_SUBSAM50_FISHP';
% MATDIR = 'F:\ML\genosdat\APOE2x3x4x\APOE2x3x4x_SUBSAM50_FISHP';
FILES.w = what(MATDIR);
FILES.w.mat



LOOPDATA.HILO_TRMEAN = zeros(200,10);
LOOPDATA.HILO_TRHIMU = zeros(200,10);
LOOPDATA.HILO_TRPOPU = zeros(200,10);
LOOPDATA.HILO_HOMEAN = zeros(200,10);
LOOPDATA.HILO_HOHIMU = zeros(200,10);
LOOPDATA.HILO_HOPOPU = zeros(200,10);
LOOPDATA.LOHI_TRMEAN = zeros(200,10);
LOOPDATA.LOHI_TRHIMU = zeros(200,10);
LOOPDATA.LOHI_TRPOPU = zeros(200,10);
LOOPDATA.LOHI_HOMEAN = zeros(200,10);
LOOPDATA.LOHI_HOHIMU = zeros(200,10);
LOOPDATA.LOHI_HOPOPU = zeros(200,10);


%==========================================================================
%% RUN MAIN LOOP
%==========================================================================
for IJ = 1:10
clearvars -except LOOPDATA FILES MATDAT ADSP LOCI PHEN CASE CTRL USNP IJ



MATDAT = load([FILES.w.path filesep FILES.w.mat{IJ}]);



VLOCI     = MATDAT.LOCI;
VCASE     = CASE;
VCTRL     = CTRL;
VUSNP     = USNP;
VTRCASE   = MATDAT.TRCASE;
VTRCTRL   = MATDAT.TRCTRL;
VTECASE   = MATDAT.TECASE;
VTECTRL   = MATDAT.TECTRL;


% SET MAIN FISHP TO TRAINING GROUP FISHP
VLOCI.FISHP      = VLOCI.TRFISHP;
VLOCI.FISHOR     = VLOCI.TRFISHOR;
VLOCI.CASEREF    = VLOCI.TRCASEREF;
VLOCI.CASEALT    = VLOCI.TRCASEALT;
VLOCI.CTRLREF    = VLOCI.TRCTRLREF;
VLOCI.CTRLALT    = VLOCI.TRCTRLALT;


% SORT VARIANTS BY EITHER TRFISHP|CHRPOS
[~,j]  = sort(VLOCI.TRFISHP);
VLOCI  = VLOCI(j,:);
VCASE  = VCASE(j);
VCTRL  = VCTRL(j);
VUSNP  = VUSNP(j);


%==========================================================================
%%                               HILO
%==========================================================================
hiloTRMEAN = zeros(200,1);
hiloTRHIMU = zeros(200,1);
hiloTRPOPU = zeros(200,1);
hiloHOMEAN = zeros(200,1);
hiloHOHIMU = zeros(200,1);
hiloHOPOPU = zeros(200,1);
%==========================================================================
for vi = 1:200

% REMOVE VARS FROM HIGH-TO-LOW P-VALUE
VI = fliplr(1:200);
HI2LOW = 1;
SNPn = VI(vi);
SNPi = 1:SNPn;


% EXTRACT TOP-N NUMBER OF VARIANTS
XLOCI  = VLOCI(SNPi,:);
XCASE  = VCASE(SNPi);
XCTRL  = VCTRL(SNPi);
XUSNP  = VUSNP(SNPi);


%==========================================================================
%      MAKE  RECTANGLE  NEURAL NET  VARIANT MATRIX
%==========================================================================

TRPHE = [VTRCASE; VTRCTRL];
TEPHE = [VTECASE; VTECTRL];



% SCRAMBLE TRAINING PHENOTYPE ORDER
NVARS  = size(TRPHE,1);         % Total number of people
k      = randperm(NVARS)';      % Get N random ints in range 1:N
TRPHE  = TRPHE(k,:);            % Scramble Phenotype table

% SCRAMBLE TESTING PHENOTYPE ORDER
NVARS  = size(TEPHE,1);         % Total number of people
k      = randperm(NVARS)';      % Get N random ints in range 1:N
TEPHE  = TEPHE(k,:);            % Scramble Phenotype table



% MAKE THE NEURAL NET TRAINING & TESTING MATRICES
VTRX = mlmx_mex(XCASE,XCTRL,XUSNP,...
    TRPHE.SRR,TRPHE.AD,TRPHE.COHORTNUM,TRPHE.AGE,TRPHE.APOE,TRPHE.BRAAK);

VTEX = mlmx_mex(XCASE,XCTRL,XUSNP,...
    TEPHE.SRR,TEPHE.AD,TEPHE.COHORTNUM,TEPHE.AGE,TEPHE.APOE,TEPHE.BRAAK);


% [VTRX, TRX, TRL] = makeomicsnet(VLOCI,VCASE,VCTRL,VUSNP,TRPHE,[-1 -0 1 3]);
% [VTEX, TEX, TEL] = makeomicsnet(VLOCI,VCASE,VCTRL,VUSNP,TEPHE,[-1 -0 1 3]);

%==========================================================================
%                  LOGISTIC REGRESSION 
%==========================================================================


TXX = VTRX(:,7:end);
HXX = VTEX(:,7:end);
TX = [ones(size(TXX,1),1) TXX]; % ADD AN INTERCEPT COLUMN
HX = [ones(size(HXX,1),1) HXX]; % ADD AN INTERCEPT COLUMN


TL = VTRX(:,2);
HL = VTEX(:,2);




% PERFORM THE SO-CALLED MACHINE LEARNING STEP
BETA = normaleq_mex(TX,TL);
%BETA = pinv(TX' * TX) * (TX' * TL);

fprintf('\n Solved GLM OLS for %0.f beta coefficients. \n\n',size(BETA,1));  




% GET TRAINED LINEAR MODEL PREDICTIONS (ACTIVATIONS)
TRfx = nansum( TX .* BETA' ,2);
HOfx = nansum( HX .* BETA' ,2);


% DESCRETIZE PREDICTIONS
TRy = round(TRfx);
HOy = round(HOfx);


% GRADE THE PREDICTIONS
TRmu = nanmean(TRy == TL);
HOmu = nanmean(HOy == HL);




% [TRAINED] HIGH CONFIDENCE PREDICTIONS
TRhi = (TRfx>.8) | (TRfx<.2);
TRhiN  = TRy(TRhi) == TL(TRhi);
TRhiMu = nanmean(TRhiN);
TRhiPop = nanmean(TRhi);


% [HOLDOUT] HIGH CONFIDENCE PREDICTIONS
HOhi = (HOfx>.8) | (HOfx<.2);
HOhiN  = HOy(HOhi) == HL(HOhi);
HOhiMu = nanmean(HOhiN);
HOhiPop = nanmean(HOhi);




hiloTRMEAN(vi) = TRmu;
hiloTRHIMU(vi) = TRhiMu;
hiloTRPOPU(vi) = TRhiPop;
hiloHOMEAN(vi) = HOmu;
hiloHOHIMU(vi) = HOhiMu;
hiloHOPOPU(vi) = HOhiPop;


clc;
disp('%=================================================================');
disp('IJ | vi | min(SNPi) | max(SNPi) | numel(BETA):'); 
disp([IJ vi min(SNPi) max(SNPi) numel(BETA)]);
disp('----------  TRAINING SET  ----------')
fprintf('TRAIN correct %0.0f%% overall(pop:100%%)\n' ,(TRmu .* 100))
fprintf('TRAIN correct %0.0f%% hicon  (pop: %0.0f%%)\n\n',TRhiMu.*100,TRhiPop.*100)

disp('----------  TESTING SET  ----------')
fprintf('TEST correct %0.0f%% overall(pop:100%%)\n' ,(HOmu .* 100))
fprintf('TEST correct %0.0f%% hicon  (pop: %0.0f%%)\n\n',HOhiMu.*100,HOhiPop.*100)
disp('%=================================================================');



end
%===========================================

% FILL IN NAN VALUES WIHT NANMEAN
hiloTRMEAN(isnan(hiloTRMEAN)) = nanmean(hiloTRMEAN);
hiloTRHIMU(isnan(hiloTRHIMU)) = nanmean(hiloTRHIMU);
hiloTRPOPU(isnan(hiloTRPOPU)) = nanmean(hiloTRPOPU);
hiloHOMEAN(isnan(hiloHOMEAN)) = nanmean(hiloHOMEAN);
hiloHOHIMU(isnan(hiloHOHIMU)) = nanmean(hiloHOHIMU);
hiloHOPOPU(isnan(hiloHOPOPU)) = nanmean(hiloHOPOPU);

%==========================================================================
%%                               LOHI
%==========================================================================
lohiTRMEAN = zeros(200,1);
lohiTRHIMU = zeros(200,1);
lohiTRPOPU = zeros(200,1);
lohiHOMEAN = zeros(200,1);
lohiHOHIMU = zeros(200,1);
lohiHOPOPU = zeros(200,1);
%==========================================================================
for vi = 1:200

% REMOVE VARS FROM LOW-TO-HIGH P-VALUE
VI = fliplr(1:200);
HI2LOW = 0;
SNPn = vi;
SNPi = SNPn:200;

% EXTRACT TOP-N NUMBER OF VARIANTS
XLOCI  = VLOCI(SNPi,:);
XCASE  = VCASE(SNPi);
XCTRL  = VCTRL(SNPi);
XUSNP  = VUSNP(SNPi);


%==========================================================================
%      MAKE  RECTANGLE  NEURAL NET  VARIANT MATRIX
%==========================================================================
TRPHE = [VTRCASE; VTRCTRL];
TEPHE = [VTECASE; VTECTRL];

% SCRAMBLE TRAINING PHENOTYPE ORDER
NVARS  = size(TRPHE,1);         % Total number of people
k      = randperm(NVARS)';      % Get N random ints in range 1:N
TRPHE  = TRPHE(k,:);            % Scramble Phenotype table

% SCRAMBLE TESTING PHENOTYPE ORDER
NVARS  = size(TEPHE,1);         % Total number of people
k      = randperm(NVARS)';      % Get N random ints in range 1:N
TEPHE  = TEPHE(k,:);            % Scramble Phenotype table



% MAKE THE NEURAL NET TRAINING & TESTING MATRICES
% [VTRX, TRX, TRL] = makeomicsnet(VLOCI,VCASE,VCTRL,VUSNP,TRPHE,[-1 -0 1 3]);
% [VTEX, TEX, TEL] = makeomicsnet(VLOCI,VCASE,VCTRL,VUSNP,TEPHE,[-1 -0 1 3]);

VTRX = mlmx_mex(XCASE,XCTRL,XUSNP,...
    TRPHE.SRR,TRPHE.AD,TRPHE.COHORTNUM,TRPHE.AGE,TRPHE.APOE,TRPHE.BRAAK);

VTEX = mlmx_mex(XCASE,XCTRL,XUSNP,...
    TEPHE.SRR,TEPHE.AD,TEPHE.COHORTNUM,TEPHE.AGE,TEPHE.APOE,TEPHE.BRAAK);



%==========================================================================
%                  LOGISTIC REGRESSION 
%==========================================================================


TXX = VTRX(:,7:end);
HXX = VTEX(:,7:end);
TX = [ones(size(TXX,1),1) TXX]; % ADD AN INTERCEPT COLUMN
HX = [ones(size(HXX,1),1) HXX]; % ADD AN INTERCEPT COLUMN


TL = VTRX(:,2);
HL = VTEX(:,2);




% PERFORM THE SO-CALLED MACHINE LEARNING STEP
% BETA = pinv(TX' * TX) * (TX' * TL);
BETA = normaleq_mex(TX,TL);
fprintf('\n Solved GLM OLS for %0.f beta coefficients. \n\n',size(BETA,1));  




% GET TRAINED LINEAR MODEL PREDICTIONS (ACTIVATIONS)
TRfx = nansum( TX .* BETA' ,2);
HOfx = nansum( HX .* BETA' ,2);


% DESCRETIZE PREDICTIONS
TRy = round(TRfx);
HOy = round(HOfx);


% GRADE THE PREDICTIONS
TRmu = nanmean(TRy == TL);
HOmu = nanmean(HOy == HL);




% [TRAINED] HIGH CONFIDENCE PREDICTIONS
TRhi = (TRfx>.8) | (TRfx<.2);
TRhiN  = TRy(TRhi) == TL(TRhi);
TRhiMu = nanmean(TRhiN);
TRhiPop = nanmean(TRhi);


% [HOLDOUT] HIGH CONFIDENCE PREDICTIONS
HOhi = (HOfx>.8) | (HOfx<.2);
HOhiN  = HOy(HOhi) == HL(HOhi);
HOhiMu = nanmean(HOhiN);
HOhiPop = nanmean(HOhi);




lohiTRMEAN(vi) = TRmu;
lohiTRHIMU(vi) = TRhiMu;
lohiTRPOPU(vi) = TRhiPop;
lohiHOMEAN(vi) = HOmu;
lohiHOHIMU(vi) = HOhiMu;
lohiHOPOPU(vi) = HOhiPop;


clc;
disp('%=================================================================');
disp('IJ | vi | min(SNPi) | max(SNPi) | numel(BETA):'); 
disp([IJ vi min(SNPi) max(SNPi) numel(BETA)]);
disp('----------  TRAINING SET  ----------')
fprintf('TRAIN correct %0.0f%% overall(pop:100%%)\n' ,(TRmu .* 100))
fprintf('TRAIN correct %0.0f%% hicon  (pop: %0.0f%%)\n\n',TRhiMu.*100,TRhiPop.*100)

disp('----------  TESTING SET  ----------')
fprintf('TEST correct %0.0f%% overall(pop:100%%)\n' ,(HOmu .* 100))
fprintf('TEST correct %0.0f%% hicon  (pop: %0.0f%%)\n\n',HOhiMu.*100,HOhiPop.*100)
disp('%=================================================================');


end


% FILL IN NAN VALUES WIHT NANMEAN
lohiTRMEAN(isnan(lohiTRMEAN)) = nanmean(lohiTRMEAN);
lohiTRHIMU(isnan(lohiTRHIMU)) = nanmean(lohiTRHIMU);
lohiTRPOPU(isnan(lohiTRPOPU)) = nanmean(lohiTRPOPU);
lohiHOMEAN(isnan(lohiHOMEAN)) = nanmean(lohiHOMEAN);
lohiHOHIMU(isnan(lohiHOHIMU)) = nanmean(lohiHOHIMU);
lohiHOPOPU(isnan(lohiHOPOPU)) = nanmean(lohiHOPOPU);


%==========================================================================
%%                               HILO-LOHI SAVER
%==========================================================================



LOOPDATA.HILO_TRMEAN(1:200,IJ) = hiloTRMEAN;
LOOPDATA.HILO_TRHIMU(1:200,IJ) = hiloTRHIMU;
LOOPDATA.HILO_TRPOPU(1:200,IJ) = hiloTRPOPU;
LOOPDATA.HILO_HOMEAN(1:200,IJ) = hiloHOMEAN;
LOOPDATA.HILO_HOHIMU(1:200,IJ) = hiloHOHIMU;
LOOPDATA.HILO_HOPOPU(1:200,IJ) = hiloHOPOPU;

LOOPDATA.LOHI_TRMEAN(1:200,IJ) = lohiTRMEAN;
LOOPDATA.LOHI_TRHIMU(1:200,IJ) = lohiTRHIMU;
LOOPDATA.LOHI_TRPOPU(1:200,IJ) = lohiTRPOPU;
LOOPDATA.LOHI_HOMEAN(1:200,IJ) = lohiHOMEAN;
LOOPDATA.LOHI_HOHIMU(1:200,IJ) = lohiHOHIMU;
LOOPDATA.LOHI_HOPOPU(1:200,IJ) = lohiHOPOPU;

end
%% SAVE LOOP DATA

%------------------------------------------%
pause(1)
dt=char(datetime(datetime,'Format','yyyy-MM-dd-HH-mm-ss'));
save(['F:\ML\genosdat\APOE33_P01P50_ML.mat'],'LOOPDATA','FILES');
disp('File Saved.')
pause(1)
%------------------------------------------%
