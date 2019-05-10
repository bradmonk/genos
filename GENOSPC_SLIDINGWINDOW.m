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
close all; clear; clc; rng('shuffle');
P.home = fileparts(which('GENOS.m')); cd(P.home);
P.P1 = [P.home filesep 'genos_functions'];
P.P3 = [P.P1 filesep 'genos_main_functions'];
P.P2 = [P.home filesep 'genosfunctions'];
P.P4 = [P.home filesep 'genos_other'];
addpath(join(string(struct2cell(P)),pathsep,1))
cd(P.home)



ADSP = load('GENOSDATA.mat');

disp('dataset loaded')
clearvars -except P ADSP


%% SET RUN OPTIONS & PATHS FOR DATA IMPORT/EXPORT
clc; clearvars -except P ADSP


P.Nloops = 25;
P.Nvars = 500;
P.windowSize = 50;
P.Ndots = P.Nvars-P.windowSize+1;
P.Lo2Hi = 1>0; %YES


P.mainmatfile = which('GENOSDATA.mat');
P.basedir = 'F:\GENOSDATA\APOE_SUBGROUPS\APOE_22_23_24_33_34_44';
P.importdir = [P.basedir '\APOE_22_23_24_33_34_44_FISHP_V0'];
dt=char(datetime(datetime,'Format','yyyy-MM-dd-HH-mm-ss'));
P.savepathfile = [P.basedir '\SLIDEWIN50V_' dt '.mat'];

P.FILES.w = what(P.importdir);
P.Nmatfiles = numel(P.FILES.w.mat);

disp(P.FILES.w.mat);
disp(P.Nmatfiles)



if P.Nloops > P.Nmatfiles
disp('ABORTING: NOT ENOUGH MAT FILES TO RUN THAT MANY LOOPS');
return; 
end

clearvars -except P ADSP


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


clc; clearvars -except P ADSP PHEN LOCI CASE CTRL USNP
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
clearvars -except P ADSP PHEN LOCI CASE CTRL USNP

% CASES
%----------------------------------------
LOOPDATA.CATRMEAN = zeros(P.Ndots,P.Nloops);
LOOPDATA.CATRLOMU = zeros(P.Ndots,P.Nloops);
LOOPDATA.CATRHIMU = zeros(P.Ndots,P.Nloops);
LOOPDATA.CATRLOPO = zeros(P.Ndots,P.Nloops);
LOOPDATA.CATRHIPO = zeros(P.Ndots,P.Nloops);

LOOPDATA.CAHOMEAN = zeros(P.Ndots,P.Nloops);
LOOPDATA.CAHOLOMU = zeros(P.Ndots,P.Nloops);
LOOPDATA.CAHOHIMU = zeros(P.Ndots,P.Nloops);
LOOPDATA.CAHOLOPO = zeros(P.Ndots,P.Nloops);
LOOPDATA.CAHOHIPO = zeros(P.Ndots,P.Nloops);


% CTRLS
%----------------------------------------
LOOPDATA.COTRMEAN = zeros(P.Ndots,P.Nloops);
LOOPDATA.COTRLOMU = zeros(P.Ndots,P.Nloops);
LOOPDATA.COTRHIMU = zeros(P.Ndots,P.Nloops);
LOOPDATA.COTRLOPO = zeros(P.Ndots,P.Nloops);
LOOPDATA.COTRHIPO = zeros(P.Ndots,P.Nloops);

LOOPDATA.COHOMEAN = zeros(P.Ndots,P.Nloops);
LOOPDATA.COHOLOMU = zeros(P.Ndots,P.Nloops);
LOOPDATA.COHOHIMU = zeros(P.Ndots,P.Nloops);
LOOPDATA.COHOLOPO = zeros(P.Ndots,P.Nloops);
LOOPDATA.COHOHIPO = zeros(P.Ndots,P.Nloops);


%==========================================================================
%% RUN MAIN LOOP
%==========================================================================
for IJ = 1:P.Nloops
clearvars -except P ADSP PHEN LOCI CASE CTRL USNP LOOPDATA IJ



MATDAT = load([P.FILES.w.path filesep P.FILES.w.mat{IJ}]);



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


% CASES
%----------------------------------------
CATRMEAN = zeros(P.Ndots,1);
CATRLOMU = zeros(P.Ndots,1);
CATRHIMU = zeros(P.Ndots,1);
CATRLOPO = zeros(P.Ndots,1);
CATRHIPO = zeros(P.Ndots,1);

CAHOMEAN = zeros(P.Ndots,1);
CAHOLOMU = zeros(P.Ndots,1);
CAHOHIMU = zeros(P.Ndots,1);
CAHOLOPO = zeros(P.Ndots,1);
CAHOHIPO = zeros(P.Ndots,1);


% CTRLS
%----------------------------------------
COTRMEAN = zeros(P.Ndots,1);
COTRLOMU = zeros(P.Ndots,1);
COTRHIMU = zeros(P.Ndots,1);
COTRLOPO = zeros(P.Ndots,1);
COTRHIPO = zeros(P.Ndots,1);

COHOMEAN = zeros(P.Ndots,1);
COHOLOMU = zeros(P.Ndots,1);
COHOHIMU = zeros(P.Ndots,1);
COHOLOPO = zeros(P.Ndots,1);
COHOHIPO = zeros(P.Ndots,1);


%==========================================================================
for vi = 1:P.Ndots


if P.Lo2Hi
    SNPi = vi:(vi+P.windowSize-1);
else
    SNPi = (P.Nvars-P.windowSize+2-vi):(P.Nvars+1-vi);
end






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
% VTRX = mlmx_mex(XCASE,XCTRL,XUSNP,...
%     TRPHE.SRR,TRPHE.AD,TRPHE.COHORTNUM,TRPHE.AGE,TRPHE.APOE,TRPHE.BRAAK);
% 
% VTEX = mlmx_mex(XCASE,XCTRL,XUSNP,...
%     TEPHE.SRR,TEPHE.AD,TEPHE.COHORTNUM,TEPHE.AGE,TEPHE.APOE,TEPHE.BRAAK);

[VTRX, ~, ~] = mkmx(XCASE,XCTRL,XUSNP,TRPHE,[-1 -0 2 3]);
[VTEX, ~, ~] = mkmx(XCASE,XCTRL,XUSNP,TEPHE,[-1 -0 2 3]);


%==========================================================================
%                  LOGISTIC REGRESSION 
%==========================================================================

TXX = VTRX(:,10:end);
HXX = VTEX(:,10:end);
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


keyboard

% CASES
%----------------------------------------
% CASE: GET PREDICTIONS
CATRfx = TRfx(TL==1);
CAHOfx = HOfx(HL==1);

% CASE: GRADE ALL PREDICTIONS
case_tr_mu = nanmean(CATRfx > .5);
case_ho_mu = nanmean(CAHOfx > .5);

% CASE: GRADE HI CONF PREDICTIONS
QTR  = quantile(TRfx,[.2 .8]);
case_tr_lo_mu  = nanmean(CATRfx(CATRfx<QTR(1)));
case_tr_hi_mu  = nanmean(CATRfx(CATRfx>QTR(2)));
case_tr_lo_pop = nanmean(CATRfx(CATRfx<QTR(1)));
case_tr_hi_pop = nanmean(CATRfx(CATRfx>QTR(2)));
case_ho_lo_mu  = nanmean(CAHOfx(CAHOfx<QTR(1)));
case_ho_hi_mu  = nanmean(CAHOfx(CAHOfx>QTR(2)));
case_ho_lo_pop = nanmean(CAHOfx(CAHOfx<QTR(1)));
case_ho_hi_pop = nanmean(CAHOfx(CAHOfx>QTR(2)));


% CTRLS
%----------------------------------------
% CTRL: GET PREDICTIONS
COTRfx = TRfx(TL==0);
COHOfx = HOfx(HL==0);

% CTRL: GRADE ALL PREDICTIONS
ctrl_tr_mu = nanmean(COTRfx > .5);
ctrl_ho_mu = nanmean(COHOfx > .5);

% CTRL: GRADE HI CONF PREDICTIONS
QTR  = quantile(TRfx,[.2 .8]);
ctrl_tr_lo_mu  = nanmean(COTRfx(COTRfx<QTR(1)));
ctrl_tr_hi_mu  = nanmean(COTRfx(COTRfx>QTR(2)));
ctrl_tr_lo_pop = nanmean(COTRfx(COTRfx<QTR(1)));
ctrl_tr_hi_pop = nanmean(COTRfx(COTRfx>QTR(2)));
ctrl_ho_lo_mu  = nanmean(COHOfx(COHOfx<QTR(1)));
ctrl_ho_hi_mu  = nanmean(COHOfx(COHOfx>QTR(2)));
ctrl_ho_lo_pop = nanmean(COHOfx(COHOfx<QTR(1)));
ctrl_ho_hi_pop = nanmean(COHOfx(COHOfx>QTR(2)));













% [TRAINED] HIGH CONFIDENCE PREDICTIONS
QTR    = quantile(TRfx,[.2 .8]);
TRlo   = HOfx<QHO(1);
TRhi   = HOfx>QHO(2);
TRhilo = TRlo | TRhi;
TRloN  = TRy(TRlo) == TL(TRlo);
TRhiN  = TRy(TRhi) == TL(TRhi);
TRloMu  = nanmean(TRloN);
TRhiMu  = nanmean(TRhiN);
TRloPop = nanmean(TRlo);
TRhiPop = nanmean(TRhi);


% [HOLDOUT] HIGH CONFIDENCE PREDICTIONS
QHO = quantile(HOfx,[.2 .8]);
HOlo = HOfx<QHO(1);
HOhi = HOfx>QHO(2);
HOhilo = HOlo | HOhi;
HOloN  = HOy(HOlo) == TL(HOlo);
HOhiN  = HOy(HOhi) == TL(HOhi);
HOloMu  = nanmean(HOloN);
HOhiMu  = nanmean(HOhiN);
HOloPop = nanmean(HOlo);
HOhiPop = nanmean(HOhi);




% CTRLS
%----------------------------------------
CATRMEAN(vi) = case_tr_mu;
CATRLOMU(vi) = case_tr_lo_mu;
CATRHIMU(vi) = case_tr_hi_mu;
CATRLOPO(vi) = case_tr_lo_pop;
CATRHIPO(vi) = case_tr_hi_pop;

CAHOMEAN(vi) = case_ho_mu;
CAHOLOMU(vi) = case_ho_lo_mu;
CAHOHIMU(vi) = case_ho_hi_mu;
CAHOLOPO(vi) = case_ho_lo_pop;
CAHOHIPO(vi) = case_ho_hi_pop;


% CTRLS
%----------------------------------------
COTRMEAN(vi) = ctrl_tr_mu;
COTRLOMU(vi) = ctrl_tr_lo_mu;
COTRHIMU(vi) = ctrl_tr_lo_mu;
COTRLOPO(vi) = ctrl_tr_lo_mu;
COTRHIPO(vi) = ctrl_tr_lo_mu;

COHOMEAN(vi) = ctrl_ho_mu;
COHOLOMU(vi) = ctrl_ho_lo_mu;
COHOHIMU(vi) = ctrl_ho_lo_mu;
COHOLOPO(vi) = ctrl_ho_lo_mu;
COHOHIPO(vi) = ctrl_ho_lo_mu;



clc;
disp('%=================================================================');
disp('IJ | vi | min(SNPi) | max(SNPi) | numel(BETA):'); 
disp([IJ vi min(SNPi) max(SNPi) numel(BETA)]);
disp('----------  TRAINING SET  ----------')
fprintf('TRAIN correct %0.0f%% overall(pop:100%%)\n' ,(TRmu .* 100))
fprintf('TRAIN correct %0.0f%% locon  (pop: %0.0f%%)\n\n',TRloMu.*100,TRloPop.*100)
fprintf('TRAIN correct %0.0f%% hicon  (pop: %0.0f%%)\n\n',TRhiMu.*100,TRhiPop.*100)

disp('----------  TESTING SET  ----------')
fprintf('TEST correct %0.0f%% overall(pop:100%%)\n' ,(HOmu .* 100))
fprintf('TEST correct %0.0f%% locon  (pop: %0.0f%%)\n\n',HOloMu.*100,HOloPop.*100)
fprintf('TEST correct %0.0f%% hicon  (pop: %0.0f%%)\n\n',HOhiMu.*100,HOhiPop.*100)
disp('%=================================================================');


end
%===========================================

% FILL IN NAN VALUES WIHT NANMEAN
% CASES
%----------------------------------------
CATRMEAN(isnan(CATRMEAN)) = nanmean(CATRMEAN);
CATRLOMU(isnan(CATRLOMU)) = nanmean(CATRLOMU);
CATRHIMU(isnan(CATRHIMU)) = nanmean(CATRHIMU);
CATRLOPO(isnan(CATRLOPO)) = nanmean(CATRLOPO);
CATRHIPO(isnan(CATRHIPO)) = nanmean(CATRHIPO);

CAHOMEAN(isnan(CAHOMEAN)) = nanmean(CAHOMEAN);
CAHOLOMU(isnan(CAHOLOMU)) = nanmean(CAHOLOMU);
CAHOHIMU(isnan(CAHOHIMU)) = nanmean(CAHOHIMU);
CAHOLOPO(isnan(CAHOLOPO)) = nanmean(CAHOLOPO);
CAHOHIPO(isnan(CAHOHIPO)) = nanmean(CAHOHIPO);


% CTRLS
%----------------------------------------
COTRMEAN(isnan(COTRMEAN)) = nanmean(COTRMEAN);
COTRLOMU(isnan(COTRLOMU)) = nanmean(COTRLOMU);
COTRHIMU(isnan(COTRHIMU)) = nanmean(COTRHIMU);
COTRLOPO(isnan(COTRLOPO)) = nanmean(COTRLOPO);
COTRHIPO(isnan(COTRHIPO)) = nanmean(COTRHIPO);

COHOMEAN(isnan(COHOMEAN)) = nanmean(COHOMEAN);
COHOLOMU(isnan(COHOLOMU)) = nanmean(COHOLOMU);
COHOHIMU(isnan(COHOHIMU)) = nanmean(COHOHIMU);
COHOLOPO(isnan(COHOLOPO)) = nanmean(COHOLOPO);
COHOHIPO(isnan(COHOHIPO)) = nanmean(COHOHIPO);


%==========================================================================
%%                               HILO-LOHI SAVER
%==========================================================================


% CASES
%----------------------------------------
LOOPDATA.CATRMEAN(1:P.Ndots,IJ) = CATRMEAN;
LOOPDATA.CATRLOMU(1:P.Ndots,IJ) = CATRLOMU;
LOOPDATA.CATRHIMU(1:P.Ndots,IJ) = CATRHIMU;
LOOPDATA.CATRLOPO(1:P.Ndots,IJ) = CATRLOPO;
LOOPDATA.CATRHIPO(1:P.Ndots,IJ) = CATRHIPO;

LOOPDATA.CAHOMEAN(1:P.Ndots,IJ) = CAHOMEAN;
LOOPDATA.CAHOLOMU(1:P.Ndots,IJ) = CAHOLOMU;
LOOPDATA.CAHOHIMU(1:P.Ndots,IJ) = CAHOHIMU;
LOOPDATA.CAHOLOPO(1:P.Ndots,IJ) = CAHOLOPO;
LOOPDATA.CAHOHIPO(1:P.Ndots,IJ) = CAHOHIPO;


% CTRLS
%----------------------------------------
LOOPDATA.COTRMEAN(1:P.Ndots,IJ) = COTRMEAN;
LOOPDATA.COTRLOMU(1:P.Ndots,IJ) = COTRLOMU;
LOOPDATA.COTRHIMU(1:P.Ndots,IJ) = COTRHIMU;
LOOPDATA.COTRLOPO(1:P.Ndots,IJ) = COTRLOPO;
LOOPDATA.COTRHIPO(1:P.Ndots,IJ) = COTRHIPO;

LOOPDATA.COHOMEAN(1:P.Ndots,IJ) = COHOMEAN;
LOOPDATA.COHOLOMU(1:P.Ndots,IJ) = COHOLOMU;
LOOPDATA.COHOHIMU(1:P.Ndots,IJ) = COHOHIMU;
LOOPDATA.COHOLOPO(1:P.Ndots,IJ) = COHOLOPO;
LOOPDATA.COHOHIPO(1:P.Ndots,IJ) = COHOHIPO;



end
%% SAVE LOOP DATA


%------------------------------------------%
pause(1)
% dt=char(datetime(datetime,'Format','yyyy-MM-dd-HH-mm-ss'));
save(P.savepathfile,'LOOPDATA','P');
disp('File Saved.')
pause(1)
%------------------------------------------%
