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


P.mainmatfile = which('GENOSDATA.mat');
disp('dataset loaded')
clearvars -except P ADSP INFO


%% SET RUN OPTIONS & PATHS FOR DATA IMPORT/EXPORT
clc; clearvars -except P ADSP INFO


P.Nloops = 50;
P.Nvars = 200;
P.windowSize = 0;
P.Ndots = P.Nvars;
P.Lo2Hi = 1>0; %YES



P.basedir = 'F:\GENOSDATA\APOE_SUBGROUPS';
P.importdir = [P.basedir '\APOE_22_23_24_33_34_44\APOE_22_23_24_33_34_44_FISHP_FINAL'];
P.APOES = '22_23_24_33_34_44';
INFO.APOE = [22 23 24 33 34 44];

% P.basedir = 'F:\GENOSDATA\APOE_SUBGROUPS';
% P.importdir = [P.basedir '\APOE_22_23_24_34_44\APOE_22_23_24_34_44_FISHP_FINAL'];
% P.APOES = '22_23_24_34_44';
% INFO.APOE = [22 23 24 34 44];

% P.basedir = 'F:\GENOSDATA\APOE_SUBGROUPS';
% P.importdir = [P.basedir '\APOE_33\APOE_33_FISHP_FINAL'];
% P.APOES = '33';
% INFO.APOE = [33];

clearvars -except P ADSP INFO



%% VALIDATE IMPORT OPTIONS & GET PATHS TO EACH FISHP.MAT FILE

P.FILES.w = what(P.importdir);
P.Nmatfiles = numel(P.FILES.w.mat);
disp(P.FILES.w.mat); disp(P.Nmatfiles);



if P.Nloops > P.Nmatfiles
disp('ABORTING: NOT ENOUGH MAT FILES TO RUN THAT MANY LOOPS');
return; 
end

clearvars -except P ADSP INFO


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


clc; clearvars -except P ADSP INFO PHEN LOCI CASE CTRL USNP
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
clearvars -except P ADSP INFO PHEN LOCI CASE CTRL USNP


% NEURAL NETWORKS CONFUSION STATS
%----------------------------------------
LOOPDATA.TR_ALL_STATS = zeros(P.Ndots,9,P.Nloops);
LOOPDATA.HO_ALL_STATS = zeros(P.Ndots,9,P.Nloops);
LOOPDATA.TR_TOP_STATS = zeros(P.Ndots,9,P.Nloops);
LOOPDATA.HO_TOP_STATS = zeros(P.Ndots,9,P.Nloops);



% CASES NEURAL NETWORKS
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


% CTRLS NEURAL NETWORKS
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


% CASES LOGISTIC REGRESSION
%----------------------------------------
LOOPDATA.CATRMEAN_LR = zeros(P.Ndots,P.Nloops);
LOOPDATA.CATRLOMU_LR = zeros(P.Ndots,P.Nloops);
LOOPDATA.CATRHIMU_LR = zeros(P.Ndots,P.Nloops);
LOOPDATA.CATRLOPO_LR = zeros(P.Ndots,P.Nloops);
LOOPDATA.CATRHIPO_LR = zeros(P.Ndots,P.Nloops);

LOOPDATA.CAHOMEAN_LR = zeros(P.Ndots,P.Nloops);
LOOPDATA.CAHOLOMU_LR = zeros(P.Ndots,P.Nloops);
LOOPDATA.CAHOHIMU_LR = zeros(P.Ndots,P.Nloops);
LOOPDATA.CAHOLOPO_LR = zeros(P.Ndots,P.Nloops);
LOOPDATA.CAHOHIPO_LR = zeros(P.Ndots,P.Nloops);


% CTRLS LOGISTIC REGRESSION
%----------------------------------------
LOOPDATA.COTRMEAN_LR = zeros(P.Ndots,P.Nloops);
LOOPDATA.COTRLOMU_LR = zeros(P.Ndots,P.Nloops);
LOOPDATA.COTRHIMU_LR = zeros(P.Ndots,P.Nloops);
LOOPDATA.COTRLOPO_LR = zeros(P.Ndots,P.Nloops);
LOOPDATA.COTRHIPO_LR = zeros(P.Ndots,P.Nloops);

LOOPDATA.COHOMEAN_LR = zeros(P.Ndots,P.Nloops);
LOOPDATA.COHOLOMU_LR = zeros(P.Ndots,P.Nloops);
LOOPDATA.COHOHIMU_LR = zeros(P.Ndots,P.Nloops);
LOOPDATA.COHOLOPO_LR = zeros(P.Ndots,P.Nloops);
LOOPDATA.COHOHIPO_LR = zeros(P.Ndots,P.Nloops);


%==========================================================================
%% RUN MAIN LOOP
%==========================================================================
for IJ = 1:P.Nloops
%%
clearvars -except P ADSP INFO PHEN LOCI CASE CTRL USNP LOOPDATA IJ



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



% Save some data for the record
INFO.TOPLOCI{IJ}   = VLOCI(1:500,:);
INFO.PHETRCASE{IJ} = VTRCASE;
INFO.PHETRCTRL{IJ} = VTRCTRL;
INFO.PHETECASE{IJ} = VTECASE;
INFO.PHETECTRL{IJ} = VTECTRL;


%==========================================================================
%                       NEURAL NET
%==========================================================================


STAT.TR_ALL_STATS = zeros(P.Ndots,9);
STAT.HO_ALL_STATS = zeros(P.Ndots,9);
STAT.TR_TOP_STATS = zeros(P.Ndots,9);
STAT.HO_TOP_STATS = zeros(P.Ndots,9);


% CASES NEURAL NETS
%----------------------------------------
STAT.CATRMEAN = zeros(P.Ndots,1);
STAT.CATRLOMU = zeros(P.Ndots,1);
STAT.CATRHIMU = zeros(P.Ndots,1);
STAT.CATRLOPO = zeros(P.Ndots,1);
STAT.CATRHIPO = zeros(P.Ndots,1);

STAT.CAHOMEAN = zeros(P.Ndots,1);
STAT.CAHOLOMU = zeros(P.Ndots,1);
STAT.CAHOHIMU = zeros(P.Ndots,1);
STAT.CAHOLOPO = zeros(P.Ndots,1);
STAT.CAHOHIPO = zeros(P.Ndots,1);


% CTRLS NEURAL NETS
%----------------------------------------
STAT.COTRMEAN = zeros(P.Ndots,1);
STAT.COTRLOMU = zeros(P.Ndots,1);
STAT.COTRHIMU = zeros(P.Ndots,1);
STAT.COTRLOPO = zeros(P.Ndots,1);
STAT.COTRHIPO = zeros(P.Ndots,1);

STAT.COHOMEAN = zeros(P.Ndots,1);
STAT.COHOLOMU = zeros(P.Ndots,1);
STAT.COHOHIMU = zeros(P.Ndots,1);
STAT.COHOLOPO = zeros(P.Ndots,1);
STAT.COHOHIPO = zeros(P.Ndots,1);




% CASES LOGISTIC REGRESSION
%----------------------------------------
STAT.CATRMEAN_LR = zeros(P.Ndots,1);
STAT.CATRLOMU_LR = zeros(P.Ndots,1);
STAT.CATRHIMU_LR = zeros(P.Ndots,1);
STAT.CATRLOPO_LR = zeros(P.Ndots,1);
STAT.CATRHIPO_LR = zeros(P.Ndots,1);

STAT.CAHOMEAN_LR = zeros(P.Ndots,1);
STAT.CAHOLOMU_LR = zeros(P.Ndots,1);
STAT.CAHOHIMU_LR = zeros(P.Ndots,1);
STAT.CAHOLOPO_LR = zeros(P.Ndots,1);
STAT.CAHOHIPO_LR = zeros(P.Ndots,1);


% CTRLS LOGISTIC REGRESSION
%----------------------------------------
STAT.COTRMEAN_LR = zeros(P.Ndots,1);
STAT.COTRLOMU_LR = zeros(P.Ndots,1);
STAT.COTRHIMU_LR = zeros(P.Ndots,1);
STAT.COTRLOPO_LR = zeros(P.Ndots,1);
STAT.COTRHIPO_LR = zeros(P.Ndots,1);

STAT.COHOMEAN_LR = zeros(P.Ndots,1);
STAT.COHOLOMU_LR = zeros(P.Ndots,1);
STAT.COHOHIMU_LR = zeros(P.Ndots,1);
STAT.COHOLOPO_LR = zeros(P.Ndots,1);
STAT.COHOHIPO_LR = zeros(P.Ndots,1);


%% =====================  NEURAL NET LOOP =================================
for vi = 1:P.Ndots
%%
% clc; close all; 
% clearvars -except P ADSP INFO PHEN LOCI CASE CTRL USNP LOOPDATA IJ...
% vi VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL INFO STAT

    SNPi = 1:vi;



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



    [PVTR, ~, ~] = mkmxwin(XCASE,XCTRL,XUSNP,TRPHE,[-1 -0 2 3]);
    [PVHO, ~, ~] = mkmxwin(XCASE,XCTRL,XUSNP,TEPHE,[-1 -0 2 3]);
    %==========================================================================


    % PVMX(: , 1)  =  PHE.SRR;        % COL1: ID
    % PVMX(: , 2)  =  PHE.AD;         % COL2: AD
    % PVMX(: , 3)  =  PHE.COHORTNUM;  % COL3: COHORT
    % PVMX(: , 4)  =  PHE.AGE;        % COL4: AGE
    % PVMX(: , 5)  =  PHE.APOE;       % COL5: APOE
    % PVMX(: , 6)  =  PHE.SEX;        % COL7: SEX
    % PVMX(: , 7)  =  PHE.BRAAK;      % COL6: BRAAK
    % PVMX(: , 8)  =  PHE.BRAAK;      % COL6: BRAAK
    % PVMX(: , 9)  =  PHE.BRAAK;      % COL6: BRAAK
    % PVTR(:,end-5:end) = PVTR(:,[1 3 4 5 6 7]);
    % PVHO(:,end-5:end) = PVHO(:,[1 3 4 5 6 7]);


    %==========================================================================
    %       TRAIN NEURAL NETWORK CLASSIFIER
    %==========================================================================

    TR_MX = PVTR(:,10:end);
    HO_MX = PVHO(:,10:end);

    TR_DL = dummyvar(categorical((PVTR(:,2)==1)==0));
    HO_DL = dummyvar(categorical((PVHO(:,2)==1)==0));
    TR_VL = TR_DL(:,1);
    HO_VL = HO_DL(:,1);


    NN = patternnet([50 20],'trainscg','crossentropy');

    NN.trainParam.max_fail = 50;
    NN.trainParam.showWindow = 0;
    NN.performParam.regularization = 0.1;
    NN.performParam.normalization = 'none';




    % TRAIN NEURAL NETWORK
    net = train(NN,TR_MX',TR_DL');
    %Mdl = fitrsvm(TX',TLL');
    
%%


    TR_ACT = net(TR_MX');
    HO_ACT = net(HO_MX');

    TR_TOPq = quantile(abs(TR_ACT(:)-.5),.75);
    HO_TOPq = quantile(abs(HO_ACT(:)-.5),.75);
    TR_TOPi = abs(TR_ACT(1,:)-.5)>TR_TOPq;
    HO_TOPi = abs(HO_ACT(1,:)-.5)>HO_TOPq;

    [TR_ALL_ERR, TR_ALL_CONXN,~, TR_ALL_CONXP] = confusion(TR_DL',TR_ACT);
    TR_ALL_COR = 1-TR_ALL_ERR;

    [TR_TOP_ERR, TR_TOP_CONXN,~, TR_TOP_CONXP] = confusion(...
    TR_DL(TR_TOPi,:)',TR_ACT(:,TR_TOPi)); TR_TOP_COR = 1-TR_TOP_ERR;

    [HO_ALL_ERR, HO_ALL_CONXN,~, HO_ALL_CONXP] = confusion(HO_DL',HO_ACT);
    HO_ALL_COR = 1-HO_ALL_ERR;

    [HO_TOP_ERR, HO_TOP_CONXN,~, HO_TOP_CONXP] = confusion(...
    HO_DL(HO_TOPi,:)',HO_ACT(:,HO_TOPi)); HO_TOP_COR = 1-HO_TOP_ERR;


    TR_ALL_STATS = [TR_ALL_COR TR_ALL_CONXN(:)' TR_ALL_CONXP(1,:)];
    HO_ALL_STATS = [HO_ALL_COR HO_ALL_CONXN(:)' HO_ALL_CONXP(1,:)];
    TR_TOP_STATS = [TR_TOP_COR TR_TOP_CONXN(:)' TR_TOP_CONXP(1,:)];
    HO_TOP_STATS = [HO_TOP_COR HO_TOP_CONXN(:)' HO_TOP_CONXP(1,:)];

    STAT.TR_ALL_STATS(vi,:) = TR_ALL_STATS;
    STAT.HO_ALL_STATS(vi,:) = HO_ALL_STATS;
    STAT.TR_TOP_STATS(vi,:) = TR_TOP_STATS;
    STAT.HO_TOP_STATS(vi,:) = HO_TOP_STATS;



    %=====================================================
    %      NEURAL NETWORK MODEL PERFORMANCE
    %=====================================================

    % GET TRAINED NEURAL NETWORK ACTIVATIONS
    TR_HACT = TR_ACT - .5;
    HO_HACT = HO_ACT - .5;


    % ALL: SPLIT PREDICTIONS FOR CASES AND CONTROLS
    ACT_TR_CASE = TR_HACT(1,TR_VL==1) .*  1;
    ACT_TR_CTRL = TR_HACT(1,TR_VL==0) .* -1;
    ACT_HO_CASE = HO_HACT(1,HO_VL==1) .*  1;
    ACT_HO_CTRL = HO_HACT(1,HO_VL==0) .* -1;


    % ALL: DESCRETIZE PREDICTIONS
    TF_TR_CASE = ACT_TR_CASE >= 0;
    TF_TR_CTRL = ACT_TR_CTRL >= 0;
    TF_HO_CASE = ACT_HO_CASE >= 0;
    TF_HO_CTRL = ACT_HO_CTRL >= 0;


    % ALL: PROPORTION CORRECT
    MU_TR_CASE = nanmean(TF_TR_CASE);
    MU_TR_CTRL = nanmean(TF_TR_CTRL);
    MU_HO_CASE = nanmean(TF_HO_CASE);
    MU_HO_CTRL = nanmean(TF_HO_CTRL);



    % ALL/HIGH: HAS HIGH ACTIVATION?
    HIPO_TR_CASE = quantile(abs(ACT_TR_CASE),.75);
    HIPO_TR_CTRL = quantile(abs(ACT_TR_CTRL),.75);
    HIPO_HO_CASE = quantile(abs(ACT_HO_CASE),.75);
    HIPO_HO_CTRL = quantile(abs(ACT_HO_CTRL),.75);
    ISHI_TR_CASE = abs(ACT_TR_CASE) > HIPO_TR_CASE;
    ISHI_TR_CTRL = abs(ACT_TR_CTRL) > HIPO_TR_CTRL;
    ISHI_HO_CASE = abs(ACT_HO_CASE) > HIPO_HO_CASE;
    ISHI_HO_CTRL = abs(ACT_HO_CTRL) > HIPO_HO_CTRL;


    % HIGH: PROPORTION CORRECT
    HIMU_TR_CASE = nanmean(ACT_TR_CASE(ISHI_TR_CASE) > 0);
    HIMU_TR_CTRL = nanmean(ACT_TR_CTRL(ISHI_TR_CTRL) > 0);
    HIMU_HO_CASE = nanmean(ACT_HO_CASE(ISHI_HO_CASE) > 0);
    HIMU_HO_CTRL = nanmean(ACT_HO_CTRL(ISHI_HO_CTRL) > 0);


    % CASE
    %----------------------------------------
    STAT.CATRMEAN(vi) = MU_TR_CASE;
    STAT.CATRHIMU(vi) = HIMU_TR_CASE;
    STAT.CATRHIPO(vi) = HIPO_TR_CASE;

    STAT.CAHOMEAN(vi) = MU_HO_CASE;
    STAT.CAHOHIMU(vi) = HIMU_HO_CASE;
    STAT.CAHOHIPO(vi) = HIPO_HO_CASE;


    % CTRL
    %----------------------------------------
    STAT.COTRMEAN(vi) = MU_TR_CTRL;
    STAT.COTRHIMU(vi) = HIMU_TR_CTRL;
    STAT.COTRHIPO(vi) = HIPO_TR_CTRL;

    STAT.COHOMEAN(vi) = MU_HO_CTRL;
    STAT.COHOHIMU(vi) = HIMU_HO_CTRL;
    STAT.COHOHIPO(vi) = HIPO_HO_CTRL;

    






    %==========================================================================
    %%       TRAIN LOGISTIC REGRESSION CLASSIFIER
    %==========================================================================

    TXX_LR = PVTR(:,10:end);
    HXX_LR = PVHO(:,10:end);
    TX_LR = [ones(size(TXX_LR,1),1) TXX_LR]; % ADD AN INTERCEPT COLUMN
    HX_LR = [ones(size(HXX_LR,1),1) HXX_LR]; % ADD AN INTERCEPT COLUMN
    TL_LR = PVTR(:,2);
    HL_LR = PVHO(:,2);

    % PERFORM THE SO-CALLED MACHINE LEARNING STEP
    BETA = pinv(TX_LR' * TX_LR) * (TX_LR' * TL_LR);


    %=====================================================
    %      LOGISTIC REGRESSION MODEL PERFORMANCE
    %=====================================================

    % GET TRAINED LINEAR MODEL PREDICTIONS (ACTIVATIONS)
    TRfx_LR = nansum( TX_LR .* BETA' ,2) - .5;
    HOfx_LR = nansum( HX_LR .* BETA' ,2) - .5;


    % ALL: SPLIT PREDICTIONS FOR CASES AND CONTROLS
    ACT_TR_CASE_LR = TRfx_LR(TL_LR==1) .*  1;
    ACT_TR_CTRL_LR = TRfx_LR(TL_LR==0) .* -1;
    ACT_HO_CASE_LR = HOfx_LR(HL_LR==1) .*  1;
    ACT_HO_CTRL_LR = HOfx_LR(HL_LR==0) .* -1;


    % ALL: DESCRETIZE PREDICTIONS
    TF_TR_CASE_LR = ACT_TR_CASE_LR >= 0;
    TF_TR_CTRL_LR = ACT_TR_CTRL_LR >  0;
    TF_HO_CASE_LR = ACT_HO_CASE_LR >= 0;
    TF_HO_CTRL_LR = ACT_HO_CTRL_LR >  0;


    % ALL: PROPORTION CORRECT
    MU_TR_CASE_LR = nanmean(TF_TR_CASE_LR);
    MU_TR_CTRL_LR = nanmean(TF_TR_CTRL_LR);
    MU_HO_CASE_LR = nanmean(TF_HO_CASE_LR);
    MU_HO_CTRL_LR = nanmean(TF_HO_CTRL_LR);



    % ALL/HIGH: HAS HIGH ACTIVATION?
    HIPO_TR_CASE_LR = quantile(abs(ACT_TR_CASE_LR),.75);
    HIPO_TR_CTRL_LR = quantile(abs(ACT_TR_CTRL_LR),.75);
    HIPO_HO_CASE_LR = quantile(abs(ACT_HO_CASE_LR),.75);
    HIPO_HO_CTRL_LR = quantile(abs(ACT_HO_CTRL_LR),.75);
    ISHI_TR_CASE_LR = abs(ACT_TR_CASE_LR) > HIPO_TR_CASE_LR;
    ISHI_TR_CTRL_LR = abs(ACT_TR_CTRL_LR) > HIPO_TR_CTRL_LR;
    ISHI_HO_CASE_LR = abs(ACT_HO_CASE_LR) > HIPO_HO_CASE_LR;
    ISHI_HO_CTRL_LR = abs(ACT_HO_CTRL_LR) > HIPO_HO_CTRL_LR;


    % HIGH: PROPORTION CORRECT
    HIMU_TR_CASE_LR = mean(ACT_TR_CASE_LR(ISHI_TR_CASE_LR) > 0);
    HIMU_TR_CTRL_LR = mean(ACT_TR_CTRL_LR(ISHI_TR_CTRL_LR) > 0);
    HIMU_HO_CASE_LR = mean(ACT_HO_CASE_LR(ISHI_HO_CASE_LR) > 0);
    HIMU_HO_CTRL_LR = mean(ACT_HO_CTRL_LR(ISHI_HO_CTRL_LR) > 0);



    % CASE
    %----------------------------------------
    STAT.CATRMEAN_LR(vi) = MU_TR_CASE_LR;
    STAT.CATRHIMU_LR(vi) = HIMU_TR_CASE_LR;
    STAT.CATRHIPO_LR(vi) = HIPO_TR_CASE_LR;

    STAT.CAHOMEAN_LR(vi) = MU_HO_CASE_LR;
    STAT.CAHOHIMU_LR(vi) = HIMU_HO_CASE_LR;
    STAT.CAHOHIPO_LR(vi) = HIPO_HO_CASE_LR;


    % CTRL
    %----------------------------------------
    STAT.COTRMEAN_LR(vi) = MU_TR_CTRL_LR;
    STAT.COTRHIMU_LR(vi) = HIMU_TR_CTRL_LR;
    STAT.COTRHIPO_LR(vi) = HIPO_TR_CTRL_LR;

    STAT.COHOMEAN_LR(vi) = MU_HO_CTRL_LR;
    STAT.COHOHIMU_LR(vi) = HIMU_HO_CTRL_LR;
    STAT.COHOHIPO_LR(vi) = HIPO_HO_CTRL_LR;




    clc;
    disp('%=================================================================');
    disp('IJ | vi | min(SNPi) | max(SNPi) | numel(BETA):'); 
    disp([IJ vi min(SNPi) max(SNPi) numel(SNPi)]);

    disp('REGRESSION TRAINING'); disp('---------------------------');
    fprintf(['  CASE ALL: %0.0f%% \n  CTRL ALL: %0.0f%% \n'...
             '  CASE TOP: %0.0f%% (pop: %0.0f%%) \n  CTRL TOP: %0.0f%% (pop: %0.0f%%)\n\n'] ,...
        MU_TR_CASE_LR.*100, MU_TR_CTRL_LR.*100,...
        HIMU_TR_CASE_LR.*100,HIPO_TR_CASE_LR.*100,HIMU_TR_CTRL_LR.*100,HIPO_TR_CTRL_LR.*100)

    disp('REGRESSION HOLDOUT'); disp('---------------------------');
    fprintf(['  CASE ALL: %0.0f%% \n  CTRL ALL: %0.0f%% \n'...
             '  CASE TOP: %0.0f%% (pop: %0.0f%%) \n  CTRL TOP: %0.0f%% (pop: %0.0f%%)\n\n'] ,...
        MU_HO_CASE_LR.*100, MU_HO_CTRL_LR.*100,...
        HIMU_HO_CASE_LR.*100,HIPO_HO_CASE_LR.*100,HIMU_HO_CTRL_LR.*100,HIPO_HO_CTRL_LR.*100)

    disp('NEURAL NET TRAINING'); disp('---------------------------');
    fprintf(['  CASE ALL: %0.0f%% \n  CTRL ALL: %0.0f%% \n'...
             '  CASE TOP: %0.0f%% (pop: %0.0f%%) \n  CTRL TOP: %0.0f%% (pop: %0.0f%%)\n\n'] ,...
        MU_TR_CASE.*100, MU_TR_CTRL.*100,...
        HIMU_TR_CASE.*100,HIPO_TR_CASE.*100,HIMU_TR_CTRL.*100,HIPO_TR_CTRL.*100)

    disp('NEURAL NET HOLDOUT'); disp('---------------------------');
    fprintf(['  CASE ALL: %0.0f%% \n  CTRL ALL: %0.0f%% \n'...
             '  CASE TOP: %0.0f%% (pop: %0.0f%%) \n  CTRL TOP: %0.0f%% (pop: %0.0f%%)\n'] ,...
        MU_HO_CASE.*100, MU_HO_CTRL.*100,...
        HIMU_HO_CASE.*100,HIPO_HO_CASE.*100,HIMU_HO_CTRL.*100,HIPO_HO_CTRL.*100)
disp('%=================================================================');


end
%===========================================

% FILL IN NAN VALUES WITH NANMEAN

% CASES NEURAL NETS
%----------------------------------------
STAT.CATRMEAN(isnan(STAT.CATRMEAN)) = nanmean(STAT.CATRMEAN);
STAT.CATRHIMU(isnan(STAT.CATRHIMU)) = nanmean(STAT.CATRHIMU);
STAT.CATRHIPO(isnan(STAT.CATRHIPO)) = nanmean(STAT.CATRHIPO);

STAT.CAHOMEAN(isnan(STAT.CAHOMEAN)) = nanmean(STAT.CAHOMEAN);
STAT.CAHOHIMU(isnan(STAT.CAHOHIMU)) = nanmean(STAT.CAHOHIMU);
STAT.CAHOHIPO(isnan(STAT.CAHOHIPO)) = nanmean(STAT.CAHOHIPO);


% CTRLS NEURAL NETS
%----------------------------------------
STAT.COTRMEAN(isnan(STAT.COTRMEAN)) = nanmean(STAT.COTRMEAN);
STAT.COTRHIMU(isnan(STAT.COTRHIMU)) = nanmean(STAT.COTRHIMU);
STAT.COTRHIPO(isnan(STAT.COTRHIPO)) = nanmean(STAT.COTRHIPO);

STAT.COHOMEAN(isnan(STAT.COHOMEAN)) = nanmean(STAT.COHOMEAN);
STAT.COHOHIMU(isnan(STAT.COHOHIMU)) = nanmean(STAT.COHOHIMU);
STAT.COHOHIPO(isnan(STAT.COHOHIPO)) = nanmean(STAT.COHOHIPO);



% CASES LOGISTIC REGRESSION
%----------------------------------------
STAT.CATRMEAN_LR(isnan(STAT.CATRMEAN_LR)) = nanmean(STAT.CATRMEAN_LR);
STAT.CATRHIMU_LR(isnan(STAT.CATRHIMU_LR)) = nanmean(STAT.CATRHIMU_LR);
STAT.CATRHIPO_LR(isnan(STAT.CATRHIPO_LR)) = nanmean(STAT.CATRHIPO_LR);

STAT.CAHOMEAN_LR(isnan(STAT.CAHOMEAN_LR)) = nanmean(STAT.CAHOMEAN_LR);
STAT.CAHOHIMU_LR(isnan(STAT.CAHOHIMU_LR)) = nanmean(STAT.CAHOHIMU_LR);
STAT.CAHOHIPO_LR(isnan(STAT.CAHOHIPO_LR)) = nanmean(STAT.CAHOHIPO_LR);


% CTRLS LOGISTIC REGRESSION
%----------------------------------------
STAT.COTRMEAN_LR(isnan(STAT.COTRMEAN_LR)) = nanmean(STAT.COTRMEAN_LR);
STAT.COTRHIMU_LR(isnan(STAT.COTRHIMU_LR)) = nanmean(STAT.COTRHIMU_LR);
STAT.COTRHIPO_LR(isnan(STAT.COTRHIPO_LR)) = nanmean(STAT.COTRHIPO_LR);

STAT.COHOMEAN_LR(isnan(STAT.COHOMEAN_LR)) = nanmean(STAT.COHOMEAN_LR);
STAT.COHOHIMU_LR(isnan(STAT.COHOHIMU_LR)) = nanmean(STAT.COHOHIMU_LR);
STAT.COHOHIPO_LR(isnan(STAT.COHOHIPO_LR)) = nanmean(STAT.COHOHIPO_LR);


%==========================================================================
%%                               HILO-LOHI SAVER
%==========================================================================



% NEURAL NETWORKS CONFUSION STATS
%----------------------------------------
LOOPDATA.TR_ALL_STATS(1:P.Ndots,1:9,IJ) = STAT.TR_ALL_STATS;
LOOPDATA.HO_ALL_STATS(1:P.Ndots,1:9,IJ) = STAT.HO_ALL_STATS;
LOOPDATA.TR_TOP_STATS(1:P.Ndots,1:9,IJ) = STAT.TR_TOP_STATS;
LOOPDATA.HO_TOP_STATS(1:P.Ndots,1:9,IJ) = STAT.HO_TOP_STATS;



% CASES NEURAL NETS
%----------------------------------------
LOOPDATA.CATRMEAN(1:P.Ndots,IJ) = STAT.CATRMEAN;
LOOPDATA.CATRHIMU(1:P.Ndots,IJ) = STAT.CATRHIMU;
LOOPDATA.CATRHIPO(1:P.Ndots,IJ) = STAT.CATRHIPO;

LOOPDATA.CAHOMEAN(1:P.Ndots,IJ) = STAT.CAHOMEAN;
LOOPDATA.CAHOHIMU(1:P.Ndots,IJ) = STAT.CAHOHIMU;
LOOPDATA.CAHOHIPO(1:P.Ndots,IJ) = STAT.CAHOHIPO;


% CTRLS NEURAL NETS
%----------------------------------------
LOOPDATA.COTRMEAN(1:P.Ndots,IJ) = STAT.COTRMEAN;
LOOPDATA.COTRHIMU(1:P.Ndots,IJ) = STAT.COTRHIMU;
LOOPDATA.COTRHIPO(1:P.Ndots,IJ) = STAT.COTRHIPO;

LOOPDATA.COHOMEAN(1:P.Ndots,IJ) = STAT.COHOMEAN;
LOOPDATA.COHOHIMU(1:P.Ndots,IJ) = STAT.COHOHIMU;
LOOPDATA.COHOHIPO(1:P.Ndots,IJ) = STAT.COHOHIPO;


% CASES LOGISTIC REGRESSION
%----------------------------------------
LOOPDATA.CATRMEAN_LR(1:P.Ndots,IJ) = STAT.CATRMEAN_LR;
LOOPDATA.CATRHIMU_LR(1:P.Ndots,IJ) = STAT.CATRHIMU_LR;
LOOPDATA.CATRHIPO_LR(1:P.Ndots,IJ) = STAT.CATRHIPO_LR;

LOOPDATA.CAHOMEAN_LR(1:P.Ndots,IJ) = STAT.CAHOMEAN_LR;
LOOPDATA.CAHOHIMU_LR(1:P.Ndots,IJ) = STAT.CAHOHIMU_LR;
LOOPDATA.CAHOHIPO_LR(1:P.Ndots,IJ) = STAT.CAHOHIPO_LR;


% CTRLS LOGISTIC REGRESSION
%----------------------------------------
LOOPDATA.COTRMEAN_LR(1:P.Ndots,IJ) = STAT.COTRMEAN_LR;
LOOPDATA.COTRHIMU_LR(1:P.Ndots,IJ) = STAT.COTRHIMU_LR;
LOOPDATA.COTRHIPO_LR(1:P.Ndots,IJ) = STAT.COTRHIPO_LR;

LOOPDATA.COHOMEAN_LR(1:P.Ndots,IJ) = STAT.COHOMEAN_LR;
LOOPDATA.COHOHIMU_LR(1:P.Ndots,IJ) = STAT.COHOHIMU_LR;
LOOPDATA.COHOHIPO_LR(1:P.Ndots,IJ) = STAT.COHOHIPO_LR;



end
%% SAVE LOOP DATA


%------------------------------------------%
pause(1)
% dt=char(datetime(datetime,'Format','yyyy-MM-dd-HH-mm-ss'));
P.savepathfile = [P.basedir filesep 'APOE_' P.APOES '_NNHILO.mat'];
save(P.savepathfile,'LOOPDATA','P','INFO');
disp('File saved...'); disp(P.savepathfile)
pause(1)
%------------------------------------------%
