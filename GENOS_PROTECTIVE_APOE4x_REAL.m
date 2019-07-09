%==========================================================================
%% STEP-1: LOAD THE DATASET
%==========================================================================
close all; clear; clc; rng('shuffle');
P.home = fileparts(which('GENOS.m')); cd(P.home);
P.funs = [P.home filesep 'genos_functions'];
P.mfuns = [P.funs filesep 'genos_main_functions'];
P.other = [P.home filesep 'genos_other'];
%P.data = [P.home filesep 'genos_data'];
P.data = 'F:\GENOSDATA';
addpath(join(string(struct2cell(P)),pathsep,1))
cd(P.home); P.f = filesep;


ADSP = load('GENOSDATAFINAL.mat');


P.mainmatfile = which('GENOSDATAFINAL.mat');
disp('dataset loaded')
clearvars -except P ADSP INFO


%==========================================================================
%% SET RUN OPTIONS & PATHS FOR DATA IMPORT/EXPORT
%==========================================================================
clc; clearvars -except P ADSP INFO


P.doPRO = 1;
P.doPLO = 0;
P.doORLO = 0;
P.doORHI = 0;
P.doSYN = 0;
P.NGeneStart = 1;
P.NGeneEnd = 500;
P.NGenes = P.NGeneEnd - P.NGeneStart + 1;
P.Nloops = 10;
P.FileStart = 1;
P.Nvars = 200;
P.windowSize = 50;
P.Ndots = 5;
P.Lo2Hi = true;
P.RemoveGenesByName = false;
P.RemoveBadGenes = false;
P.TargetSNP_Path = [P.data P.f 'GENOX' P.f 'CSVSNP'];
P.Output_Path = [P.data P.f 'GENOX'];



clearvars -except P ADSP INFO
%==========================================================================
%% SET PATH TO (APOE DEPENDENT) SUBSAMPLE FOLDERS AND GET MAT PATHS
%==========================================================================
clc; clearvars -except P ADSP INFO



P.APOExx_path = [P.data P.f 'APOE' P.f 'APOE_22_23_24_33_34_44' P.f 'APOE_22_23_24_33_34_44_FISHP'];
P.APOExxFILES.sets = what(P.APOExx_path);
P.APOExxN = numel(P.APOExxFILES.sets.mat);


P.APOE4x_path = [P.data P.f 'APOE' P.f 'APOE_34_44' P.f 'APOE_34_44_FISHP'];
P.APOE4xFILES.sets = what(P.APOE4x_path);
P.APOE4xN = numel(P.APOE4xFILES.sets.mat);



fprintf('\n APOExx has N mat files: %0.f \n',P.APOExxN)
fprintf('\n APOE4x has N mat files: %0.f \n',P.APOE4xN)




clearvars -except P ADSP INFO
%==========================================================================
%% IMPORT TOP VARIANTS FROM EXCEL SHEET
%==========================================================================
clc; clearvars -except P ADSP INFO

PLOops  = detectImportOptions([P.TargetSNP_Path P.f 'PLO.csv']);
ORLOops = detectImportOptions([P.TargetSNP_Path P.f 'ORLO.csv']);
ORHIops = detectImportOptions([P.TargetSNP_Path P.f 'ORHI.csv']);
SYNops  = detectImportOptions([P.TargetSNP_Path P.f 'SYN.csv']);
PROops  = detectImportOptions([P.TargetSNP_Path P.f 'PRO.csv']);

SNPTAB.PLO  = readtable([P.TargetSNP_Path P.f 'PLO.csv'],PLOops);
SNPTAB.ORLO = readtable([P.TargetSNP_Path P.f 'ORLO.csv'],ORLOops);
SNPTAB.ORHI = readtable([P.TargetSNP_Path P.f 'ORHI.csv'],ORHIops);
SNPTAB.SYN  = readtable([P.TargetSNP_Path P.f 'SYN.csv'],SYNops);
SNPTAB.PRO  = readtable([P.TargetSNP_Path P.f 'PRO.csv'],PROops);

SNPTAB.PLO.CHRPOS  = uint64(SNPTAB.PLO.CHRPOS);
SNPTAB.ORLO.CHRPOS = uint64(SNPTAB.ORLO.CHRPOS);
SNPTAB.ORHI.CHRPOS = uint64(SNPTAB.ORHI.CHRPOS);
SNPTAB.SYN.CHRPOS  = uint64(SNPTAB.SYN.CHRPOS);
SNPTAB.PRO.CHRPOS  = uint64(SNPTAB.PRO.CHRPOS);


ADSP.SNPTAB = SNPTAB;

if P.doPLO == 1
    ADSP.SNP = ADSP.SNPTAB.PLO;
elseif P.doORLO == 1
    ADSP.SNP = ADSP.SNPTAB.ORLO;
elseif P.doORHI == 1
     ADSP.SNP = ADSP.SNPTAB.ORHI;
elseif P.doSYN == 1
     ADSP.SNP = ADSP.SNPTAB.SYN;
elseif P.doPRO == 1
     ADSP.SNP = ADSP.SNPTAB.PRO;
end



disp(head(ADSP.SNP));
clearvars -except P ADSP INFO
%==========================================================================
%%   CARBON COPY MAIN VARIABLES FROM ADSP.STRUCT
%==========================================================================
clc; clearvars -except P ADSP INFO SNPTAB



LOCI = ADSP.LOCI;
CASE = ADSP.CASE;
CTRL = ADSP.CTRL;
USNP = ADSP.USNP;
PHEN = ADSP.PHEN;



disp(head(PHEN)); disp(head(LOCI));
clearvars -except P ADSP INFO PHEN LOCI CASE CTRL USNP
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
clc; clearvars -except P ADSP INFO PHEN LOCI CASE CTRL USNP


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



% TRAINED MACHINE LEARNER
%----------------------------------------
LOOPDATA.NETS = cell(P.Nloops);
LOOPDATA.NETQ = cell(P.Nloops);
LOOPDATA.NETD = cell(P.Nloops);
LOOPDATA.SVMQ = nan(50,2,P.Nloops);
LOOPDATA.SVMD = nan(50,2,P.Nloops);


% NEURAL NETWORK & CONFUSION MATRIX
%----------------------------------------
LOOPDATA.CONFU_TR = nan(4,4,P.Nloops,6);
LOOPDATA.CONFU_HO = nan(4,4,P.Nloops,6);
LOOPDATA.PERF_TR = nan(P.Nloops,6,6);
LOOPDATA.PERF_HO = nan(P.Nloops,6,6);
LOOPDATA.AREA_TR = nan(50,5,P.Nloops,6);
LOOPDATA.AREA_HO = nan(50,5,P.Nloops,6);



% TARGET GENE & CHRPOS
%----------------------------------------
LOOPDATA.GENECHRPOS = cell(P.Nloops,2,P.NGenes);








clearvars -except P ADSP INFO PHEN LOCI CASE CTRL USNP LOOPDATA
%% LOOP OVER EACH GENE-CHRPOS TO PERTURB
%##########################################################################
%
% 
%
for kk = P.NGeneStart:P.NGeneEnd
%
%
% 
%##########################################################################
close all; clc;






%% LOOP OVER EACH OF 50 UNIQUE PARTICIPANT SUBSETS
%==========================================================================
% 
% 
% 
for ij = 1:P.Nloops
% 
% 
% 
%========================================================================== 
clearvars -except P ADSP INFO SNPTAB PHEN LOCI CASE CTRL USNP LOOPDATA kk ij
close all; clc; fprintf('\n\n | GENE LOOP: %.0f  \n | SUBSET LOOP: %.0f \n\n',kk,ij)


    % LOAD MAT DATA CONTAINING UNIQUE PARTICIPANT SUBSET

    MATDAT = load([P.APOExxFILES.sets.path  P.f  P.APOExxFILES.sets.mat{randi(50)} ]);


    



    % PULL OUT VARIABLES FROM MATDAT PARTICIPANT SUBSET
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




    %======================================================================
    % REMOVE SELECT GENES BY NAME
    %======================================================================
    %{
    if P.RemoveGenesByName
    ALZGENES = string(["APOE","TOMM40","CR1","BIN1","INPP5D","HLA-DRB1","TREM2",...
        "CD2AP","NYAP1g","EPHA1","PTK2B","CLU","SPI1h","MS4A2","PICALM","SORL1",...
        "FERMT2","SLC24A4","ABCA7","CASS4","ECHDC3","ACE","MEF2C","NME8","TYRO3"]);
    for nn = 1:numel(ALZGENES)
        x = strcmp(VLOCI.GENE,ALZGENES(nn));
        VLOCI(x,:) = [];
        VCASE(x) = [];
        VCTRL(x) = [];
        VUSNP(x) = [];
    end
    end
    if P.RemoveBadGenes
    BADGENES = string(["TYRO3","TOMM40"]);
    for nn = 1:numel(BADGENES)
        x = strcmp(VLOCI.GENE,BADGENES(nn));
        VLOCI(x,:) = [];
        VCASE(x) = [];
        VCTRL(x) = [];
        VUSNP(x) = [];
    end
    end
    VLOCI.VID  = (1:size(VLOCI,1))';


%     BADGENES = string(["TYRO3","TOMM40"]);
%     for nn = 1:numel(BADGENES)
%         x = strcmp(VLOCI.GENE,BADGENES(nn));
%         VLOCI(x,:) = [];
%         VCASE(x) = [];
%         VCTRL(x) = [];
%         VUSNP(x) = [];
%     end
    %}




    %==========================================================================
    % GET GENE-CHRLOC TO PERTURB AND ENSURE IT'S #1 AMONG SNPS
    %==========================================================================

    GENE   = ADSP.SNP.GENE{kk};
    CHRPOS = ADSP.SNP.CHRPOS(kk);

    VLOCI.TRFISHP(VLOCI.CHRPOS==CHRPOS) = 0;

    LOOPDATA.GENECHRPOS{ij,1} = GENE;
    LOOPDATA.GENECHRPOS{ij,2} = CHRPOS;




    %======================================================================
    % SORT VARIANTS BY EITHER TRFISHP|CHRPOS
    %======================================================================
    [~,j]  = sort(VLOCI.TRFISHP);
    VLOCI  = VLOCI(j,:);
    VCASE  = VCASE(j);
    VCTRL  = VCTRL(j);
    VUSNP  = VUSNP(j);


    disp(GENE);
    disp(VLOCI(1:5,1:9)); 
    pause(2);


    % Save some data for the record
    INFO.TOPLOCI{ij}   = VLOCI(1:500,:);
    INFO.PHETRCASE{ij} = VTRCASE;
    INFO.PHETRCTRL{ij} = VTRCTRL;
    INFO.PHETECASE{ij} = VTECASE;
    INFO.PHETECTRL{ij} = VTECTRL;





    %======================================================================
    %      MAKE  RECTANGLE  NEURAL NET  VARIANT MATRIX
    %======================================================================
    SNPi = 1:P.windowSize;


    % EXTRACT TOP-N NUMBER OF VARIANTS
    XLOCI  = VLOCI(SNPi,:);
    XCASE  = VCASE(SNPi);
    XCTRL  = VCTRL(SNPi);
    XUSNP  = VUSNP(SNPi);


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




    INFO.NNwts = [-1 -0 2 3];
    % [PVXt, VXt, DVXt, DXt, Yt, YDt] = mkdx(XLOCI,XCASE,XCTRL,XUSNP,TRPHE,INFO.NNwts);
    % [PVXh, VXh, DVXh, DXh, Yh, YDh] = mkdx(XLOCI,XCASE,XCTRL,XUSNP,TEPHE,INFO.NNwts);
    [PVXt, VXt, DVXt, DXt, Yt, YDt] = vxdx(XLOCI,XCASE,XCTRL,XUSNP,TRPHE,INFO.NNwts);
    [PVXh, VXh, DVXh, DXh, Yh, YDh] = vxdx(XLOCI,XCASE,XCTRL,XUSNP,TEPHE,INFO.NNwts);



%{
close all;
fh1 = figure('Units','normalized','Position',[.02 .04 .95 .85],'Color','w');
ax1 = axes('Position',[.05 .53 .93 .45],'Color','none');
ax2 = axes('Position',[.05 .04 .93 .45],'Color','none');
axes(ax1); bar( ([sum(VXt==-1);sum(VXt==0);sum(VXt==2);sum(VXt==3)]') ,'stacked');axis tight
axes(ax2); bar( ([sum(DXt==-1);sum(DXt==0);sum(DXt==2);sum(DXt==3)]') ,'stacked');axis tight
%}

%======================================================================
%%       TRAIN NEURAL NETWORK CLASSIFIER
%======================================================================
clc; clearvars -except P ADSP INFO PHEN LOCI CASE CTRL USNP LOOPDATA kk ij...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL SNPi XLOCI...
XCASE XCTRL XUSNP TRPHE TEPHE PVXt DVXt PVXh DVXh VXt VXh Yt Yh...
YDt YDh DXt DXh
%----------------------------------------------------------------------



TRAINVARS = {PVXt, VXt, DVXt, DXt, Yt, YDt};
TESTVARS  = {PVXh, VXh, DVXh, DXh, Yh, YDh};

PARAMS.LOOPS = 3;
PARAMS.doDX  = 0;
PARAMS.NNEURONS = [50,20];

net = nettrain(TRAINVARS,TESTVARS,PARAMS);

LOOPDATA.NETS{ij,1} = net;



[CONFU_TR, PERF_TR, AREA_TR] = confusionmx(Yt, VXt, net );
[CONFU_HO, PERF_HO, AREA_HO] = confusionmx(Yh, VXh, net ,1);
% [CONFU_TR, PERF_TR, PCUT_TR, AREA_TR] = nnperform(Yt', DXt', PVXt, net, 'cut',.5,'plt',0);
% [CONFU_HO, PERF_HO, PCUT_HO, AREA_HO] = nnperform(Yh', DXh', PVXh, net, 'cut',.5,'plt',0);





%% SYSTEMATICALLY PERTURB NEURAL NET BY ALTERING (1) SNP IN ALL PARTICIPANTS
% clc; clearvars -except P ADSP INFO PHEN LOCI CASE CTRL USNP LOOPDATA kk ij...
% VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL SNPi XLOCI...
% XCASE XCTRL XUSNP TRPHE TEPHE PVXt DVXt PVXh DVXh VXt VXh Yt Yh...
% YDt YDh DXt DXh netq netd net 
%----------------------------------------------------------------------
clearvars -except P ADSP INFO SNPTAB PHEN LOCI CASE CTRL USNP LOOPDATA kk ij...
VTRCASE VTRCTRL VTECASE VTECTRL VLOCI VCASE VCTRL VUSNP...
XLOCI XCASE XCTRL XUSNP TRPHE TEPHE net



GENE   = ADSP.SNP.GENE{kk}
CHRPOS = ADSP.SNP.CHRPOS(kk);




%------------------------------------------------------------------
% USING ONE-HOT-TRAINING-MATRIX (DVXt, DXt)
%-------------------------------------------
% TRAINING APOES:  22 23 24 33 34 44
% HOLDOUT  APOE:   {22|23|24|33|34|44}
% HOLDOUT TARSNP:  ALT|ALT
%------------------------------------------------------------------
PHEt = TRPHE(TRPHE.APOE ~= 11 ,:);
PHEh = TEPHE(TEPHE.APOE ~= 11 ,:);

[PVXt, VXt, DVXt, DXt, Yt, YDt] = vxdx(XLOCI,XCASE,XCTRL,XUSNP,PHEt,INFO.NNwts);
[PVXh, VXh, DVXh, DXh, Yh, YDh] = vxdx(XLOCI,XCASE,XCTRL,XUSNP,PHEh,INFO.NNwts);


% [[ APOE ]] UPDATE APOE STATUS IN HOLDOUT GROUP FOR ONE-HOT DUAL-VAR-MX (DVXh)
%-----------------------------------------------------
%{
% %--- SET E2{REF|REF}
% e2i = find(XLOCI.CHRPOS==190045412079);
% DVXh(:,(e2i*2+8))  = INFO.NNwts(1);
% DVXh(:,(e2i*2+9))  = INFO.NNwts(1);
% DXh(:,e2i*2-1)     = INFO.NNwts(1);
% DXh(:,e2i*2-0)     = INFO.NNwts(1);
% PVXh(:,(9+e2i))    = INFO.NNwts(1);
% VXh(:,e2i)         = INFO.NNwts(1);
% 
% %--- SET E4{ALT|REF}
% e4i = find(XLOCI.CHRPOS==190045411941);
% DVXh(:,(e4i*2+8))  = INFO.NNwts(3);
% DVXh(:,(e4i*2+9))  = INFO.NNwts(1);
% DXh(:,e4i*2-1)     = INFO.NNwts(3);
% DXh(:,e4i*2-0)     = INFO.NNwts(1);
% PVXh(:,(9+e4i))    = INFO.NNwts(3);
% VXh(:,e4i)         = INFO.NNwts(3);
%}


% [[ SNP ]] UPDATE TARGET SNP OF HOLDOUT GROUP IN ONE-HOT DUAL-VAR-MX (DVXh)
%-----------------------------------------------------
vi = find(XLOCI.CHRPOS==CHRPOS);
DVXh(:,(vi*2+8))  = INFO.NNwts(3);
DVXh(:,(vi*2+9))  = INFO.NNwts(4);
DXh(:,vi*2-1)     = INFO.NNwts(3);
DXh(:,vi*2-0)     = INFO.NNwts(4);
PVXh(:,(9+vi))    = INFO.NNwts(4);
VXh(:,vi)         = INFO.NNwts(4);



% NET: GET ACTIVATIONS FOR TRAINING AND HOLDOUT
[CONFU_TR, PERF_TR, AREA_TR] = confusionmx(Yt, VXt, net );
[CONFU_HO, PERF_HO, AREA_HO] = confusionmx(Yh, VXh, net );

LOOPDATA.CONFU_TR(:,:,ij, 1) = CONFU_TR;
LOOPDATA.CONFU_HO(:,:,ij, 1) = CONFU_HO;
LOOPDATA.PERF_TR(ij,:,    1) = PERF_TR;
LOOPDATA.PERF_HO(ij,:,    1) = PERF_HO;
LOOPDATA.AREA_TR(:,:,ij,  1) = AREA_TR;
LOOPDATA.AREA_HO(:,:,ij,  1) = AREA_HO;
%------------------------------------------------------------------













%------------------------------------------------------------------
% USING ONE-HOT-TRAINING-MATRIX (DVXt, DXt)
%-------------------------------------------
% TRAINING APOES:  22 23 24 33 34 44
% HOLDOUT  APOE:   {33}
% HOLDOUT TARSNP:  ALT|ALT
%------------------------------------------------------------------
PHEt = TRPHE(TRPHE.APOE == 33 ,:);
PHEh = TEPHE(TEPHE.APOE == 33 ,:);

[PVXt, VXt, DVXt, DXt, Yt, YDt] = vxdx(XLOCI,XCASE,XCTRL,XUSNP,PHEt,INFO.NNwts);
[PVXh, VXh, DVXh, DXh, Yh, YDh] = vxdx(XLOCI,XCASE,XCTRL,XUSNP,PHEh,INFO.NNwts);


% [[ APOE ]] UPDATE APOE STATUS IN HOLDOUT GROUP FOR ONE-HOT DUAL-VAR-MX (DVXh)
%-----------------------------------------------------
%{
%--- SET E2{REF|REF}
e2i = find(XLOCI.CHRPOS==190045412079);
DVXh(:,(e2i*2+8))  = INFO.NNwts(1);
DVXh(:,(e2i*2+9))  = INFO.NNwts(1);
DXh(:,e2i*2-1)     = INFO.NNwts(1);
DXh(:,e2i*2-0)     = INFO.NNwts(1);
PVXh(:,(9+e2i))    = INFO.NNwts(1);
VXh(:,e2i)         = INFO.NNwts(1);

%--- SET E4{ALT|REF}
e4i = find(XLOCI.CHRPOS==190045411941);
DVXh(:,(e4i*2+8))  = INFO.NNwts(3);
DVXh(:,(e4i*2+9))  = INFO.NNwts(1);
DXh(:,e4i*2-1)     = INFO.NNwts(3);
DXh(:,e4i*2-0)     = INFO.NNwts(1);
PVXh(:,(9+e4i))    = INFO.NNwts(3);
VXh(:,e4i)         = INFO.NNwts(3);
%}


% [[ SNP ]] UPDATE TARGET SNP OF HOLDOUT GROUP IN ONE-HOT DUAL-VAR-MX (DVXh)
%-----------------------------------------------------
vi = find(XLOCI.CHRPOS==CHRPOS);
DVXh(:,(vi*2+8))  = INFO.NNwts(3);
DVXh(:,(vi*2+9))  = INFO.NNwts(4);
DXh(:,vi*2-1)     = INFO.NNwts(3);
DXh(:,vi*2-0)     = INFO.NNwts(4);
PVXh(:,(9+vi))    = INFO.NNwts(4);
VXh(:,vi)         = INFO.NNwts(4);



% NET: GET ACTIVATIONS FOR TRAINING AND HOLDOUT
[CONFU_TR, PERF_TR, AREA_TR] = confusionmx(Yt, VXt, net );
[CONFU_HO, PERF_HO, AREA_HO] = confusionmx(Yh, VXh, net );

LOOPDATA.CONFU_TR(:,:,ij, 2) = CONFU_TR;
LOOPDATA.CONFU_HO(:,:,ij, 2) = CONFU_HO;
LOOPDATA.PERF_TR(ij,:,    2) = PERF_TR;
LOOPDATA.PERF_HO(ij,:,    2) = PERF_HO;
LOOPDATA.AREA_TR(:,:,ij,  2) = AREA_TR;
LOOPDATA.AREA_HO(:,:,ij,  2) = AREA_HO;
%------------------------------------------------------------------












%------------------------------------------------------------------
% USING ONE-HOT-TRAINING-MATRIX (DVXt, DXt)
%-------------------------------------------
% TRAINING APOES:  22 23 24 33 34 44
% HOLDOUT  APOE:   {24|34|44}
% HOLDOUT TARSNP:  ALT|ALT
%------------------------------------------------------------------
PHEt = TRPHE( sum( TRPHE.APOE==[24,34,44] ,2)>0 ,:);
PHEh = TEPHE( sum( TEPHE.APOE==[24,34,44] ,2)>0 ,:);

[PVXt, VXt, DVXt, DXt, Yt, YDt] = vxdx(XLOCI,XCASE,XCTRL,XUSNP,PHEt,INFO.NNwts);
[PVXh, VXh, DVXh, DXh, Yh, YDh] = vxdx(XLOCI,XCASE,XCTRL,XUSNP,PHEh,INFO.NNwts);


% [[ APOE ]] UPDATE APOE STATUS IN HOLDOUT GROUP FOR ONE-HOT DUAL-VAR-MX (DVXh)
%-----------------------------------------------------
%{
%--- SET E2{REF|REF}
e2i = find(XLOCI.CHRPOS==190045412079);
DVXh(:,(e2i*2+8))  = INFO.NNwts(1);
DVXh(:,(e2i*2+9))  = INFO.NNwts(1);
DXh(:,e2i*2-1)     = INFO.NNwts(1);
DXh(:,e2i*2-0)     = INFO.NNwts(1);
PVXh(:,(9+e2i))    = INFO.NNwts(1);
VXh(:,e2i)         = INFO.NNwts(1);

%--- SET E4{ALT|REF}
e4i = find(XLOCI.CHRPOS==190045411941);
DVXh(:,(e4i*2+8))  = INFO.NNwts(3);
DVXh(:,(e4i*2+9))  = INFO.NNwts(4);
DXh(:,e4i*2-1)     = INFO.NNwts(3);
DXh(:,e4i*2-0)     = INFO.NNwts(4);
PVXh(:,(9+e4i))    = INFO.NNwts(4);
VXh(:,e4i)         = INFO.NNwts(4);
%}


% [[ SNP ]] UPDATE TARGET SNP OF HOLDOUT GROUP IN ONE-HOT DUAL-VAR-MX (DVXh)
%-----------------------------------------------------
vi = find(XLOCI.CHRPOS==CHRPOS);
DVXh(:,(vi*2+8))  = INFO.NNwts(3);
DVXh(:,(vi*2+9))  = INFO.NNwts(4);
DXh(:,vi*2-1)     = INFO.NNwts(3);
DXh(:,vi*2-0)     = INFO.NNwts(4);
PVXh(:,(9+vi))    = INFO.NNwts(4);
VXh(:,vi)         = INFO.NNwts(4);



% NET: GET ACTIVATIONS FOR TRAINING AND HOLDOUT
[CONFU_TR, PERF_TR, AREA_TR] = confusionmx(Yt, VXt, net );
[CONFU_HO, PERF_HO, AREA_HO] = confusionmx(Yh, VXh, net );

LOOPDATA.CONFU_TR(:,:,ij, 3) = CONFU_TR;
LOOPDATA.CONFU_HO(:,:,ij, 3) = CONFU_HO;
LOOPDATA.PERF_TR(ij,:,    3) = PERF_TR;
LOOPDATA.PERF_HO(ij,:,    3) = PERF_HO;
LOOPDATA.AREA_TR(:,:,ij,  3) = AREA_TR;
LOOPDATA.AREA_HO(:,:,ij,  3) = AREA_HO;
%------------------------------------------------------------------










%------------------------------------------------------------------
% USING ONE-HOT-TRAINING-MATRIX (DVXt, DXt)
%-------------------------------------------
% TRAINING APOES:  22 23 24 33 34 44
% HOLDOUT  APOE:   {34|44}
% HOLDOUT TARSNP:  ALT|ALT   NNwts(3)|NNwts(4)
%------------------------------------------------------------------
PHEt = TRPHE(TRPHE.APOE >= 34 ,:);
PHEh = TEPHE(TEPHE.APOE >= 34 ,:);

[PVXt, VXt, DVXt, DXt, Yt, YDt] = vxdx(XLOCI,XCASE,XCTRL,XUSNP,PHEt,INFO.NNwts);
[PVXh, VXh, DVXh, DXh, Yh, YDh] = vxdx(XLOCI,XCASE,XCTRL,XUSNP,PHEh,INFO.NNwts);


% [[ APOE ]] UPDATE APOE STATUS IN HOLDOUT GROUP FOR ONE-HOT DUAL-VAR-MX (DVXh)
%-----------------------------------------------------
%{
%--- SET E2{REF|REF}
e2i = find(XLOCI.CHRPOS==190045412079);
DVXh(:,(e2i*2+8))  = INFO.NNwts(1);
DVXh(:,(e2i*2+9))  = INFO.NNwts(1);
DXh(:,e2i*2-1)     = INFO.NNwts(1);
DXh(:,e2i*2-0)     = INFO.NNwts(1);
PVXh(:,(9+e2i))    = INFO.NNwts(1);
VXh(:,e2i)         = INFO.NNwts(1);

%--- SET E4{ALT|REF}
e4i = find(XLOCI.CHRPOS==190045411941);
DVXh(:,(e4i*2+8))  = INFO.NNwts(3);
DVXh(:,(e4i*2+9))  = INFO.NNwts(4);
DXh(:,e4i*2-1)     = INFO.NNwts(3);
DXh(:,e4i*2-0)     = INFO.NNwts(4);
PVXh(:,(9+e4i))    = INFO.NNwts(4);
VXh(:,e4i)         = INFO.NNwts(4);
%}


% [[ SNP ]] UPDATE TARGET SNP OF HOLDOUT GROUP IN ONE-HOT DUAL-VAR-MX (DVXh)
%-----------------------------------------------------
vi = find(XLOCI.CHRPOS==CHRPOS);
DVXh(:,(vi*2+8))  = INFO.NNwts(3);
DVXh(:,(vi*2+9))  = INFO.NNwts(4);
DXh(:,vi*2-1)     = INFO.NNwts(3);
DXh(:,vi*2-0)     = INFO.NNwts(4);
PVXh(:,(9+vi))    = INFO.NNwts(4);
VXh(:,vi)         = INFO.NNwts(4);



% NET: GET ACTIVATIONS FOR TRAINING AND HOLDOUT
[CONFU_TR, PERF_TR, AREA_TR] = confusionmx(Yt, VXt, net );
[CONFU_HO, PERF_HO, AREA_HO] = confusionmx(Yh, VXh, net ,1);

LOOPDATA.CONFU_TR(:,:,ij, 4) = CONFU_TR;
LOOPDATA.CONFU_HO(:,:,ij, 4) = CONFU_HO;
LOOPDATA.PERF_TR(ij,:,    4) = PERF_TR;
LOOPDATA.PERF_HO(ij,:,    4) = PERF_HO;
LOOPDATA.AREA_TR(:,:,ij,  4) = AREA_TR;
LOOPDATA.AREA_HO(:,:,ij,  4) = AREA_HO;
%------------------------------------------------------------------











%==========================================================================
%==========================================================================
% 
end
%
%==========================================================================
%==========================================================================
% clc; clearvars -except P ADSP INFO PHEN LOCI CASE CTRL USNP LOOPDATA kk ij...
% VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL SNPi XLOCI...
% XCASE XCTRL XUSNP TRPHE TEPHE PVXt DVXt PVXh DVXh VXt VXh Yt Yh...
% YDt YDh DXt DXh netq netd net 
%----------------------------------------------------------------------
clearvars -except P ADSP INFO SNPTAB PHEN LOCI CASE CTRL USNP LOOPDATA kk ij...
XLOCI XCASE XCTRL XUSNP TRPHE TEPHE net









%% SAVE LOOP DATA


GENE   = ADSP.SNP.GENE{kk};
CHRPOS = ADSP.SNP.CHRPOS(kk);


% SAVE MAT FILE
%------------------------------------------%
P.SAVEPATH = [P.Output_Path P.f GENE '_' num2str(CHRPOS) '.mat'];
save(P.SAVEPATH,'LOOPDATA','P','INFO');
disp('File saved...'); disp(P.SAVEPATH)
%------------------------------------------%



P.FINAME1 = '_SNP44_APOE_22_23_24_33_34_44';
P.FINAME2 = '_SNP44_APOE_33';
P.FINAME3 = '_SNP44_APOE_24_34_44';
P.FINAME4 = '_SNP44_APOE_34_44';



%------------------------------------------------------------------
% TRAINING APOE:  22 23 24 33 34 44
% HOLDOUT  APOE2: -1 (REF/REF)
% HOLDOUT  APOE4:  2 (REF/ALT)
% HOLDOUT  VARIS:  2 (REF/ALT)
%------------------------------------------------------------------
fh01 = figure('Units','normalized','OuterPosition',[.01 .06 .55 .37],'Color','w');
ax1 = axes('Position',[.07 .16 .61 .80],'Color','none');
ax2 = axes('Position',[.78 .16 .21 .80],'Color','none');

%AREA:[X, CASE_OK, CTRL_OK, CASE_BAD, CTRL_BAD] maybe?
AREAS = mean(LOOPDATA.AREA_HO(:,:,:,1)  ,3);

%PERF:[MU_CASE, MU_CTRL, MU_ALL, HIMU_CASE, HIMU_CTRL, HIMU_ALL]
PERFS = mean(LOOPDATA.PERF_HO(:,:,1)) .* 100;

%---PLOT AREAGRAM ---------------------
    axes(ax1) 
area(AREAS(:,1),AREAS(:,2:5))
    legend({'Case Miss','Ctrl Miss',...
            'Case Hit','Ctrl Hit'},...
            'Location','best');
    ylabel('Count')
    ax1.FontSize=16;
    %title('Classifier Performance Area');

%---BAR GRAPH ------------------------------------------
    axes(ax2)
bar(([ PERFS(4) , PERFS(5) , PERFS(6) ]),.5, 'FaceColor',[.31 .31 .31]); 
    hold on; 
bar(([ PERFS(1) , PERFS(2) , PERFS(3) ]),.20,'FaceColor',[.95 .85 .50]);
    grid on; ylabel('Pct. Correct')
    legend({'Top 25%','All'},'Location','Northwest','NumColumns',2)
    ax2.YLim = [0 116]; 
    ax2.YTick = [0 25 50 75 100]; 
    ax2.XTickLabels = {'CASE','CTRL','ALL'};
    ax2.XTickLabelRotation = 33;
    ax2.FontSize=16;
    %title('Performance Summary')
%----------------------------------------------------------------------
pause(1)
set(gcf, 'PaperPositionMode', 'auto');
%dt=char(datetime(datetime,'Format','yyyy-MM-dd-HH-mm-ss'));
saveas(gcf, [P.Output_Path P.f GENE '_' num2str(CHRPOS) P.FINAME1 '.png']);
pause(1)
%---------------------------------------------------------------------- 




%------------------------------------------------------------------
% TRAINING APOE:  22 23 24 33 34 44
% HOLDOUT  APOE2: -1 (REF/REF)
% HOLDOUT  APOE4:  2 (REF/ALT)
% HOLDOUT  VARIS:  3 (ALT/ALT)
%------------------------------------------------------------------
fh02 = figure('Units','normalized','OuterPosition',[.01 .06 .55 .37],'Color','w');
ax1 = axes('Position',[.07 .16 .61 .80],'Color','none');
ax2 = axes('Position',[.78 .16 .21 .80],'Color','none');

%AREA:[X, CASE_OK, CTRL_OK, CASE_BAD, CTRL_BAD] maybe?
AREAS = mean(LOOPDATA.AREA_HO(:,:,:,2)  ,3);

%PERF:[MU_CASE, MU_CTRL, MU_ALL, HIMU_CASE, HIMU_CTRL, HIMU_ALL]
PERFS = mean(LOOPDATA.PERF_HO(:,:,2)) .* 100;

%---PLOT AREAGRAM ---------------------
    axes(ax1) 
area(AREAS(:,1),AREAS(:,2:5))
    legend({'Case Miss','Ctrl Miss',...
            'Case Hit','Ctrl Hit'},...
            'Location','best');
    ylabel('Count')
    ax1.FontSize=16;
    %title('Classifier Performance Area');

%---BAR GRAPH ------------------------------------------
    axes(ax2)
bar(([ PERFS(4) , PERFS(5) , PERFS(6) ]),.5, 'FaceColor',[.31 .31 .31]); 
    hold on; 
bar(([ PERFS(1) , PERFS(2) , PERFS(3) ]),.20,'FaceColor',[.95 .85 .50]);
    grid on; ylabel('Pct. Correct')
    legend({'Top 25%','All'},'Location','Northwest','NumColumns',2)
    ax2.YLim = [0 116]; 
    ax2.YTick = [0 25 50 75 100]; 
    ax2.XTickLabels = {'CASE','CTRL','ALL'};
    ax2.XTickLabelRotation = 33;
    ax2.FontSize=16;
    %title('Performance Summary')
%----------------------------------------------------------------------
pause(1)
set(gcf, 'PaperPositionMode', 'auto');
%dt=char(datetime(datetime,'Format','yyyy-MM-dd-HH-mm-ss'));
saveas(gcf, [P.Output_Path P.f GENE '_' num2str(CHRPOS) P.FINAME2 '.png']);
pause(1)
%---------------------------------------------------------------------- 


%------------------------------------------------------------------
% TRAINING APOE:  22 23 24 33 34 44
% HOLDOUT  APOE2: -1 (REF/REF)
% HOLDOUT  APOE4:  3 (ALT/ALT)
% HOLDOUT  VARIS:  2 (REF/ALT)
%------------------------------------------------------------------
fh02 = figure('Units','normalized','OuterPosition',[.01 .06 .55 .37],'Color','w');
ax1 = axes('Position',[.07 .16 .61 .80],'Color','none');
ax2 = axes('Position',[.78 .16 .21 .80],'Color','none');

%AREA:[X, CASE_OK, CTRL_OK, CASE_BAD, CTRL_BAD] maybe?
AREAS = mean(LOOPDATA.AREA_HO(:,:,:,3)  ,3);

%PERF:[MU_CASE, MU_CTRL, MU_ALL, HIMU_CASE, HIMU_CTRL, HIMU_ALL]
PERFS = mean(LOOPDATA.PERF_HO(:,:,3)) .* 100;

%---PLOT AREAGRAM ---------------------
    axes(ax1) 
area(AREAS(:,1),AREAS(:,2:5))
    legend({'Case Miss','Ctrl Miss',...
            'Case Hit','Ctrl Hit'},...
            'Location','best');
    ylabel('Count')
    ax1.FontSize=16;
    %title('Classifier Performance Area');

%---BAR GRAPH ------------------------------------------
    axes(ax2)
bar(([ PERFS(4) , PERFS(5) , PERFS(6) ]),.5, 'FaceColor',[.31 .31 .31]); 
    hold on; 
bar(([ PERFS(1) , PERFS(2) , PERFS(3) ]),.20,'FaceColor',[.95 .85 .50]);
    grid on; ylabel('Pct. Correct')
    legend({'Top 25%','All'},'Location','Northwest','NumColumns',2)
    ax2.YLim = [0 116]; 
    ax2.YTick = [0 25 50 75 100]; 
    ax2.XTickLabels = {'CASE','CTRL','ALL'};
    ax2.XTickLabelRotation = 33;
    ax2.FontSize=16;
    %title('Performance Summary')
%----------------------------------------------------------------------
pause(1)
set(gcf, 'PaperPositionMode', 'auto');
%dt=char(datetime(datetime,'Format','yyyy-MM-dd-HH-mm-ss'));
saveas(gcf, [P.Output_Path P.f GENE '_' num2str(CHRPOS) P.FINAME3 '.png']);
pause(1)
%---------------------------------------------------------------------- 





%------------------------------------------------------------------
% TRAINING APOE:  22 23 24 33 34 44
% HOLDOUT  APOE2: -1 (REF/REF)
% HOLDOUT  APOE4:  3 (ALT/ALT)
% HOLDOUT  VARIS:  3 (ALT/ALT)
%------------------------------------------------------------------
fh02 = figure('Units','normalized','OuterPosition',[.01 .06 .55 .37],'Color','w');
ax1 = axes('Position',[.07 .16 .61 .80],'Color','none');
ax2 = axes('Position',[.78 .16 .21 .80],'Color','none');

%AREA:[X, CASE_OK, CTRL_OK, CASE_BAD, CTRL_BAD] maybe?
AREAS = mean(LOOPDATA.AREA_HO(:,:,:,4)  ,3);

%PERF:[MU_CASE, MU_CTRL, MU_ALL, HIMU_CASE, HIMU_CTRL, HIMU_ALL]
PERFS = mean(LOOPDATA.PERF_HO(:,:,4)) .* 100;

%---PLOT AREAGRAM ---------------------
    axes(ax1) 
area(AREAS(:,1),AREAS(:,2:5))
    legend({'Case Miss','Ctrl Miss',...
            'Case Hit','Ctrl Hit'},...
            'Location','best');
    ylabel('Count')
    ax1.FontSize=16;
    %title('Classifier Performance Area');

%---BAR GRAPH ------------------------------------------
    axes(ax2)
bar(([ PERFS(4) , PERFS(5) , PERFS(6) ]),.5, 'FaceColor',[.31 .31 .31]); 
    hold on; 
bar(([ PERFS(1) , PERFS(2) , PERFS(3) ]),.20,'FaceColor',[.95 .85 .50]);
    grid on; ylabel('Pct. Correct')
    legend({'Top 25%','All'},'Location','Northwest','NumColumns',2)
    ax2.YLim = [0 116]; 
    ax2.YTick = [0 25 50 75 100]; 
    ax2.XTickLabels = {'CASE','CTRL','ALL'};
    ax2.XTickLabelRotation = 33;
    ax2.FontSize=16;
    %title('Performance Summary')
%----------------------------------------------------------------------
pause(1)
set(gcf, 'PaperPositionMode', 'auto');
%dt=char(datetime(datetime,'Format','yyyy-MM-dd-HH-mm-ss'));
saveas(gcf, [P.Output_Path P.f GENE '_' num2str(CHRPOS) P.FINAME4 '.png']);
pause(1)
%---------------------------------------------------------------------- 











%##########################################################################
%##########################################################################
% 
end
%
%##########################################################################
%##########################################################################