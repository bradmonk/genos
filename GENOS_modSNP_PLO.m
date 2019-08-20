%% GENOS_FISHP
%{
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
%}
%==========================================================================
%% STEP-1: LOAD THE DATASET
%==========================================================================
close all; clear; clc; rng('shuffle');
P.home = fileparts(which('GENOS.m')); cd(P.home);
P.funs = [P.home filesep 'genos_functions'];
P.mfuns = [P.funs filesep 'genos_main_functions'];
P.other = [P.home filesep 'genos_other'];
P.data = [P.home filesep 'genos_data'];
P.data = 'F:\GENOSDATA';
addpath(join(string(struct2cell(P)),pathsep,1))
cd(P.home); P.f = filesep;


ADSP = load('GENOSDATA.mat');


P.mainmatfile = which('GENOSDATA.mat');
disp('dataset loaded')
clearvars -except P ADSP INFO
%==========================================================================
%% SET RUN OPTIONS & PATHS FOR DATA IMPORT/EXPORT
%==========================================================================
clc; clearvars -except P ADSP INFO



P.nSNPsets = 4000;
P.TargetSNPs = 1;
P.TopSNPs = 49;
P.Nloops = 5;
P.Nsnps = P.TopSNPs + P.TargetSNPs;
%P.Ngenes = P.NgeneEnd - P.NgeneStart + 1;
P.Output_Path = [P.data P.f 'modSNP' P.f 'RUN10_PLO_NOAPOE' P.f 'MAT'];
P.TargetSNP_Path = [P.data P.f 'CSVSNP' P.f 'PLO.csv'];


clearvars -except P ADSP INFO
%==========================================================================
%% SET PATH TO (APOE DEPENDENT) SUBSAMPLE FOLDERS AND GET MAT PATHS
%==========================================================================
clc; clearvars -except P ADSP INFO



P.APOExx_path = [P.data P.f 'SUBSETS' P.f 'APOE_22_23_24_33_34_44' P.f 'APOE_22_23_24_33_34_44_FISHP'];
P.APOExxFILES.sets = what(P.APOExx_path);
P.APOExxN = numel(P.APOExxFILES.sets.mat);


P.APOE4x_path = [P.data P.f 'SUBSETS' P.f 'APOE_34_44' P.f 'APOE_34_44_FISHP'];
P.APOE4xFILES.sets = what(P.APOE4x_path);
P.APOE4xN = numel(P.APOE4xFILES.sets.mat);



fprintf('\n APOExx has N mat files: %0.f \n',P.APOExxN)
fprintf('\n APOE4x has N mat files: %0.f \n',P.APOE4xN)




clearvars -except P ADSP INFO
%==========================================================================
%% IMPORT TOP VARIANTS FROM EXCEL SHEET
%==========================================================================
clc; clearvars -except P ADSP INFO


ops  = detectImportOptions([P.TargetSNP_Path]);
SNPTAB  = readtable([P.TargetSNP_Path],ops);
SNPTAB.CHRPOS  = uint64(SNPTAB.CHRPOS);
ADSP.SNP = SNPTAB;



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
%% LOOP OVER EACH GENE-CHRPOS TO PERTURB
%##########################################################################
%
% 
%
for kk = 1:P.nSNPsets
%
%
% 
%##########################################################################
close all; clc;
clearvars -except P ADSP INFO SNPTAB PHEN LOCI CASE CTRL USNP kk ij


% NEURAL NETWORKS CONFUSION STATS
%----------------------------------------
LOOPDATA.TR_ALL_STATS = zeros(9,P.Nloops);
LOOPDATA.HO_ALL_STATS = zeros(9,P.Nloops);
LOOPDATA.TR_TOP_STATS = zeros(9,P.Nloops);
LOOPDATA.HO_TOP_STATS = zeros(9,P.Nloops);



% CASES NEURAL NETWORKS
%----------------------------------------
LOOPDATA.CATRMEAN = zeros(P.Nloops);
LOOPDATA.CATRLOMU = zeros(P.Nloops);
LOOPDATA.CATRHIMU = zeros(P.Nloops);
LOOPDATA.CATRLOPO = zeros(P.Nloops);
LOOPDATA.CATRHIPO = zeros(P.Nloops);

LOOPDATA.CAHOMEAN = zeros(P.Nloops);
LOOPDATA.CAHOLOMU = zeros(P.Nloops);
LOOPDATA.CAHOHIMU = zeros(P.Nloops);
LOOPDATA.CAHOLOPO = zeros(P.Nloops);
LOOPDATA.CAHOHIPO = zeros(P.Nloops);



% TRAINED MACHINE LEARNER
%----------------------------------------
LOOPDATA.NETS = cell(P.Nloops);



% NEURAL NETWORK & CONFUSION MATRIX
%----------------------------------------
LOOPDATA.CONFU_TR = nan(4,4,P.Nloops,4);
LOOPDATA.CONFU_HO = nan(4,4,P.Nloops,4);
LOOPDATA.PERF_TR  = nan(P.Nloops,6,4);
LOOPDATA.PERF_HO  = nan(P.Nloops,6,4);
LOOPDATA.AREA_TR  = nan(50,5,P.Nloops,4);
LOOPDATA.AREA_HO  = nan(50,5,P.Nloops,4);



% TARGET GENE & CHRPOS
%----------------------------------------
LOOPDATA.GENECHRPOS = cell(P.Nloops,2);


% TRAINING MATRIX INFO
%----------------------------------------
LOOPDATA.XLOCI{P.Nloops} = LOCI(1:P.Nsnps,:);
LOOPDATA.PVXt       = nan(7000,(P.TargetSNPs+9),P.Nloops);
LOOPDATA.PVXh       = nan(4000,(P.TargetSNPs+9),P.Nloops);
LOOPDATA.NNNATt     = nan(7000,(P.TargetSNPs+9),P.Nloops);
LOOPDATA.NNNATh     = nan(4000,(P.TargetSNPs+9),P.Nloops);
LOOPDATA.NNREFt     = nan(7000,(P.TargetSNPs+9),P.Nloops);
LOOPDATA.NNREFh     = nan(4000,(P.TargetSNPs+9),P.Nloops);
LOOPDATA.NNALTt     = nan(7000,(P.TargetSNPs+9),P.Nloops);
LOOPDATA.NNALTh     = nan(4000,(P.TargetSNPs+9),P.Nloops);



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
    % REMOVE TOMM40 FROM CONTENTION UNLESS A TARGET SNP
    %======================================================================

    x = strcmp(VLOCI.GENE,"TOMM40");
    VLOCI.TRFISHP( x ) = 1;
    VLOCI.TEFISHP( x ) = 1;
    VLOCI.FISHP( x )   = 1;



    %======================================================================
    % REMOVE APOE FROM THE TRAINING MATRIX UNLESS A TARGET SNP
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

    
    x = strcmp(VLOCI.GENE,"APOE");
    VLOCI.TRFISHP( x ) = 1;
    VLOCI.TEFISHP( x ) = 1;
    VLOCI.FISHP( x )   = 1;
    

    %==========================================================================
    % IDENTIFY TARGET SNPS AND MAKE THEIR P-VALUE = 0
    %==========================================================================
    %keyboard

    s = (kk*P.TargetSNPs-(P.TargetSNPs-1)):(kk*P.TargetSNPs);

    GENE   = string(ADSP.SNP.GENE(s));
    CHRPOS = ADSP.SNP.CHRPOS(s);

    % MARK ROW OF TARGET SNPs WITH LOGICAL 1
    g = sum(VLOCI.CHRPOS == (CHRPOS') ,2)>0;
    
    if sum(g) ~= P.TargetSNPs; keyboard; end


    VLOCI.TRFISHP( g ) = 0;

    LOOPDATA.GENECHRPOS{ij,1} = GENE;
    LOOPDATA.GENECHRPOS{ij,2} = CHRPOS;



    %======================================================================
    % MAKE 50 RANDOM SNPS HAVE VERY LOW P-VALUES
    %======================================================================
    %{
    i = randperm(size(LOCI,1),50);
    VLOCI.TRFISHP( i ) = 1e-100;
    %}



    %======================================================================
    % SORT VARIANTS BY TRFISHP
    %======================================================================
    [~,j]  = sort(VLOCI.TRFISHP);
    VLOCI  = VLOCI(j,:);
    VCASE  = VCASE(j);
    VCTRL  = VCTRL(j);
    VUSNP  = VUSNP(j);


    clc;
    fprintf('\nTarget SNP Genes: \n\n');
    disp(GENE')
    disp(VLOCI(1:12,[1:8 31 32 41 42 43 44])); 
    pause(1);







    %======================================================================
    %      MAKE  RECTANGLE  NEURAL NET  VARIANT MATRIX
    %======================================================================
    SNPi = 1:P.Nsnps;


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



    % Save some data for the record
    LOOPDATA.XLOCI{ij} = XLOCI;
    LOOPDATA.PVXt(1:(size(PVXt,1)),1:(size(PVXt,2)),ij)  = PVXt;
    LOOPDATA.PVXh(1:(size(PVXh,1)),1:(size(PVXh,2)),ij)  = PVXh;



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



% [CONFU_TR, PERF_TR, AREA_TR] = confusionmx(Yt, VXt, net);
% [CONFU_HO, PERF_HO, AREA_HO] = confusionmx(Yh, VXh, net,1);
% [CONFU_TR, PERF_TR, PCUT_TR, AREA_TR] = nnperform(Yt', DXt', PVXt, net, 'cut',.5,'plt',0);
% [CONFU_HO, PERF_HO, PCUT_HO, AREA_HO] = nnperform(Yh', DXh', PVXh, net, 'cut',.5,'plt',0);




%======================================================================
%% SERIAL SNP TEST
%======================================================================
clearvars -except P ADSP INFO SNPTAB PHEN LOCI CASE CTRL USNP LOOPDATA kk ij...
XLOCI XCASE XCTRL XUSNP TRPHE TEPHE net



NNNATt = nan(size(TRPHE,1),P.TargetSNPs);
NNNATh = nan(size(TEPHE,1),P.TargetSNPs);
NNREFt = nan(size(TRPHE,1),P.TargetSNPs);
NNREFh = nan(size(TEPHE,1),P.TargetSNPs);
NNALTt = nan(size(TRPHE,1),P.TargetSNPs);
NNALTh = nan(size(TEPHE,1),P.TargetSNPs);


% ITERATE OVER EACH SNP IN THE NN MATRIX
%-------------------------------------------
for vi = 1:P.TargetSNPs

% PVXt(SRR,AD,COH,AGE,APOE,BRAGEz,BRAAK,AGEz,BRAGEz ...)
[PVXt, VXt, DVXt, DXt, Yt, YDt] = vxdx(XLOCI,XCASE,XCTRL,XUSNP,TRPHE,INFO.NNwts);
[PVXh, VXh, DVXh, DXh, Yh, YDh] = vxdx(XLOCI,XCASE,XCTRL,XUSNP,TEPHE,INFO.NNwts);


NNNATt(:,vi)  = net(VXt')';
NNNATh(:,vi)  = net(VXh')';

VXt(:,vi)     = INFO.NNwts(1);
VXh(:,vi)     = INFO.NNwts(1);

NNREFt(:,vi)  = net(VXt')';
NNREFh(:,vi)  = net(VXh')';

VXt(:,vi)     = INFO.NNwts(4);
VXh(:,vi)     = INFO.NNwts(4);

NNALTt(:,vi)  = net(VXt')';
NNALTh(:,vi)  = net(VXh')';

end
%-------------------------------------------


NNNATt = [PVXt(:,1:9) , (NNNATt-.5)];
NNNATh = [PVXh(:,1:9) , (NNNATh-.5)];
NNREFt = [PVXt(:,1:9) , (NNREFt-.5)];
NNREFh = [PVXh(:,1:9) , (NNREFh-.5)];
NNALTt = [PVXt(:,1:9) , (NNALTt-.5)];
NNALTh = [PVXh(:,1:9) , (NNALTh-.5)];


clearvars -except P ADSP INFO SNPTAB PHEN LOCI CASE CTRL USNP LOOPDATA kk ij...
XLOCI XCASE XCTRL XUSNP TRPHE TEPHE net NNNATt NNNATh NNREFt NNREFh NNALTt NNALTh PVXt PVXh


%% SAVE LOOP DATA

LOOPDATA.XLOCI{ij}  = XLOCI;

LOOPDATA.PVXt(1:size(PVXt,1),:,ij)     = PVXt;
LOOPDATA.PVXh(1:size(PVXh,1),:,ij)     = PVXh;

LOOPDATA.NNNATt(1:size(NNNATt,1),:,ij) = NNNATt;
LOOPDATA.NNNATt(1:size(NNNATh,1),:,ij) = NNNATh;
LOOPDATA.NNREFt(1:size(NNREFt,1),:,ij) = NNREFt;
LOOPDATA.NNREFh(1:size(NNREFh,1),:,ij) = NNREFh;
LOOPDATA.NNALTt(1:size(NNALTt,1),:,ij) = NNALTt;
LOOPDATA.NNALTh(1:size(NNALTh,1),:,ij) = NNALTh;



%==========================================================================
%==========================================================================
% 
end
%
%==========================================================================
%==========================================================================
clearvars -except P ADSP INFO SNPTAB PHEN LOCI CASE CTRL USNP LOOPDATA kk ij...
XLOCI XCASE XCTRL XUSNP TRPHE TEPHE net NNREFt NNREFh NNALTt NNALTh PVXt PVXh





%% SAVE LOOP DATA


%------------------------------------------%
%dt=char(datetime(datetime,'Format','yyyy-MM-dd-HH-mm-ss'));
s1 = sprintf('%04d',kk*10-9);
s2 = sprintf('%04d',kk*10);
P.SAVEPATH = [P.Output_Path P.f 'SNPS_' s1 '_' s2 '.mat'];
save(P.SAVEPATH,'LOOPDATA','P','INFO');
disp('File saved...'); disp(P.SAVEPATH)
%------------------------------------------%


%##########################################################################
%##########################################################################
% 
end
%
%##########################################################################
%##########################################################################