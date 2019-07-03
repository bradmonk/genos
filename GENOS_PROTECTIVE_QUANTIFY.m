%% GENOS: 
%{
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
% APOE4 CHRPOS: 190045411941
% APOE2 CHRPOS: 190045412079
%
% APOE  190045411941  190045412079
% e22:      -1            5
% e23:      -1            2
% e24:       2            2
% e33:      -1           -1
% e34:       2           -1
% e44:       5           -1
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
%}
%--------------------------------------------------------------------------

%==========================================================================
%% STEP-1: LOAD THE DATASET
%==========================================================================
close all; clear; clc; rng('shuffle');
P.home = fileparts(which('GENOS.m')); cd(P.home);
P.P1 = [P.home filesep 'genos_functions'];
P.P3 = [P.P1 filesep 'genos_main_functions'];
P.P4 = [P.home filesep 'genos_other'];
addpath(join(string(struct2cell(P)),pathsep,1))
cd(P.home)


ADSP = load('GENOSDATAFINAL.mat');


P.mainmatfile = which('GENOSDATAFINAL.mat');
disp('dataset loaded')
clearvars -except P ADSP INFO



%==========================================================================
%% APOE related options
%==========================================================================
clc; clearvars -except P ADSP INFO

P.basedir = '/Users/bradleymonk/Documents/MATLAB/GIT/genomics/genos/genos_data/APOE';
P.f = filesep;


P.APOES = '22_23_24_33_34_44'; INFO.APOE = [22 23 24 33 34 44];
P.importdir = [P.basedir P.f 'APOE_22_23_24_33_34_44' P.f 'APOE_22_23_24_33_34_44_FISHP'];

% P.APOES = '22_23_24_34_44'; INFO.APOE = [22 23 24 34 44];
% P.importdir = [P.basedir P.f 'APOE_22_23_24_34_44' P.f 'APOE_22_23_24_34_44_FISHP'];

% P.APOES = '33'; INFO.APOE = [33];
% P.importdir = [P.basedir P.f 'APOE_33' P.f 'APOE_33_FISHP'];

% P.APOES = '34_44'; INFO.APOE = [34 44];
% P.importdir = [P.basedir P.f 'APOE_34_44' P.f 'APOE_34_44_FISHP'];

% P.APOES = '33_44'; INFO.APOE = [34 44];
% P.importdir = [P.basedir P.f 'APOE_34_44' P.f 'APOE_34_44_FISHP'];


clearvars -except P ADSP INFO





%==========================================================================
%% (GENOX) PERTURB TEST OPTIONS & FOLDERS 
%==========================================================================
clc; clearvars -except P ADSP INFO


% SET PERTURB TEST OPTIONS
%------------------------------------------------------
P.f = filesep;
P.NGenes = 100;
P.Nloops = 10;
P.FileStart = 1;
P.Nvars = 200;
P.windowSize = 50;
P.Ndots = 5;
P.Lo2Hi = true;
P.RemoveGenesByName = false;
P.RemoveBadGenes = false;



% SET DEFAULT FOLDERS 
%------------------------------------------------------
% P.genox.DIRroot = [P.home P.f 'genos_data' P.f 'GENOX'];
% P.genox.DIRmat  = [P.genox.DIRroot P.f 'PLO_MAT'];
% P.genox.DIRimg  = [P.genox.DIRroot P.f 'PLO_IMG'];
% P.genox.DIRstatsmat = [P.genox.DIRroot P.f 'PLO_STATS_MAT'];
% P.genox.DIRstatsimg = [P.genox.DIRroot P.f 'PLO_STATS_IMG'];
% 
% 
% P.genox.DIRroot = [P.home P.f 'genos_data' P.f 'GENOX'];
% P.genox.DIRmat  = [P.genox.DIRroot P.f 'ORLO_MAT'];
% P.genox.DIRimg  = [P.genox.DIRroot P.f 'ORLO_IMG'];
% P.genox.DIRstatsmat = [P.genox.DIRroot P.f 'ORLO_STATS_MAT'];
% P.genox.DIRstatsimg = [P.genox.DIRroot P.f 'ORLO_STATS_IMG'];
% 
% 
% P.genox.DIRroot = [P.home P.f 'genos_data' P.f 'GENOX'];
% P.genox.DIRmat  = [P.genox.DIRroot P.f 'ORHI_MAT'];
% P.genox.DIRimg  = [P.genox.DIRroot P.f 'ORHI_IMG'];
% P.genox.DIRstatsmat = [P.genox.DIRroot P.f 'ORHI_STATS_MAT'];
% P.genox.DIRstatsimg = [P.genox.DIRroot P.f 'ORHI_STATS_IMG'];
% 
% 
% P.genox.DIRroot = [P.home P.f 'genos_data' P.f 'GENOX'];
% P.genox.DIRmat  = [P.genox.DIRroot P.f 'SYN_MAT'];
% P.genox.DIRimg  = [P.genox.DIRroot P.f 'SYN_IMG'];
% P.genox.DIRstatsmat = [P.genox.DIRroot P.f 'SYN_STATS_MAT'];
% P.genox.DIRstatsimg = [P.genox.DIRroot P.f 'SYN_STATS_IMG'];

P.genox.DIRroot = [P.home P.f 'genos_data' P.f 'GENOX'];
P.genox.DIRmat  = [P.genox.DIRroot P.f 'PRO/PRO_MAT'];
P.genox.DIRimg  = [P.genox.DIRroot P.f 'PRO/PRO_MAT'];
P.genox.DIRstatsmat = [P.genox.DIRroot P.f 'PRO/PRO_STATS_MAT'];
P.genox.DIRstatsimg = [P.genox.DIRroot P.f 'PRO/PRO_STATS_IMG'];




% IMPORT MAT FILES
%------------------------------------------------------
P.PERT.w = what(P.genox.DIRmat);
P.PERT.finfo = dir(P.PERT.w.path);
P.PERT.finames = {P.PERT.finfo.name};
c=~cellfun(@isempty,regexp(P.PERT.finames,'((\S)+(\.mat+))'));
P.PERT.finames = string(P.PERT.finames(c)');
P.PERT.folder = P.PERT.finfo.folder;
P.PERT.fipaths = fullfile(P.PERT.folder,P.PERT.finames);
disp(P.PERT.fipaths); disp(P.PERT.finames);


ij=2;
PERT = load(P.PERT.fipaths{ij});

clearvars -except P ADSP INFO PERT


%==========================================================================
%% GET PERTURB STATS
%==========================================================================
clc; clearvars -except P ADSP INFO PERT


for ij = 1:numel(P.PERT.fipaths)
clc; clearvars -except P ADSP INFO PER ij
disp(ij);
disp(P.PERT.finames{ij})


PERT = load(P.PERT.fipaths{ij});

PER(1).GENE   = PERT.INFO.TOPLOCI{1,1}.GENE(1);
PER(1).CHRPOS = PERT.INFO.TOPLOCI{1,1}.CHRPOS(1);


% NAT/NAT PERFORMANCE
%-----------------------------------
PER(1).area = mean(PERT.LOOPDATA.AREA_HO(:,:,1:10,1) ,3);
PER(1).caseMu = midgravity(PER(1).area(:,1), PER(1).area(:,2), PER(1).area(:,4));
PER(1).ctrlMu = midgravity(PER(1).area(:,1), PER(1).area(:,3), PER(1).area(:,5));
PER(1).cacoZd = 0;


% A = mean(PERT.LOOPDATA.AREA_HO(:,:,1:10,1) ,3);
% close all
% area(A(:,1),A(:,2:end))


% REF/REF PERFORMANCE
%-----------------------------------
PER(2).area = mean(PERT.LOOPDATA.AREA_HO(:,:,1:10,2) ,3);
PER(2).caseMu = midgravity(PER(2).area(:,1), PER(2).area(:,2), PER(2).area(:,4));
PER(2).ctrlMu = midgravity(PER(2).area(:,1), PER(2).area(:,3), PER(2).area(:,5));
PER(2).cacoZd = (abs(PER(2).caseMu - PER(1).caseMu)+...
                 abs(PER(2).ctrlMu - PER(1).ctrlMu))/2;


% REF/ALT PERFORMANCE
%-----------------------------------
PER(3).area = mean(PERT.LOOPDATA.AREA_HO(:,:,1:10,3) ,3);
PER(3).caseMu = midgravity(PER(3).area(:,1), PER(3).area(:,2), PER(3).area(:,4));
PER(3).ctrlMu = midgravity(PER(3).area(:,1), PER(3).area(:,3), PER(3).area(:,5));
PER(3).cacoZd = (abs(PER(3).caseMu - PER(1).caseMu)+...
                 abs(PER(3).ctrlMu - PER(1).ctrlMu))/2;





% ALT/ALT PERFORMANCE
%-----------------------------------
PER(4).area = mean(PERT.LOOPDATA.AREA_HO(:,:,1:10,4) ,3);
PER(4).caseMu = midgravity(PER(4).area(:,1), PER(4).area(:,2), PER(4).area(:,4));
PER(4).ctrlMu = midgravity(PER(4).area(:,1), PER(4).area(:,3), PER(4).area(:,5));
PER(4).cacoZd = (abs(PER(4).caseMu - PER(1).caseMu)+...
                 abs(PER(4).ctrlMu - PER(1).ctrlMu))/2;



%------------------------------------------%
P.SAVEPATH = [P.genox.DIRstatsmat P.f ...
              PER(1).GENE{1} '_' num2str(PER(1).CHRPOS) '.mat'];
% save(P.SAVEPATH,'P','PERT','PER');
save(P.SAVEPATH,'PER');
disp('File saved...'); disp(P.SAVEPATH)
%------------------------------------------%

end


clc; clearvars -except P ADSP INFO PERT PER
%==========================================================================
%% PLOT CONFUSION PERFORMANCE HISTOGRAMS
%==========================================================================


% IMPORT PEROX FILES
%-----------------------------------
clc; clearvars -except P ADSP INFO PERT PER kk

P.PERX.w = what(P.genox.DIRstatsmat);
P.PERX.finfo = dir(P.PERX.w.path);
P.PERX.finames = {P.PERX.finfo.name};
c=~cellfun(@isempty,regexp(P.PERX.finames,'((\S)+(\.mat+))'));
P.PERX.finames = string(P.PERX.finames(c)');
P.PERX.folder = P.PERX.finfo.folder;
P.PERX.fipaths = fullfile(P.PERX.folder,P.PERX.finames);
disp(P.PERX.fipaths); disp(P.PERX.finames);


cm=cmbyr(15); cmap = cm([12,5,14,2,15,1],:);
set(groot,'defaultAxesColorOrder',cmap)
% edges = (-.49:.02:.49);




for ij = 1:numel(P.PERX.fipaths)
clc; clearvars -except P ADSP INFO PERT PER ij
disp(ij);



load(P.PERX.fipaths{ij});


f = char(P.PERX.finames(ij));
gnam = f(1:end-4);



close all
fh01 = figure('Units','normalized','OuterPosition',[.01 .06 .97 .88],'Color','w');

ax10 = axes('Position',[.05 .59 .31 .37],'Color','none');
ax11 = axes('Position',[.42 .59 .06 .37],'Color','none');

ax20 = axes('Position',[.54 .59 .31 .37],'Color','none');
ax21 = axes('Position',[.91 .59 .06 .37],'Color','none');

ax30 = axes('Position',[.05 .09 .31 .37],'Color','none');
ax31 = axes('Position',[.42 .09 .06 .37],'Color','none');

ax40 = axes('Position',[.54 .09 .31 .37],'Color','none');
ax41 = axes('Position',[.91 .09 .06 .37],'Color','none');




%--------------------------------------------------------------------------
%---   NAT/NAT | PLOT CONFUSION PERFORMANCE HISTOGRAM
%--------------------------------------------------------------------------

%---STACKED HISTOGRAM ---------------------
    axes(ax10)
ph0=area(PER(1).area(:,1),PER(1).area(:,2:5)); hold on;
    ax10.ColorOrderIndex = 1;
ph1=bar( PER(1).area(:,1), PER(1).area(:,2:5) ,'stacked');
    ylabel('Count'); 
    ax10.FontSize=16; 
    lh=legend([ph1],{'Case Miss','Ctrl Miss','Case Hit','Ctrl Hit'},...
    'Location','best');
    title([gnam ' | APOE34 | SNP REF/ALT'],'interpreter','none')
    hold on
    %pause(.1); 


%---MID LINES ON HISTOGRAM ---------------------
    Y = ax10.YLim(1):50:ax10.YLim(2);
    Xa = repmat(PER(1).caseMu,1,numel(Y));
    Xb = repmat(PER(1).ctrlMu,1,numel(Y));
ph2=plot(Xa,Y,'.-',Xb,Y,'.-','MarkerSize',30,'LineWidth',6);
    ph2(1).Color = [.5 .1 .1];
    ph2(2).Color = [.1 .2 .6];
    lh.String(5:6) = {'Case Mu','Ctrl Mu'};
    hold on
    %pause(.2); 


%---BAR GRAPH ---------------------
    axes(ax11)
ph3=bar([.6 2 3.4],[ PER(1).caseMu , PER(1).ctrlMu ,  PER(1).cacoZd],...
    .5,'FaceColor',[.31 .31 .31]); 
    %---BAR GRAPH FORMATTING ------
    ylabel('deviation from zero')
    ax11.YGrid = 'on';
    ax11.YLim = [-.5 .5]; 
    ax11.XTickLabels = {'CASE','CTRL','AbsMd'};
    ax11.XTickLabelRotation = 45;
    ax11.FontSize=16;
    % ax11.YTick = [0 25 50 75 100]; 
    %title('Performance Summary')





%--------------------------------------------------------------------------
%---   REF/REF | PLOT CONFUSION PERFORMANCE HISTOGRAM
%--------------------------------------------------------------------------


%---STACKED HISTOGRAM ---------------------
    axes(ax20)
ph1=area( PER(2).area(:,1), PER(2).area(:,2:5) ); hold on;
    ax20.ColorOrderIndex = 1;
ph1=bar( PER(2).area(:,1), PER(2).area(:,2:5) ,'stacked');
    ylabel('Count'); 
    ax20.FontSize=16; 
    lh=legend([ph1],{'Case Miss','Ctrl Miss','Case Hit','Ctrl Hit'},...
    'Location','best');
    title([gnam ' | APOE34 | SNP ALT/ALT'],'interpreter','none')
    hold on
    %pause(.1); 


%---MID LINES ON HISTOGRAM ---------------------
    Y = ax20.YLim(1):50:ax20.YLim(2);
    Xa = repmat(PER(2).caseMu,1,numel(Y));
    Xb = repmat(PER(2).ctrlMu,1,numel(Y));
ph2=plot(Xa,Y,'.-',Xb,Y,'.-','MarkerSize',30,'LineWidth',6);
    ph2(1).Color = [.5 .1 .1];
    ph2(2).Color = [.1 .2 .6];
    lh.String(5:6) = {'Case Mu','Ctrl Mu'};
    hold on
    %pause(.2); 



%---BAR GRAPH ---------------------
    axes(ax21)
ph3=bar([.6 2 3.4], [ PER(2).caseMu , PER(2).ctrlMu ,  PER(2).cacoZd],...
    .5,'FaceColor',[.31 .31 .31]); 
    %---BAR GRAPH FORMATTING ------
    grid on; ylabel('zero deviance weight')
    ax21.YLim = [-.5 .5]; 
    ax21.XTickLabels = {'CASE','CTRL','AbsMd'};
    ax21.XTickLabelRotation = 45;
    ax21.FontSize=16;
    % ax21.YTick = [0 25 50 75 100]; 
    %title('Performance Summary')







%--------------------------------------------------------------------------
%---   REF/ALT | PLOT CONFUSION PERFORMANCE HISTOGRAM
%--------------------------------------------------------------------------

%---STACKED HISTOGRAM ---------------------
    axes(ax30)
ph1=area(PER(3).area(:,1),PER(3).area(:,2:5)); hold on;
    ax30.ColorOrderIndex = 1;
ph1=bar( PER(3).area(:,1), PER(3).area(:,2:5) ,'stacked');
    ylabel('Count'); 
    ax30.FontSize=16; 
    lh=legend([ph1],{'Case Miss','Ctrl Miss','Case Hit','Ctrl Hit'},...
    'Location','best'); 
    title([gnam ' | APOE44 | SNP REF/ALT'],'interpreter','none')
    hold on
    %pause(.1); 
    

%---MID LINES ON HISTOGRAM ---------------------
    Y = ax30.YLim(1):50:ax30.YLim(2);
    Xa = repmat(PER(3).caseMu,1,numel(Y));
    Xb = repmat(PER(3).ctrlMu,1,numel(Y));
ph2=plot(Xa,Y,'.-',Xb,Y,'.-','MarkerSize',30,'LineWidth',6);
    ph2(1).Color = [.5 .1 .1];
    ph2(2).Color = [.1 .2 .6];
    lh.String(5:6) = {'Case Mu','Ctrl Mu'};
    hold on
    %pause(.2); 


%---BAR GRAPH ---------------------
    axes(ax31)
ph3=bar([.6 2 3.4], [ PER(3).caseMu , PER(3).ctrlMu ,  PER(3).cacoZd],...
    .5,'FaceColor',[.31 .31 .31]); 
    %---BAR GRAPH FORMATTING ------
    grid on; ylabel('zero deviance weight')
    ax31.YLim = [-.5 .5]; 
    ax31.XTickLabels = {'CASE','CTRL','AbsMd'};
    ax31.XTickLabelRotation = 45;
    ax31.FontSize=16;
    % ax31.YTick = [0 25 50 75 100]; 
    %title('Performance Summary')





%--------------------------------------------------------------------------
%---   ALT/ALT | PLOT CONFUSION PERFORMANCE HISTOGRAM
%--------------------------------------------------------------------------

%---STACKED HISTOGRAM ---------------------
    axes(ax40)
ph1=area(PER(4).area(:,1),PER(4).area(:,2:5)); hold on;
    ax40.ColorOrderIndex = 1;
ph1=bar( PER(4).area(:,1), PER(4).area(:,2:5) ,'stacked');
    ylabel('Count'); 
    ax40.FontSize=16; 
    lh=legend([ph1],{'Case Miss','Ctrl Miss','Case Hit','Ctrl Hit'},...
    'Location','best');
    title([gnam ' | APOE44 | SNP ALT/ALT'],'interpreter','none');
    hold on
    %pause(.1); 


%---MID LINES ON HISTOGRAM ---------------------
    Y = ax40.YLim(1):50:ax40.YLim(2);
    Xa = repmat(PER(4).caseMu,1,numel(Y));
    Xb = repmat(PER(4).ctrlMu,1,numel(Y));
ph2=plot(Xa,Y,'.-',Xb,Y,'.-','MarkerSize',30,'LineWidth',6);
    ph2(1).Color = [.5 .1 .1];
    ph2(2).Color = [.1 .2 .6];
    lh.String(5:6) = {'Case Mu','Ctrl Mu'};
    hold on
    %pause(.2); 


%---BAR GRAPH ---------------------
    axes(ax41)
ph3=bar([.6 2 3.4], [ PER(4).caseMu , PER(4).ctrlMu ,  PER(4).cacoZd],...
    .5,'FaceColor',[.31 .31 .31]); 
    %---BAR GRAPH FORMATTING ------
    grid on; ylabel('zero deviance weight')
    ax41.YLim = [-.5 .5]; 
    ax41.XTickLabels = {'CASE','CTRL','AbsMd'};
    ax41.XTickLabelRotation = 45;
    ax41.FontSize=16;
    % ax41.YTick = [0 25 50 75 100]; 
    %title('Performance Summary')
%----------------------------------------------------------------------
pause(1)
set(gcf, 'PaperPositionMode', 'auto');
finam = char(P.PERX.finames(ij));
saveas(gcf, [P.genox.DIRstatsimg P.f finam(1:end-4) '.png']);
pause(1)
%---------------------------------------------------------------------- 


end




%==========================================================================
%% xxx
%==========================================================================
