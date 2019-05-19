%% GENOS
%--------------------------------------------------------------------------
% SUMMARY TABLE OF THE 24 COHORTS
%--------------------------------------------------------------------------
%{
% 
% ID  CONSOR  STUDY      COHORT  CASE   CTRL    TOT  %CASE   EQL BRAK ID OK
% 01  DGC     Adult_Chng  ACT     323    945   1268    25    323   1  01  1
% 02  ADGC    AD_Centers  ADC    2438    817   3255    75    817   0  02  1
% 03  CHARGE  Athrosclro  ASC      39     18     57    68     18   0  03  0
% 04  CHARGE  Aus_Stroke  SKE     121      5    126    96      5   0  04  0
% 05  ADGC    ChiT_Aging  CHA      27    204    231    12     27   0  05  0
% 06  CHARGE  CardioHlth  CHS     250    583    833    30    250   1  06  1
% 07  ADGC    Hispanic_S  HSP     160    171    331    48    160   0  07  X
% 08  CHARGE  Erasmus_Er  ERF      45      0     45   100      0   0  08  0
% 09  CHARGE  Framingham  FHS     157    424    581    27    157   1  09  1
% 10  ADGC    Gene_Diffs  GDF     111     96    207    54     96   1  10  1
% 11  ADGC    NIA_LOAD    LOD     367    109    476    77    109   1  11  1
% 12  ADGC    Aging_Proj  MAP     138    277    415    33    138   1  12  1
% 13  ADGC    Mayo_Clini  MAY     250     99    349    72     99   1  13  1
% 14  ADGC    Miami_Univ  MIA     186     14    200    93     14   1  14  0
% 15  ADGC    AD_Genetic  MIR     316     15    331    95     15   0  15  0
% 16  ADGC    Mayo_cl_PD  MPD       0     20     20     0      0   1  16  0
% 17  ADGC    NationC_AD  NCA     160      0    160   100      0   1  17  0
% 18  ADGC    Wash_Unive  RAS      46      0     46   100      0   1  18  0
% 19  ADGC    Relig_Ordr  ROS     154    197    351    44    154   1  19  1
% 20  CHARGE  RotterdamS  RDS     276    813   1089    25    276   0  20  1
% 21  ADGC    Texas_AD_S  TAR     132     12    144    92     12   0  21  0
% 22  ADGC    Un_Toronto  TOR       9      0      9   100      0   0  22  0
% 23  ADGC    Vanderbilt  VAN     210     26    236    89     26   1  23  1
% 24  ADGC    WashNY_Age  WCA      34    116    150    23     34   0  24  1
% 
% 
% ADSPCOHORTS = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]
% GOODCOHORTS = [1 2       6 7   9 10 11 12 13                19 20       23 24]
% BRAKCOHORTS = [1         6     9 10 11 12 13 14    16 17 18 19          23   ]
% BAADCOHORTS = [    3 4 5     8               14 15 16 17 18       21 22      ]
% 
%}
%==========================================================================
%% STEP-1: LOAD THE DATASET
%==========================================================================
clc; close all; clear; rng('shuffle'); f = filesep;
P.root  = [f 'Users' f 'bradleymonk' f 'Documents' f 'MATLAB'];
P.home =  [P.root f 'GIT' f 'genomics' f 'genos'];
P.funs  = [P.home f 'genos_functions'];
P.data  = [P.home f 'genos_data'];
P.figs  = [P.home f 'genos_figures'];
P.mat1  = [P.data f 'APOE'];
addpath(join(string(struct2cell(P)),pathsep,1))
cd(P.home); clearvars -except P


% ADSP = load('GENOSDATAFINAL.mat');




% APOE 22 23 24 33 34 44
%-------------------------------
% 
% 50-VARIANT WINDOW
%---------------------
% P.datapath = ['/Users/bradleymonk/Documents/MATLAB/GIT/genomics/genos/genos_data/APOE'...
% '/APOE_22_23_24_33_34_44/APOE_22_23_24_33_34_44_WIN50_p01-25.mat'];
% P.alleles = 'APOE_22_23_24_33_34_44';
% P.APOE = [22 23 24 33 34 44];
% 
% 
% 20-VARIANT WINDOW
%---------------------
% P.datapath = ['/Users/bradleymonk/Documents/MATLAB/GIT/genomics/genos/genos_data/APOE'...
% '/APOE_22_23_24_33_34_44/APOE_22_23_24_33_34_44_WIN20_p01-10.mat'];
% P.alleles = 'APOE_22_23_24_33_34_44';
% P.APOE = [22 23 24 33 34 44];




% APOE 22 23 24 34 44
%-------------------------------
% 
% 50-VARIANT WINDOW
%---------------------
% P.datapath = ['/Users/bradleymonk/Documents/MATLAB/GIT/genomics/genos/genos_data/APOE'...
% '/APOE_22_23_24_34_44/APOE_22_23_24_34_44_WIN50_p16-25.mat'];
% P.alleles = 'APOE_22_23_24_34_44';
% P.APOE = [22 23 24 34 44];
% 
% 
% 20-VARIANT WINDOW
%---------------------
% P.datapath = ['/Users/bradleymonk/Documents/MATLAB/GIT/genomics/genos/genos_data/APOE'...
% '/APOE_22_23_24_34_44/APOE_22_23_24_34_44_WIN20_p01-10.mat'];
% P.alleles = 'APOE_22_23_24_34_44';
% P.APOE = [22 23 24 34 44];



% APOE 33
%-------------------------------
% 
% 50-VARIANT WINDOW
%---------------------
% P.datapath = ['/Users/bradleymonk/Documents/MATLAB/GIT/genomics/genos/genos_data/APOE'...
% '/APOE_33/APOE_33_WIN50_p01-10.mat'];
% P.alleles = 'APOE_33';
% P.APOE = [33];
% 
% 
% 20-VARIANT WINDOW
%---------------------
% P.datapath = ['/Users/bradleymonk/Documents/MATLAB/GIT/genomics/genos/genos_data/APOE'...
% '/APOE_33/APOE_33_WIN20_p01-25.mat'];
% P.alleles = 'APOE_33';
% P.APOE = [33];





ADSP = load([P.datapath]);


clc; disp('DATASETS LOADED');
who; clearvars -except P ADSP INFO


%==========================================================================
%%   CARBON COPY LOOP DATA FROM ADSP.STRUCT
%==========================================================================


LOOPDATA = ADSP.LOOPDATA;

INFO = ADSP.INFO;


clc; clearvars -except P ADSP INFO LOOPDATA




%==========================================================================
%%   CHECK FOR USAGE OF BAD COHORTS OR BAD APOE ALLELES
%==========================================================================


TRCASE = INFO.PHETRCASE{1,1};
TRCTRL = INFO.PHETRCTRL{1,1};
TECASE = INFO.PHETECASE{1,1};
TECTRL = INFO.PHETECTRL{1,1};

TRCASE.TRTECACO = zeros(size(TRCASE,1),1) + 1;
TRCTRL.TRTECACO = zeros(size(TRCTRL,1),1) + 2;
TECASE.TRTECACO = zeros(size(TECASE,1),1) + 3;
TECTRL.TRTECACO = zeros(size(TECTRL,1),1) + 4;


PHE = [TRCASE; TRCTRL; TECASE; TECTRL];

COH = PHE.COHORTNUM;

APO = PHE.APOE;


ADSPCOHORTS = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24];
GOODCOHORTS = [1 2       6 7   9 10 11 12 13                19 20       23 24];
BRAKCOHORTS = [1         6     9 10 11 12 13 14    16 17 18 19          23   ];
BAADCOHORTS = [    3 4 5     8               14 15 16 17 18       21 22      ];



clc;
fprintf('\n APOE alleles found in dataset: %s \n', num2str(unique(APO')) );



fprintf('\n N bad APOE alleles in dataset: %s \n',...
    num2str(any(~sum(APO == P.APOE,2)))  );



fprintf('\n COHORTS found in dataset: %s \n', num2str(unique(COH')) );



fprintf('\n N bad COHORTS in dataset: %s \n',...
    num2str(any(~sum(COH == GOODCOHORTS,2)))  );







%==========================================================================
%% GET MEAN & STDEV & SEM STATS FOR LOOP
%==========================================================================
clc; clearvars -except P ADSP INFO LOOPDATA


INFO.NMAXV = size(LOOPDATA.CATRMEAN,1);
INFO.NLOOPS  = size(LOOPDATA.CATRMEAN,2);

INFO.NLOOPSi = 1:INFO.NLOOPS;
INFO.VARi    = 1:INFO.NMAXV;



STATS.muCATRMEAN = nanmean(LOOPDATA.CATRMEAN,2);
STATS.muCATRHIMU = nanmean(LOOPDATA.CATRHIMU,2);
STATS.muCATRHIPO = nanmean(LOOPDATA.CATRHIPO,2);
STATS.muCAHOMEAN = nanmean(LOOPDATA.CAHOMEAN,2);
STATS.muCAHOHIMU = nanmean(LOOPDATA.CAHOHIMU,2);
STATS.muCAHOHIPO = nanmean(LOOPDATA.CAHOHIPO,2);

STATS.sdCATRMEAN = nanstd(LOOPDATA.CATRMEAN,[],2);
STATS.sdCATRHIMU = nanstd(LOOPDATA.CATRHIMU,[],2);
STATS.sdCATRHIPO = nanstd(LOOPDATA.CATRHIPO,[],2);
STATS.sdCAHOMEAN = nanstd(LOOPDATA.CAHOMEAN,[],2);
STATS.sdCAHOHIMU = nanstd(LOOPDATA.CAHOHIMU,[],2);
STATS.sdCAHOHIPO = nanstd(LOOPDATA.CAHOHIPO,[],2);

STATS.seCATRMEAN = nanstd(LOOPDATA.CATRMEAN,[],2) ./ sqrt(INFO.NLOOPS);
STATS.seCATRHIMU = nanstd(LOOPDATA.CATRHIMU,[],2) ./ sqrt(INFO.NLOOPS);
STATS.seCATRHIPO = nanstd(LOOPDATA.CATRHIPO,[],2) ./ sqrt(INFO.NLOOPS);
STATS.seCAHOMEAN = nanstd(LOOPDATA.CAHOMEAN,[],2) ./ sqrt(INFO.NLOOPS);
STATS.seCAHOHIMU = nanstd(LOOPDATA.CAHOHIMU,[],2) ./ sqrt(INFO.NLOOPS);
STATS.seCAHOHIPO = nanstd(LOOPDATA.CAHOHIPO,[],2) ./ sqrt(INFO.NLOOPS);



STATS.muCOTRMEAN = nanmean(LOOPDATA.COTRMEAN,2);
STATS.muCOTRHIMU = nanmean(LOOPDATA.COTRHIMU,2);
STATS.muCOTRHIPO = nanmean(LOOPDATA.COTRHIPO,2);
STATS.muCOHOMEAN = nanmean(LOOPDATA.COHOMEAN,2);
STATS.muCOHOHIMU = nanmean(LOOPDATA.COHOHIMU,2);
STATS.muCOHOHIPO = nanmean(LOOPDATA.COHOHIPO,2);

STATS.sdCOTRMEAN = nanstd(LOOPDATA.COTRMEAN,[],2);
STATS.sdCOTRHIMU = nanstd(LOOPDATA.COTRHIMU,[],2);
STATS.sdCOTRHIPO = nanstd(LOOPDATA.COTRHIPO,[],2);
STATS.sdCOHOMEAN = nanstd(LOOPDATA.COHOMEAN,[],2);
STATS.sdCOHOHIMU = nanstd(LOOPDATA.COHOHIMU,[],2);
STATS.sdCOHOHIPO = nanstd(LOOPDATA.COHOHIPO,[],2);

STATS.seCOTRMEAN = nanstd(LOOPDATA.COTRMEAN,[],2) ./ sqrt(INFO.NLOOPS);
STATS.seCOTRHIMU = nanstd(LOOPDATA.COTRHIMU,[],2) ./ sqrt(INFO.NLOOPS);
STATS.seCOTRHIPO = nanstd(LOOPDATA.COTRHIPO,[],2) ./ sqrt(INFO.NLOOPS);
STATS.seCOHOMEAN = nanstd(LOOPDATA.COHOMEAN,[],2) ./ sqrt(INFO.NLOOPS);
STATS.seCOHOHIMU = nanstd(LOOPDATA.COHOHIMU,[],2) ./ sqrt(INFO.NLOOPS);
STATS.seCOHOHIPO = nanstd(LOOPDATA.COHOHIPO,[],2) ./ sqrt(INFO.NLOOPS);






clearvars -except P ADSP INFO LOOPDATA STATS
%==========================================================================
%%              SUMMARY DATA FOR PRISM
%==========================================================================
clearvars -except P ADSP INFO LOOPDATA STATS




AA.CaHoMu   = STATS.muCAHOMEAN(1:151,:);
AA.CaHoSe   = STATS.seCAHOMEAN(1:151,:);

AA.CaHoHiMu = STATS.muCAHOHIMU(1:151,:);
AA.CaHoHiSe = STATS.seCAHOHIMU(1:151,:);


AA.CoHoMu   = STATS.muCOHOMEAN(1:151,:);
AA.CoHoSe   = STATS.seCOHOMEAN(1:151,:);

AA.CoHoHiMu = STATS.muCOHOHIMU(1:151,:);
AA.CoHoHiSe = STATS.seCOHOHIMU(1:151,:);


AA.ALL = [AA.CaHoMu, AA.CaHoSe, AA.CoHoMu, AA.CoHoSe,...
          AA.CaHoHiMu, AA.CaHoHiSe, AA.CoHoHiMu, AA.CoHoHiSe];


TR_ROC = mean(LOOPDATA.TR_ALL_STATS,3);


return
%==========================================================================
%==========================================================================
%==========================================================================
%%                          HILO  4-PACK
%==========================================================================
clearvars -except P ADSP INFO LOOPDATA STATS



ESz = 2; % Error bar cap size
szO1 = 350;
szO2 = 50;
PH = flipud(jet); 
SaveFig = 1==1;
XL = [-5 INFO.NMAXV];
YL = [0 1];


%################   FOUR PACK   ################
close all
fh01 = figure('Units','normalized','OuterPosition',[.01 .05 .95 .88],'Color','w');
ax01 = axes('Position',[.06 .56 .4 .4],'Color','none');
ax02 = axes('Position',[.56 .56 .4 .4],'Color','none');
ax03 = axes('Position',[.06 .06 .4 .4],'Color','none');
ax04 = axes('Position',[.56 .06 .4 .4],'Color','none');

    axes(ax01);
    line([-10 500],[.5 .5],'LineStyle',':','LineWidth',1.5); hold on;
errorbar( INFO.VARi, STATS.muCATRMEAN, STATS.seCATRMEAN,...
    '.','MarkerSize',4,'LineStyle','none','CapSize',4,'Color',[.05 .2 .6]); hold on;
scatter(INFO.VARi, STATS.muCATRMEAN,szO1,[.05 .2 .6],'.');
    hold on;
errorbar( INFO.VARi, STATS.muCOTRMEAN, STATS.seCOTRMEAN,...
    '.','MarkerSize',4,'LineStyle','none','CapSize',4,'Color',[.6 .2 .05]); hold on;
scatter(INFO.VARi,   STATS.muCOTRMEAN,szO1,[.6 .2 .1],'.');

    title('TRAINING SUBSET (ALL CONFIDENCE)')
    xlabel('Starting index of 30-SNP sliding window')
    ylabel('Proportion correct')
    ax01.XLim=XL; ax01.YLim=YL; hold on;



    axes(ax02);
    line([-10 500],[.5 .5],'LineStyle',':','LineWidth',1.5); hold on;
errorbar(INFO.VARi, STATS.muCAHOMEAN, STATS.seCAHOMEAN,...
    '.','MarkerSize',4,'LineStyle','none','CapSize',4,'Color',[.05 .2 .6]); hold on;
scatter(INFO.VARi,STATS.muCAHOMEAN,szO1,[.05 .2 .6],'.');
    hold on;
errorbar(INFO.VARi, STATS.muCOHOMEAN, STATS.seCOHOMEAN,...
    '.','MarkerSize',4,'LineStyle','none','CapSize',4,'Color',[.6 .2 .05]); hold on;
scatter(INFO.VARi,STATS.muCOHOMEAN,szO1,[.6 .2 .05],'.');

    title('HOLDOUT SUBSET (ALL CONFIDENCE)')
    xlabel('Starting index of 30-SNP sliding window')
    ylabel('Proportion correct')
    ax02.XLim=XL; ax02.YLim=YL; hold on;





    axes(ax03);
line([-10 500],[.5 .5],'LineStyle',':','LineWidth',1.5); hold on;
errorbar( INFO.VARi, STATS.muCATRHIMU, STATS.seCATRHIMU,...
    '.','MarkerSize',4,'LineStyle','none','CapSize',4,'Color',[.05 .2 .6]); hold on;
scatter(INFO.VARi, STATS.muCATRHIMU,szO1,[.05 .2 .6],'.');
    hold on;
errorbar( INFO.VARi, STATS.muCOTRHIMU, STATS.seCOTRHIMU,...
    '.','MarkerSize',4,'LineStyle','none','CapSize',4,'Color',[.6 .2 .05]); hold on;
scatter(INFO.VARi,   STATS.muCOTRHIMU,szO1,[.6 .2 .1],'.');

    title('TRAINING SUBSET (HIGH CONFIDENCE)')
    xlabel('Starting index of 30-SNP sliding window')
    ylabel('Proportion correct')
    ax03.XLim=XL; ax03.YLim=YL; hold on;



    axes(ax04); 
line([-10 500],[.5 .5],'LineStyle',':','LineWidth',1.5); hold on;
errorbar(INFO.VARi, STATS.muCAHOHIMU, STATS.seCAHOHIMU,...
    '.','MarkerSize',4,'LineStyle','none','CapSize',4,'Color',[.05 .2 .6]); hold on;
     hold on
scatter(INFO.VARi, STATS.muCAHOHIMU,szO1,[.05 .2 .6],'.');
    hold on;
errorbar(INFO.VARi, STATS.muCOHOHIMU, STATS.seCOHOHIMU,...
    '.','MarkerSize',4,'LineStyle','none','CapSize',4,'Color',[.6 .2 .05]); hold on;
    hold on
scatter(INFO.VARi,STATS.muCOHOHIMU,szO1,[.6 .2 .05],'.');

    
    title('HOLDOUT SUBSET (HIGH CONFIDENCE)')
    xlabel('Starting index of 30-SNP sliding window')
    ylabel('Proportion correct')
    ax04.XLim=XL; ax04.YLim=YL; hold on;










%------------------------------------------%
pause(1)
if SaveFig
set(gcf, 'PaperPositionMode', 'auto');
dt=char(datetime(datetime,'Format','yyyy-MM-dd-HH-mm-ss'));
saveas(gcf, ['/Users/bradleymonk/Desktop/' P.alleles dt '.png']);
pause(1)
end
%------------------------------------------%









%==========================================================================
%% 2-PACK GRAPHS TRAINING & HOLDOUT PROPORTION CORRECT x N-VARIANT LOCI
%  STRICT AXES LIMITS
%==========================================================================
% clearvars -except P ADSP INFO PHEN LOCI CASE CTRL USNP LOOPDATA STATS VARi


clc; close all
%################   TWO PACK   ################
fh01 = figure('Units','pixels','Position',[10 35 1290 400],'Color','w');
ax01 = axes('Position',[.06 .17 .38 .78],'Color','none','FontSize',20); 
ax01.Colormap = parula; hold on;

ax02 = axes('Position',[.55 .17 .38 .78],'Color','none','FontSize',20); 
ax02.Colormap = parula; hold on;





%----- AXIS-1 50:50 CHANCE LINE
axes(ax01);
line([-10 -10 -10 INFO.NMAXV],[0 1 .5 .5],'LineStyle',':','LineWidth',1.5);
hold on



%----- TRAINING CASE ALL
    axes(ax01);
ph01 = errorbar(INFO.VARi,STATS.muCATRMEAN, STATS.seCATRMEAN,...
    '.','MarkerSize',4,'LineStyle','none','CapSize',ESz,'Color',[.5 .5 .7]); hold on;
    % 'MarkerEdgeColor','red','MarkerFaceColor','red'
scatter(INFO.VARi,STATS.muCATRMEAN,szO1,[.01 .1 .5],'.'); hold on; 

%----- TRAINING CASE HICONF
ph02 = errorbar(INFO.VARi, STATS.muCATRHIMU, STATS.seCATRHIMU,...
    '.','MarkerSize',4,'LineStyle','none','CapSize',ESz,'Color',[.5 .5 .7]); hold on;
    hold on;
scatter(INFO.VARi,STATS.muCATRHIMU,szO1,STATS.muCATRHIPO,'.'); hold on;
    set(ax01,'colormap',parula); hold on;
    ax01.XLim=XL; ax01.YLim=YL; hold on;
     





%----- TRAINING CONTROL ALL
    axes(ax01);
errorbar(INFO.VARi,STATS.muCOTRMEAN, STATS.seCOTRMEAN,...
    'o','MarkerSize',4,'LineStyle','none','CapSize',ESz,'Color',[.5 .5 .7]); hold on;
    ax01.XLim=XL; hold on; pause(.1);

scatter(INFO.VARi,STATS.muCOTRMEAN,szO2,[.01 .1 .5],'o');
    ax01.XLim=XL; hold on; pause(.1);

%----- TRAINING CONTROL HICONF
errorbar(INFO.VARi, STATS.muCOTRHIMU, STATS.seCOTRHIMU,...
    'o','MarkerSize',4,'LineStyle','none','CapSize',ESz,'Color',[.5 .5 .7]); hold on;
   ax01.XLim=XL; hold on; pause(.1);
    
scatter(INFO.VARi,STATS.muCOTRHIMU,szO2,STATS.muCOTRHIPO,'o');
    ax01.XLim=XL; hold on; pause(.1);

    set(ax01,'colormap',parula); hold on;
    ax01.XLim=XL; ax01.YLim=YL; hold on;

    %title('TRAINING SUBSET (HIGH & ALL CONFIDENCE)')
    xlabel('Sliding-window start index')
    ylabel('Proportion correct')
    axes(ax01); cb01 = colorbar;
    %cb01.Label.String = 'Proportion of Sample';
    %cb01.Label.FontSize = 12;
    cb01.Label.HorizontalAlignment = 'center';
    cb01.Label.VerticalAlignment = 'bottom';
    cb01.Label.Rotation = -90;

    set(ax01,'colormap',parula); hold on;
    ax01.XLim=XL; ax01.YLim=YL; hold on;








%----- AXIS-2 50:50 CHANCE LINE
axes(ax02);
line([-10 -10 -10 INFO.NMAXV],[0 1 .5 .5],'LineStyle',':','LineWidth',1.5);
hold on




%----- HOLDOUT CASE ALL
    axes(ax02);
ph01 = errorbar(INFO.VARi,STATS.muCAHOMEAN, STATS.seCAHOMEAN,...
    '.','MarkerSize',4,'LineStyle','none','CapSize',ESz,'Color',[.5 .5 .7]); hold on;
    hold on
scatter(INFO.VARi,STATS.muCAHOMEAN,szO1,[.01 .1 .5],'.');
    hold on; 
%----- HOLDOUT CASE HICONF
line([-10 -10 -10 INFO.NMAXV],[0 1 .5 .5],'LineStyle',':','LineWidth',1.5);
    hold on
ph02 = errorbar(INFO.VARi, STATS.muCAHOHIMU, STATS.seCAHOHIMU,...
    '.','MarkerSize',4,'LineStyle','none','CapSize',ESz,'Color',[.5 .5 .7]); hold on;
    hold on;
scatter(INFO.VARi,STATS.muCAHOHIMU,szO1,STATS.muCAHOHIPO,'.');
    set(ax02,'colormap',parula); hold on;
    ax02.XLim=XL; ax02.YLim=YL; hold on;





%----- HOLDOUT CONTROL ALL
    axes(ax02);
errorbar(INFO.VARi,STATS.muCOHOMEAN, STATS.seCOHOMEAN,...
    'o','MarkerSize',4,'LineStyle','none','CapSize',ESz,'Color',[.5 .5 .7]); hold on;
    ax02.XLim= [-10 INFO.NMAXV]; hold on; pause(.1);

scatter(INFO.VARi,STATS.muCOHOMEAN,szO2,[.01 .1 .5],'o');
    ax02.XLim= [-10 INFO.NMAXV]; hold on; pause(.1);

%----- HOLDOUT CONTROL HICONF
errorbar(INFO.VARi, STATS.muCOHOHIMU, STATS.seCOHOHIMU,...
    'o','MarkerSize',4,'LineStyle','none','CapSize',ESz,'Color',[.5 .5 .7]); hold on;
   ax02.XLim= [-10 INFO.NMAXV]; hold on; pause(.1);
    
scatter(INFO.VARi,STATS.muCOHOHIMU,szO2,STATS.muCOHOHIPO,'o');
    ax02.XLim= [-10 INFO.NMAXV]; hold on; pause(.1); hold on

    set(ax02,'colormap',parula); hold on;
    ax02.XLim=XL; ax02.YLim=YL; hold on;



% title('HOLDOUT SUBSET (HIGH & ALL CONFIDENCE)')
xlabel('Sliding-window start index')
ylabel('Proportion correct')
axes(ax02); cb01 = colorbar;
cb01.Label.String = 'High-confidence output pct.';
cb01.Label.FontSize = 20;
% cb01.FontSize = 16;
cb01.Label.HorizontalAlignment = 'center';
cb01.Label.VerticalAlignment = 'bottom';
cb01.Limits = [.15 .30];
cb01.TickLabels = {'15','20','25','30'};
cb01.Label.Rotation = -90;
set(ax02,'colormap',parula); hold on;
ax02.XLim=XL; ax02.YLim=YL; hold on;




% ax01.CLim = ax01.CLim./1.5;
% ax02.CLim = ax02.CLim./1.5;


%------------------------------------------%
pause(1)
if SaveFig
set(gcf, 'PaperPositionMode', 'auto');
dt=char(datetime(datetime,'Format','yyyy-MM-dd-HH-mm-ss'));
saveas(gcf, ['/Users/bradleymonk/Desktop/' P.alleles dt '_2PACK.png']);
pause(1)
end
%------------------------------------------%
