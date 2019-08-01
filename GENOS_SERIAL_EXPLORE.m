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
% P.data = 'F:\GENOSDATA';
addpath(join(string(struct2cell(P)),pathsep,1))
cd(P.home); P.f = filesep;



tic
ADSP = load('GENOSDATA.mat');
toc



clearvars -except P ADSP
%==========================================================================
%% IMPORT iSNP
%==========================================================================
clc; clearvars -except P ADSP
%iLOC = readtable([P.serial.DIRsnp P.f 'iLOC.csv']);
%iLOC.CHRPOS = uint64(iLOC.CHRPOS);
%iSNP = SNP.iSNP;
%LOOP = SNP.LOOP;
%save([P.serial.DIRsnp P.f 'SNP.mat' ],'iSNP','iLOC','LOOP','-v6');



RUN = 'RUN5';



P.serial.DIRroot = [P.home P.f 'genos_data' P.f 'SERIAL' P.f RUN];
P.serial.DIRmat  = [P.serial.DIRroot P.f 'MAT' ];
P.serial.DIRimg  = [P.serial.DIRroot P.f 'IMG' ];
P.serial.DIRout  = [P.serial.DIRroot P.f 'OUT' ];
P.serial.DIRsnp  = [P.serial.DIRroot P.f 'iSNP'];
SNP = load([P.serial.DIRsnp P.f 'SNP.mat' ],'iSNP','iLOC','LOOP');




clearvars -except P ADSP SNP
%==========================================================================
%% IMPORT OLD DATA
%==========================================================================
%{
clc; clearvars -except P ADSP SNP




OLD.SNP = readtable([P.serial.DIRsnp P.f 'SERIALSNPS_ALTmREFv1.csv']);
OLD.DIF = table2array(OLD.SNP(:,23:1022));
OLD.SNP(:,23:1022) = [];
OLD.SNP.DIF = OLD.DIF;
r = sum((OLD.SNP.DIF==0),2);
OLD.SNP(r>0,:) = [];

SNP.oSNP = OLD.SNP;



clearvars -except P ADSP SNP
%}







%==========================================================================
%% PLOT A HISTOGRAM OF ALL VALUES
%==========================================================================
close all; clc; clearvars -except P ADSP SNP
%----------------
iSNP = SNP.iSNP;
iLOC = SNP.iLOC;
%----------------
PVX = iSNP.PVX;
REF = iSNP.REF;
ALT = iSNP.ALT;
DIF = iSNP.DIF;
NAT = iSNP.NAT;
%----------------





close all
fh01 = figure('Units','normalized','OuterPosition',[.01 .05 .98 .90],...
              'Color','w','MenuBar','none');
ax01 = axes('Position',[.03 .56 .44 .4],'Color','none');
ax02 = axes('Position',[.53 .56 .44 .4],'Color','none');
ax03 = axes('Position',[.03 .06 .44 .4],'Color','none');
ax04 = axes('Position',[.53 .06 .44 .4],'Color','none');


axes(ax01)
h=histogram(iSNP.REF,'LineWidth',.1,'Normalization','pdf','LineWidth',.01);
axis([-.5 .5 0 4])

axes(ax02)
h=histogram(iSNP.ALT,'LineWidth',.1,'Normalization','pdf','LineWidth',.01);
axis([-.5 .5 0 4])

axes(ax03)
h=histogram(iSNP.DIF,'LineWidth',.1,'Normalization','pdf','LineWidth',.01);
axis([-.5 .5 0 7])

axes(ax04)
h=histogram(iSNP.NAT,'LineWidth',.1,'Normalization','pdf','LineWidth',.01);
axis([-.5 .5 0 4])




%==========================================================================
%% PLOT A HISTOGRAM OF ALL VALUES VS. MEAN SNP VALUES
%==========================================================================
close all; clc; clearvars -except P ADSP SNP
%----------------
iSNP = SNP.iSNP;
iLOC = SNP.iLOC;
%----------------




close all
fh01 = figure('Units','normalized','OuterPosition',[.01 .04 .7 .94],'Color','w');
ax01 = axes('Position',[.07 .56 .40 .40],'Color','none');
ax02 = axes('Position',[.56 .56 .40 .40],'Color','none');
ax03 = axes('Position',[.07 .06 .40 .40],'Color','none');
ax04 = axes('Position',[.56 .06 .40 .40],'Color','none');
edges = (-.5:.02:.5);

axes(ax01)
h=histogram(iSNP.REF,(-.5:.01:.5),'LineWidth',.1,'Normalization','probability');
ylabel('P(x)'); axis([-.5 .5 0 .05])

axes(ax02)
h=histogram(iSNP.ALT,(-.5:.01:.5),'LineWidth',.1,'Normalization','probability');
ylabel('P(x)'); axis([-.5 .5 0 .05])

axes(ax03)
h=histogram(mean(iSNP.REF),(-.5:.005:.5),'LineWidth',.1,'Normalization','probability');
ylabel('P(x)'); axis([-.5 .5 0 .5])

axes(ax04)
h=histogram(mean(iSNP.ALT),(-.5:.01:.5),'LineWidth',.1,'Normalization','probability');
ylabel('P(x)'); axis([-.5 .5 0 .5])



%################   TWO PACK   ################
fh1 = figure('Units','normalized','OuterPosition',[.01 .06 .9 .7],'Color','w');
ax1 = axes('Position',[.05 .09 .42 .83],'Color','none');
ax2 = axes('Position',[.55 .09 .42 .83],'Color','none');


axes(ax1);  axis([-.5 .5 -.5 .5]); hold on;
h=binscatter((iSNP.REF(:)),(iSNP.ALT(:)),'XLimits',[-.5 5],'YLimits',[-.5 5],'NumBins',[150 150]);
xlabel('REF'); ylabel('ALT');
axis([-.5 .5 -.5 .5])
h.ShowEmptyBins = 'on';
cmocean('thermo')
colorbar


axes(ax2); axis([-.5 .5 -.5 .5]); hold on;
h=binscatter(mean(iSNP.REF),mean(iSNP.ALT),'XLimits',[-.5 5],'YLimits',[-.5 5],'NumBins',[100 100]);
xlabel('mean(REF)'); ylabel('mean(ALT)');
axis([-.5 .5 -.5 .5]);
h.ShowEmptyBins = 'on';
cmocean('thermo')
colorbar
% ax2.Colormap = ax1.Colormap; colorbar
% h.ShowEmptyBins = 'on';
% ax2.ColorScale = 'log';







%==========================================================================
%% ARE DIF SCORES DIFFERENT BETWEEN CASES AND CONTROLS
%==========================================================================
close all; clc; clearvars -except P ADSP SNP
%----------------
iSNP = SNP.iSNP;
iLOC = SNP.iLOC;
%----------------



close all
fh01 = figure('Units','normalized','OuterPosition',[.01 .05 .6 .90],...
              'Color','w','MenuBar','none');
ax01 = axes('Position',[.06 .56 .4 .4],'Color','none');
ax02 = axes('Position',[.56 .56 .4 .4],'Color','none');
ax03 = axes('Position',[.06 .06 .4 .4],'Color','none');
ax04 = axes('Position',[.56 .06 .4 .4],'Color','none');


axes(ax01)
scatter(  mean(iSNP.REF(iSNP.AD==0,:)), mean(iSNP.REF(iSNP.AD==1,:)),  100,'.')
axis([-.4 .4 -.4 .4])
hold on; line([-.4 .4], [-.4 .4])
xlabel('\fontsize{16}  MEAN REF CTRL')
ylabel('\fontsize{16}  MEAN REF CASE')

axes(ax02)
scatter( mean(iSNP.ALT(iSNP.AD==0,:)), mean(iSNP.ALT(iSNP.AD==1,:)),  100,'.')
axis([-.4 .4 -.4 .4])
hold on; line([-.4 .4], [-.4 .4])
xlabel('\fontsize{16}   MEAN ALT CTRL')
ylabel('\fontsize{16}   MEAN ALT CASE')

axes(ax03)
scatter( mean(iSNP.DIF(iSNP.AD==0,:)), mean(iSNP.DIF(iSNP.AD==1,:)),  100,'.')
axis([-.4 .4 -.4 .4])
hold on; line([-.4 .4], [-.4 .4])
xlabel('\fontsize{16}   MEAN DIF CTRL')
ylabel('\fontsize{16}   MEAN DIF CASE')

axes(ax04)
scatter( mean(iSNP.NAT(iSNP.AD==0,:)), mean(iSNP.NAT(iSNP.AD==1,:)),  100,'.')
axis([-.4 .4 -.4 .4])
hold on; line([-.4 .4], [-.4 .4])
xlabel('\fontsize{16} MEAN NAT CTRL')
ylabel('\fontsize{16} MEAN NAT CASE')





%==========================================================================
%% ARE REF|ALT|DIF SCORES CORRELATED WITH APOE ALLELE STATUS
%==========================================================================
clc; clearvars -except P ADSP SNP
%----------------
iSNP = SNP.iSNP;
iLOC = SNP.iLOC;
%----------------


CTRLAPOE22 = mean(iSNP.DIF( (iSNP.AD==0) & (iSNP.APOE==22),:));
CASEAPOE22 = mean(iSNP.DIF( (iSNP.AD==1) & (iSNP.APOE==22),:));
CTRLAPOE23 = mean(iSNP.DIF( (iSNP.AD==0) & (iSNP.APOE==23),:));
CASEAPOE23 = mean(iSNP.DIF( (iSNP.AD==1) & (iSNP.APOE==23),:));
CTRLAPOE24 = mean(iSNP.DIF( (iSNP.AD==0) & (iSNP.APOE==24),:));
CASEAPOE24 = mean(iSNP.DIF( (iSNP.AD==1) & (iSNP.APOE==24),:));
CTRLAPOE33 = mean(iSNP.DIF( (iSNP.AD==0) & (iSNP.APOE==33),:));
CASEAPOE33 = mean(iSNP.DIF( (iSNP.AD==1) & (iSNP.APOE==33),:));
CTRLAPOE34 = mean(iSNP.DIF( (iSNP.AD==0) & (iSNP.APOE==34),:));
CASEAPOE34 = mean(iSNP.DIF( (iSNP.AD==1) & (iSNP.APOE==34),:));
CTRLAPOE44 = mean(iSNP.DIF( (iSNP.AD==0) & (iSNP.APOE==44),:));
CASEAPOE44 = mean(iSNP.DIF( (iSNP.AD==1) & (iSNP.APOE==44),:));



%################   SIX PACK   ################
close all
fh01 = figure('Units','normalized','OuterPosition',[.01 .04 .95 .95],'Color','w');
ax01 = axes('Position',[.04 .56 .25 .35],'Color','none');
ax02 = axes('Position',[.35 .56 .25 .35],'Color','none');
ax03 = axes('Position',[.70 .56 .25 .35],'Color','none');
ax04 = axes('Position',[.04 .06 .25 .35],'Color','none');
ax05 = axes('Position',[.35 .06 .25 .35],'Color','none');
ax06 = axes('Position',[.70 .06 .25 .35],'Color','none');


    axes(ax01)
scatter(CTRLAPOE22,CASEAPOE22,[],(CTRLAPOE22-CASEAPOE22),'filled')
    axis([-.4 .4 -.4 .4])
    hold on; line([-.4 .4], [-.4 .4])
    xlabel('\fontsize{16} MEAN DIF CTRL')
    ylabel('\fontsize{16} MEAN DIF CASE')
    title('APOE 22');

    axes(ax02)
scatter(CTRLAPOE23,CASEAPOE23,[],(CTRLAPOE23-CASEAPOE23),'.')
    axis([-.4 .4 -.4 .4])
    hold on; line([-.4 .4], [-.4 .4])
    xlabel('\fontsize{16} MEAN DIF CTRL')
    ylabel('\fontsize{16} MEAN DIF CASE')
    title('APOE 23');

    axes(ax03)
scatter(CTRLAPOE24,CASEAPOE24,[],(CTRLAPOE24-CASEAPOE24),'.')
    axis([-.4 .4 -.4 .4])
    hold on; line([-.4 .4], [-.4 .4])
    xlabel('\fontsize{16} MEAN DIF CTRL')
    ylabel('\fontsize{16} MEAN DIF CASE')
    title('APOE 24');


    axes(ax04)
scatter(CTRLAPOE33,CASEAPOE33,[],(CTRLAPOE33-CASEAPOE33),'.')
    axis([-.4 .4 -.4 .4])
    hold on; line([-.4 .4], [-.4 .4])
    xlabel('\fontsize{16} MEAN DIF CTRL')
    ylabel('\fontsize{16} MEAN DIF CASE')
    title('APOE 33');

    axes(ax05)
scatter(CTRLAPOE34,CASEAPOE34,[],(CTRLAPOE33-CASEAPOE33),'.')
    axis([-.4 .4 -.4 .4])
    hold on; line([-.4 .4], [-.4 .4])
    xlabel('\fontsize{16} MEAN DIF CTRL')
    ylabel('\fontsize{16} MEAN DIF CASE')
    title('APOE 34');


    axes(ax06)
scatter(CTRLAPOE44,CASEAPOE44,[],(CTRLAPOE44-CASEAPOE44),'.')
    axis([-.4 .4 -.4 .4])
    hold on; line([-.4 .4], [-.4 .4])
    xlabel('\fontsize{16} MEAN DIF CTRL')
    ylabel('\fontsize{16} MEAN DIF CASE')
    title('APOE 44');






%==========================================================================
%% GET THE AVERAGE REF & ALT & DIF & NAT FOR EACH LOCI
%==========================================================================
clc; clearvars -except P ADSP SNP
%----------------
iSNP = SNP.iSNP;
iLOC = SNP.iLOC;
%----------------




iLOC.muREF = mean(iSNP.REF)';
iLOC.muALT = mean(iSNP.ALT)';
iLOC.muDIF = mean(iSNP.DIF)';
% iLOC.muNAT = mean(iSNP.NAT)';



close all 
fh1=figure('Units','normalized','OuterPosition',[.01 .05 .95 .6],'Color','w');
t = uitable(fh1,'Data',iLOC.Properties.VariableNames,'Units','normalized','Position',[.01 .01 .98 .95]);



%==========================================================================
%% WHAT ARE THE TOP 100 DIF SHIFT DISTRIBUTIONS
%==========================================================================
%{
clc; clearvars -except P ADSP SNP
%----------------
iSNP = SNP.iSNP;
iLOC = SNP.iLOC;
%----------------

iSNP.DIF = iSNP.ALT - iSNP.REF;


muDIF = mean(iSNP.DIF)';

[B,I] = sort(muDIF,'descend');


muDIF = muDIF(I);

close all
histogram(iSNP.DIF(:,I(1)))

topLOC = iLOC(I,:);

topLOC.muDIF = muDIF;


close all 
fh1=figure('Units','normalized','OuterPosition',[.01 .05 .95 .6],'Color','w');
t = uitable(fh1,'Data',topLOC.Properties.VariableNames,'Units','normalized','Position',[.01 .01 .98 .95]);
%}






%==========================================================================
%% IS BRAAK SCORE CORRELATED WITH REF|ALT|DIF SCORES
%==========================================================================
clc; clearvars -except P ADSP SNP
%----------------
iSNP = SNP.iSNP;
iLOC = SNP.iLOC;
%----------------
SNPS = iSNP.ALT;
nR  = 20;
nL  = 20;
%----------------


% histogram(iSNP.ALT)


muSNPS = mean(SNPS);
[B,LI] = sort(muSNPS);
RI = fliplr(LI);

Ri  = RI(1:nR); 
Li  = LI(1:nL);       


RiSHIFT   = SNPS(:,Ri);
LiSHIFT   = SNPS(:,Li);




% GET AVERAGE ALT SCORE PER BRAAK GROUP
HASBRAAK = ~isnan(iSNP.BRAAK);
BRAAK = iSNP.BRAAK(HASBRAAK);
BRAAK(BRAAK==0) = 1;
BRAAK(BRAAK==2) = 1;


AD       = iSNP.AD(HASBRAAK);
RSHIFT   = RiSHIFT(HASBRAAK,:);
LSHIFT   = LiSHIFT(HASBRAAK,:);
RLSHIFT  = [RSHIFT;LSHIFT];

RL = [zeros(size(RSHIFT,1),1); ones(size(RSHIFT,1),1)];

close all;
fh1 = figure('Units','normalized','OuterPosition',[.01 .06 .35 .85],'Color','w');
ax1 = axes('Position',[.09 .09 .8 .8],'Color','none');

Y = RLSHIFT';
G = {RL',[AD;AD]',[BRAAK;BRAAK]'};

    axes(ax1)
ph1=boxplot(Y,G,...
    'OutlierSize',1,'Symbol','.','Widths',.6,...
    'FullFactors','on','FactorGap',[10,10],...
    'Jitter',0,'Notch','marker','ExtremeMode','compress');
    ax1.YLim=[-.5 .5];


LTAB = table();
RTAB = table();
LTAB.SHIFT = LSHIFT;
RTAB.SHIFT = RSHIFT;
LTAB.AD = AD;
RTAB.AD = AD;
LTAB.BRAAK = BRAAK;
RTAB.BRAAK = BRAAK;

D = RTAB.SHIFT((RTAB.BRAAK==6) & (RTAB.AD==0),:);
D = RTAB.SHIFT((RTAB.BRAAK==6) & (RTAB.AD==1),:);




% statarray = grpstats(T,{'AD','BRAAK'})

% G = string([RL, [AD;AD], [BRAAK;BRAAK]]);
% [p,t,stats] = anova1(Y,G,'off');
% [c,m,h,nms] = multcompare(stats);
% anovan(Y,G,3,3,varnames)
% 
% [~,~,stats] = anovan(Y,G,...
%     'varnames',{'RL','AD','BRAAK'});

% 'DataLim',[-.4 .4]
% 'ExtremeMode','compress'
% 'MedianStyle','target',
% ,'BoxStyle','filled'
% ,'Notch','marker'
% 'Whisker',.7,
% ,'PlotStyle','compact'

%==========================================================================
%% ANDREWS PLOT
%==========================================================================
clc; clearvars -except P ADSP SNP
% load fisheriris
%----------------
iSNP = SNP.iSNP;
iLOC = SNP.iLOC;
%----------------


% ONLY ROWS WITH BRAAK SCORES
bSNP = iSNP(~isnan(iSNP.BRAAK),:);
bSNP.BRAAK(bSNP.BRAAK==0) = 1;
bSNP.BRAAK(bSNP.BRAAK==2) = 1;


mu = mean(bSNP.ALT);
[~,LI] = sort(mu);
RI = fliplr(LI);
Ri  = RI(1:15); 
Li  = LI(1:15);


    close all;
    figure('Units','normalized','Position',[.1 .1 .5 .7],'Color','w');
    axes('Position',[.1 .1 .8 .8],'Color','none');

andrewsplot(bSNP.ALT(:,[Ri Li]),...
    'group',bSNP.BRAAK,'Standardize','PCA','quantile',.51)
    axis([0 1 -.5 .5])




S = SNP.iSNP;
S.BRAAK(isnan(S.BRAAK)) = -1;

GRP = S(:,[2:9 12 18:23]);
% GRP.SEX      =categorical(GRP.SEX);
% GRP.APOE     =categorical(GRP.APOE);
% GRP.AUTOPSY  =categorical(GRP.AUTOPSY);
% GRP.BRAAK    =categorical(GRP.BRAAK);
% GRP.PREVAD   =categorical(GRP.PREVAD);
% GRP.INCAD    =categorical(GRP.INCAD);
% GRP.AD       =categorical(GRP.AD);
% GRP.Consent  =categorical(GRP.Consent);
% GRP.COHORTNUM=categorical(GRP.COHORTNUM);
% GRP.RACE     =categorical(GRP.RACE);
% GRP.ETHNIC   =categorical(GRP.ETHNIC);
% GRP.GOODCOH  =categorical(GRP.GOODCOH);
% GRP.BRAKCOH  =categorical(GRP.BRAKCOH);
grpstats(GRP,{'AD','BRAAK'})



%==========================================================================
%% GENERAL LINEAR MODEL
%==========================================================================
close all; clc; clearvars -except P ADSP SNP LOOP
%----------------
iSNP = SNP.iSNP;
iLOC = SNP.iLOC;

bSNP = iSNP(~isnan(iSNP.BRAAK),:);
%----------------


iSNP.muREF = zscore(mean(iSNP.REF,2));
iSNP.muALT = zscore(mean(iSNP.ALT,2));
iSNP.muDIF = zscore(mean(iSNP.DIF,2));
iSNP.muPVX = zscore(mean(iSNP.PVX,2));
iSNP.muNAT = zscore(mean(iSNP.NAT,2));
% writetable( iSNP ,'/Users/bradleymonk/Desktop/iSNP.csv');



lm = fitlm(iSNP,'AD~muREF')

lm = fitlm(iSNP,'AD~muALT')

lm = fitlm(iSNP,'AD~muDIF')





bSNP.muREF = zscore(mean(bSNP.REF,2));
bSNP.muALT = zscore(mean(bSNP.ALT,2));
bSNP.muDIF = zscore(mean(bSNP.DIF,2));


lm = fitlm(bSNP,'BRAAK~muREF')

lm = fitlm(bSNP,'BRAAK~muALT')

lm = fitlm(bSNP,'BRAAK~muDIF')



%%
% lm = fitlm(iSNP,'AD~muREF+muALT')


[PCC,PCS] = pca(iSNP.ALT,'Options',ss);


scatter(PCS(:,1),PCS(:,2),100,iSNP.AD,'.')
colormap(lines(2)); colorbar


%==========================================================================
%% PCA
%==========================================================================
close all; clc; clearvars -except P ADSP SNP iSNP iLOC LOOP
%----------------
iSNP = SNP.iSNP;
iLOC = SNP.iLOC;

PVX = iSNP.PVX;
REF = iSNP.REF;
ALT = iSNP.ALT;
DIF = iSNP.DIF;
NAT = iSNP.NAT;
%----------------




% PCA
%--------------------------------
ss = statset('pca');
ss.Display = 'none';
ss.MaxIter = 100;
ss.TolFun = 1e4;
ss.TolX = 1e4;
ss.UseParallel = true;
%--------------------------------
[PCC,PCS] = pca(DIF,'Options',ss);

close all; 
fh1 = figure('Units','normalized','Position',[.1 .1 .55 .7],'Color','w');
ax1 = axes('Position',[.1 .1 .8 .8],'Color','none');

scatter(PCS(:,1),PCS(:,2),100,iSNP.AD,'.')
AdvancedColormap('bled'); colorbar
%--------------------------------




%--------------------------------
[ PCC , PCS ,~,~, PCE ] = pca(  DIF' , 'Options',ss);
close all; 
fh1 = figure('Units','normalized','Position',[.1 .1 .55 .7],'Color','w');
ax1 = axes('Position',[.1 .1 .8 .8],'Color','none');
% CM = -log(iLOC.FISHP); CM(CM>10)=10; CM(CM<3)=3;
% CM = log(iLOC.FISHOR); %CM(CM>2)=2;
CM = log(iLOC.srtVAL); CM(CM>1.3)=1.3;
scatter(PCS(:,1),PCS(:,2),700,CM,'.')
cmocean('matter'); colorbar
%--------------------------------






%==========================================================================
%% BOOSTED TREES
%==========================================================================
close all; clc; clearvars -except P ADSP SNP iSNP iLOC LOOP
%----------------
iSNP = SNP.iSNP;
iLOC = SNP.iLOC;
%----------------

[PCC,PCS] = pca(iSNP.DIF,'Options',ss);

BOOSTTREES = fitcensemble(PCS(:,1:5),iSNP.AD, ...
    'Method', 'AdaBoostM1', ...
    'NumLearningCycles', 30, ...
    'Learners', templateTree('MaxNumSplits', 20), ...
    'LearnRate', 0.1, ...
    'ClassNames', [0; 1]);

yfit = BOOSTTREES.predict(PCS(:,1:5));
mean(iSNP.AD == yfit)



%======================================================================
% NEURAL NET
%======================================================================
NN = patternnet([50 20],'trainscg','crossentropy');
NN.trainParam.max_fail = 50;
NN.trainParam.showWindow = 0;
NN.performParam.regularization = 0.1;
NN.performParam.normalization = 'none';

net = train(NN,PCS(:,1:5)',iSNP.AD');  % TRAIN NETQ
[ERR,~,~,~] = confusion(iSNP.AD',net(PCS(:,1:5)'));
net_pctCorrect = 1-ERR






%==========================================================================
%% tSNE
%==========================================================================








close all
Y = tsne([DIF(:,1),DIF(:,2),DIF(:,3)]);
gscatter(Y(:,1),Y(:,2),iSNP.AD)





rng('default') % for reproducibility
Y = tsne(meas,'Algorithm','exact','Distance','mahalanobis');
subplot(2,2,1)
gscatter(Y(:,1),Y(:,2),species)
title('Mahalanobis')

rng('default') % for fair comparison
Y = tsne(meas,'Algorithm','exact','Distance','cosine');
subplot(2,2,2)
gscatter(Y(:,1),Y(:,2),species)
title('Cosine')

rng('default') % for fair comparison
Y = tsne(meas,'Algorithm','exact','Distance','chebychev');
subplot(2,2,3)
gscatter(Y(:,1),Y(:,2),species)
title('Chebychev')

rng('default') % for fair comparison
Y = tsne(meas,'Algorithm','exact','Distance','euclidean');
subplot(2,2,4)
gscatter(Y(:,1),Y(:,2),species)
title('Euclidean')




scatter3(PCdif(:,1),PCdif(:,2),PCdif(:,3),100,iSNP.AD,'.')
colormap(lines(2))
% ,'MarkerSize',5








%==========================================================================
%% DOES iSNP.DIF HAVE DIFFERENTIAL EFFECTS FOR DIFFERENT APOE ALLELS
%==========================================================================
clc; clearvars -except P ADSP SNP iSNP iLOC LOOP



DIF = iSNP.DIF;
ALT = iSNP.ALT;

CTRL = iSNP.DIF(iSNP.AD~=1,:);
CASE = iSNP.DIF(iSNP.AD==1,:);
CTRLmu = mean(CTRL,1);
CASEmu = mean(CASE,1);

srtVAL = iLOC.srtVAL;



quantile(abs(sum(iSNP.DIF)),.25)

sumSNP = sum(iSNP.DIF ,2);

suSNP.b0 = sumSNP(iSNP.BRAAK==0);
suSNP.b1 = sumSNP(iSNP.BRAAK==1);
suSNP.b2 = sumSNP(iSNP.BRAAK==2);
suSNP.b3 = sumSNP(iSNP.BRAAK==3);
suSNP.b4 = sumSNP(iSNP.BRAAK==4);
suSNP.b5 = sumSNP(iSNP.BRAAK==5);
suSNP.b6 = sumSNP(iSNP.BRAAK==6);

Xmu = [mean(suSNP.b0); mean(suSNP.b1); mean(suSNP.b2); mean(suSNP.b3);...
       mean(suSNP.b4); mean(suSNP.b5);  mean(suSNP.b6)];

Xsd = [std(suSNP.b0); std(suSNP.b1); std(suSNP.b2); std(suSNP.b3);...
       std(suSNP.b4); std(suSNP.b5);  std(suSNP.b6)];

Xnn = [numel(suSNP.b0); numel(suSNP.b1); numel(suSNP.b2); numel(suSNP.b3);...
       numel(suSNP.b4); numel(suSNP.b5);  numel(suSNP.b6)];

Xse = Xsd./sqrt(Xnn);


close all
errorbar(Xmu,Xse)





%==========================================================================
%%
%==========================================================================
scatter([0:6],[mean(CACOBRAAKb0) mean(CACOBRAAKb0) mean(CACOBRAAKb0)...
mean(CACOBRAAKb0) mean(CACOBRAAKb0) mean(CACOBRAAKb0) mean(CACOBRAAKb0)])



SV = iLOC.srtVAL';

SVBRAAKb0 = sum(mean(SV(iSNP.BRAAK==0) ,2),2);
SVBRAAKb1 = sum(mean(SV(iSNP.BRAAK==1) ,2),2);
SVBRAAKb2 = sum(mean(SV(iSNP.BRAAK==2) ,2),2);
SVBRAAKb3 = sum(mean(SV(iSNP.BRAAK==3) ,2),2);
SVBRAAKb4 = sum(mean(SV(iSNP.BRAAK==4) ,2),2);
SVBRAAKb5 = sum(mean(SV(iSNP.BRAAK==5) ,2),2);
SVBRAAKb6 = sum(mean(SV(iSNP.BRAAK==6) ,2),2);





%################   SIX PACK   ################
close all
fh01 = figure('Units','normalized','OuterPosition',[.01 .05 .95 .90],...
              'Color','w','MenuBar','none');
ax01 = axes('Position',[.04 .56 .25 .35],'Color','none');
ax02 = axes('Position',[.35 .56 .25 .35],'Color','none');
ax03 = axes('Position',[.70 .56 .25 .35],'Color','none');
ax04 = axes('Position',[.04 .06 .25 .35],'Color','none');
ax05 = axes('Position',[.35 .06 .25 .35],'Color','none');
ax06 = axes('Position',[.70 .06 .25 .35],'Color','none');


    axes(ax01)
scatter(CTRLAPOE22,CASEAPOE22,[],(CTRLAPOE22-CASEAPOE22),'filled')
    axis([-.4 .4 -.4 .4])
    hold on; line([-.4 .4], [-.4 .4])
    xlabel('\fontsize{16} ?   ? CTRL')
    ylabel('\fontsize{16} ?   ? CTRL')
    title('APOE 22');








%==========================================================================
%% PLOT AGE BY MEAN REF VALUE, AND COLOR BY APOE STATUS
%==========================================================================
close all; clc; clearvars -except P ADSP SNP
%----------------
iSNP = SNP.iSNP;
iLOC = SNP.iLOC;
%----------------
PVX = iSNP.PVX;
REF = iSNP.REF;
ALT = iSNP.ALT;
DIF = iSNP.DIF;
NAT = iSNP.NAT;
%----------------



% LOAD RUN#3 SNP.mat INTO YOUR WORKSPACE (IT CONTAINS iSNP & iLOC)
SNP = load(['SNP.mat'],'iSNP','iLOC');


% PULL OUT iSNP
iSNP = SNP.iSNP;


% REMOVE PEOPLE UNDER 60 YEARS OLD
iSNP(iSNP.AGE<60,:) = [];


% PUT CASES & CONTROLS INTO THEIR OWN iSNP CONTAINER
iSNPca = iSNP(iSNP.AD==1,:);
iSNPco = iSNP(iSNP.AD~=1,:);



% GET THE AVERAGE CV PER PERSON
CaseCV = mean(iSNPca.REF,2);
CtrlCV = mean(iSNPco.REF,2);


% PLOT AGE X CV (COLOR BY APOE)
subplot(1,2,1); 
gscatter( iSNPca.AGE, CaseCV, iSNPca.APOE , lines(6), '.', 15)
subplot(1,2,2); 
gscatter( iSNPco.AGE, CtrlCV, iSNPco.APOE , lines(6), '.', 15)







%%






