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
P.home  = fileparts(which('GENOS.m')); cd(P.home);
P.funs  = [P.home filesep 'genos_functions'];
P.mfuns = [P.funs filesep 'genos_main_functions'];
P.other = [P.home filesep 'genos_other'];
P.data  = [P.home filesep 'genos_data'];
P.mod   = [P.data filesep 'modSNP'];
% P.data = 'F:\GENOSDATA';
addpath(join(string(struct2cell(P)),pathsep,1))
cd(P.home); P.f = filesep;



tic
ADSP = load('GENOSDATA.mat');
toc



clearvars -except P ADSP
%==========================================================================
%% IMPORT modSNP DATA
%==========================================================================
close all; clear; clc;



P.BASE = '/Users/bradleymonk/Documents/MATLAB/GIT/genomics/genos/genos_data/modSNP/';

P.RANDRUN = [P.BASE 'RUN5_RANDTARGS/DAT/SNP.mat'];
P.TARGRUN = [P.BASE 'RUN9_PLO/DAT/SNP.mat'];
% P.TARGRUN = [P.BASE 'RUN6_INCLUDE_APOE/DAT/SNP.mat'];




RSNP = load(  P.RANDRUN  ,'iSNP','iLOC');
TSNP = load(  P.TARGRUN  ,'iSNP','iLOC');





whos('TSNP'); fieldnames(TSNP)
clearvars -except P TSNP RSNP
%==========================================================================
%% COMPUTE MEAN SHIFT PER SNP & PER PERSON
%==========================================================================
clc; close all; clearvars -except P TSNP RSNP
%----------------
tSNP = TSNP.iSNP;
rSNP = RSNP.iSNP;
tLOC = TSNP.iLOC;
rLOC = RSNP.iLOC;
%----------------




TSNP.iLOC.muPVX = nanmean(tSNP.PVX,1)';
TSNP.iLOC.muREF = nanmean(tSNP.REF,1)';
TSNP.iLOC.muALT = nanmean(tSNP.ALT,1)';
TSNP.iLOC.muDIF = nanmean(tSNP.DIF,1)';
TSNP.iLOC.muNAT = nanmean(tSNP.NAT,1)';
% writetable(tLOC,'/Users/bradleymonk/Desktop/tTOPSNPS.csv')



RSNP.iLOC.muPVX = mean(rSNP.PVX,1)';
RSNP.iLOC.muREF = mean(rSNP.REF,1)';
RSNP.iLOC.muALT = mean(rSNP.ALT,1)';
RSNP.iLOC.muDIF = mean(rSNP.DIF,1)';
% writetable(rLOC,'/Users/bradleymonk/Desktop/rTOPSNPS.csv')



TSNP.iSNP.muPVX = nanmean(tSNP.PVX,2);
TSNP.iSNP.muREF = nanmean(tSNP.REF,2);
TSNP.iSNP.muALT = nanmean(tSNP.ALT,2);
TSNP.iSNP.muDIF = nanmean(tSNP.DIF,2);
TSNP.iSNP.muNAT = nanmean(tSNP.NAT,2);
% writetable(tLOC,'/Users/bradleymonk/Desktop/tTOPSNPS.csv')



RSNP.iSNP.muPVX = mean(rSNP.PVX,2);
RSNP.iSNP.muREF = mean(rSNP.REF,2);
RSNP.iSNP.muALT = mean(rSNP.ALT,2);
RSNP.iSNP.muDIF = mean(rSNP.DIF,2);
% writetable(rLOC,'/Users/bradleymonk/Desktop/rTOPSNPS.csv')



%==========================================================================
%% PLOT HISTOGRAM OF APOE STATUS
%==========================================================================
clc; close all; clearvars -except P TSNP RSNP
%----------------
tSNP = TSNP.iSNP;
rSNP = RSNP.iSNP;
tLOC = TSNP.iLOC;
rLOC = RSNP.iLOC;
%----------------



% THEORETICAL APOE SUBTYPE CV
%--------------------------------------------------------------


BINS = -.5:.01:.5;

close all
fh01 = figure('Units','pixels','Position',[10 35 900 500],'Color','w'); 
ax01 = axes('Position',[.1 .1 .8 .8],'Color','none');
set(groot,'defaultAxesFontSize',24)

histogram(mean(tSNP.ALT(tSNP.APOE == 24,:)),BINS)
hold on
histogram(mean(tSNP.ALT(tSNP.APOE == 33,:)),BINS)
hold on
histogram(mean(tSNP.ALT(tSNP.APOE == 34,:)),BINS)
hold on
histogram(mean(tSNP.ALT(tSNP.APOE == 23,:)),BINS)
hold on
histogram(mean(tSNP.ALT(tSNP.APOE == 44,:)),BINS)
hold on
histogram(mean(tSNP.ALT(tSNP.APOE == 22,:)),BINS)
% legend({'24','33','34','23','44','22',})
xlim([-.5 .5]); box off


[E22,EDG22] = histcounts(mean(tSNP.ALT(tSNP.APOE==22,:),2),BINS);
[E23,EDG23] = histcounts(mean(tSNP.ALT(tSNP.APOE==23,:),2),BINS);
[E24,EDG24] = histcounts(mean(tSNP.ALT(tSNP.APOE==24,:),2),BINS);
[E33,EDG33] = histcounts(mean(tSNP.ALT(tSNP.APOE==33,:),2),BINS);
[E34,EDG34] = histcounts(mean(tSNP.ALT(tSNP.APOE==34,:),2),BINS);
[E44,EDG44] = histcounts(mean(tSNP.ALT(tSNP.APOE==44,:),2),BINS);





%% REAL APOE SUBTYPE CV
%--------------------------------------------------------------



e22 = -.5:.025:.5;
e23 = -.5:.025:.5;
e24 = -.5:.025:.5;
e33 = -.5:.025:.5;
e34 = -.5:.025:.5;
e44 = -.5:.025:.5;

EDGES = {e22,e23,e24,e33,e34,e44};

CENTERS = {
(e22(1:end-1)+e22(2:end))./2,...
(e23(1:end-1)+e23(2:end))./2,...
(e24(1:end-1)+e24(2:end))./2,...
(e33(1:end-1)+e33(2:end))./2,...
(e34(1:end-1)+e34(2:end))./2,...
(e44(1:end-1)+e44(2:end))./2,...
};


% 'Normalization','count'
% 'Normalization','probability'
% 'Normalization','countdensity'
% 'Normalization','pdf'
% 'Normalization','cumcount'
% 'Normalization','cdf'


[N22,edges] = histcounts(mean(tSNP.REF(tSNP.APOE==22,:),2),EDGES{1},'Normalization','probability');
[N23,edges] = histcounts(mean(tSNP.REF(tSNP.APOE==23,:),2),EDGES{2},'Normalization','probability');
[N24,edges] = histcounts(mean(tSNP.REF(tSNP.APOE==24,:),2),EDGES{3},'Normalization','probability');
[N33,edges] = histcounts(mean(tSNP.REF(tSNP.APOE==33,:),2),EDGES{4},'Normalization','probability');
[N34,edges] = histcounts(mean(tSNP.REF(tSNP.APOE==34,:),2),EDGES{5},'Normalization','probability');
[N44,edges] = histcounts(mean(tSNP.REF(tSNP.APOE==44,:),2),EDGES{6},'Normalization','probability');






fh02 = figure('Units','pixels','Position',[50 90 900 500],'Color','w'); 
ax02 = axes('Position',[.1 .1 .8 .8],'Color','none');
set(groot,'defaultAxesFontSize',24)

cm=lines(6); colormap(cm)
area(CENTERS{3},N33,'FaceAlpha',.5,'LineWidth',5,'EdgeColor',[.4 .4 .4],'FaceColor',[.4 .4 .4]); hold on
area(CENTERS{4},N24,'FaceAlpha',.5,'LineWidth',5,'EdgeColor',cm(2,:),'FaceColor',cm(2,:)); hold on
area(CENTERS{5},N34,'FaceAlpha',.5,'LineWidth',5,'EdgeColor',cm(3,:),'FaceColor',cm(3,:)); hold on
area(CENTERS{2},N23,'FaceAlpha',.5,'LineWidth',5,'EdgeColor',cm(4,:),'FaceColor',cm(4,:)); hold on
area(CENTERS{6},N44,'FaceAlpha',.5,'LineWidth',5,'EdgeColor',[.1 .9 .4],'FaceColor',[.1 .9 .4]); hold on
area(CENTERS{1},N22,'FaceAlpha',.5,'LineWidth',5,'EdgeColor',cm(6,:),'FaceColor',cm(6,:)); hold on
% legend({'24','33','34','23','44','22',})
xlim([-.5 .5]); box off




%==========================================================================
%% COMPUTE ECDF_RANDTARGS & PVAL_TOPTARGS
%==========================================================================
clc; close all; clearvars -except P TSNP RSNP
%----------------
tSNP = TSNP.iSNP;
rSNP = RSNP.iSNP;
tLOC = TSNP.iLOC;
rLOC = RSNP.iLOC;
%----------------



% GET MEAN CV PER SNP
tREF = mean(tSNP.REF,1)';
tALT = mean(tSNP.ALT,1)';
tDIF = mean(tSNP.DIF,1)';

rREF = mean(rSNP.REF,1)';
rALT = mean(rSNP.ALT,1)';
rDIF = mean(rSNP.DIF,1)';



% MAKE RAND eCDF
[f,x] = ecdf(rREF);
x = x + rand(numel(x),1).*1e-10;
rCDFref = makedist('PiecewiseLinear',unique(x,'sorted'),f);

[f,x] = ecdf(rALT);
x = x + rand(numel(x),1).*1e-10;
rCDFalt = makedist('PiecewiseLinear',unique(x,'sorted'),f);

[f,x] = ecdf(rDIF);
x = x + rand(numel(x),1).*1e-10;
rCDFdif = makedist('PiecewiseLinear',unique(x,'sorted'),f);



% DETERMINE P-VALUES OF targDIF FROM eCDF OF randDIF
pCDFref = cdf(rCDFref,tREF);
pCDFalt = cdf(rCDFalt,tALT);
pCDFdif = cdf(rCDFdif,tDIF);

pCDFref(pCDFref>.5) = 1 - pCDFref(pCDFref>.5);
pCDFalt(pCDFalt>.5) = 1 - pCDFalt(pCDFalt>.5);
pCDFdif(pCDFdif>.5) = 1 - pCDFdif(pCDFdif>.5);


TSNP.iLOC.pCDFref = pCDFref;
TSNP.iLOC.pCDFalt = pCDFalt;
TSNP.iLOC.pCDFdif = pCDFdif;





%==========================================================================
%% DETERMINE WHETHER modSNP SHIFT PREDICTS BRAAK
%==========================================================================
clc; close all; clearvars -except P TSNP RSNP
%----------------
tSNP = TSNP.iSNP;
rSNP = RSNP.iSNP;
tLOC = TSNP.iLOC;
rLOC = RSNP.iLOC;
%----------------








% IS MEAN SNP PVAL BELOW ECDF ALPHA (FOR CASE TESTING)
% r = (tLOC.pCDFref<.04) & (tLOC.TRFISHOR<1); sum(r)
% a = (tLOC.pCDFalt<.01) & (tLOC.TRFISHOR<1); sum(a)
% d = (tLOC.pCDFdif<.01) & (tLOC.TRFISHOR<1); sum(d)

% IS MEAN SNP PVAL BELOW ECDF ALPHA (FOR CTRL TESTING)
% r = (tLOC.pCDFref<.04) & (tLOC.TRFISHOR>1); sum(r)
% a = (tLOC.pCDFalt<.01) & (tLOC.TRFISHOR>1); sum(a)
% d = (tLOC.pCDFdif<.01) & (tLOC.TRFISHOR>1); sum(d)


% USE FOR RANDOM RUN
% r = (tLOC.pCDFref>.04); sum(r)
% a = (tLOC.pCDFalt>.01); sum(a)
% d = (tLOC.pCDFdif>.01); sum(d)




% IS MEAN SNP PVAL BELOW eCDF ALPHA CUTOFF
r = (tLOC.pCDFref<.01); sum(r)
a = (tLOC.pCDFalt<.01); sum(a)
d = (tLOC.pCDFdif<.01); sum(d)


% GET PER-PERSON GENOTYPES FOR SNPs WITH CVs BELOW eCDF CUTOFF
tSNP.rPVX = tSNP.PVX(:,r);
tSNP.aPVX = tSNP.PVX(:,a);
tSNP.dPVX = tSNP.PVX(:,d);


% GET R.A.D. CVs FOR SNPs WITH CVs BELOW eCDF CUTOFF
tSNP.rREF = tSNP.REF(:,r);
tSNP.aALT = tSNP.ALT(:,a);
tSNP.dDIF = tSNP.DIF(:,d);


% INSERT NaN INTO R.A.D. CVs BASED ON PERSON GENOTYPE
mM = (tSNP.rPVX<0) .* tSNP.rREF; mM(mM==0) = NaN;
mR = (tSNP.rPVX>0) .* tSNP.rREF; mR(mR==0) = NaN;
mA = (tSNP.aPVX>0) .* tSNP.aALT; mA(mA==0) = NaN;
mD = (tSNP.dPVX>0) .* tSNP.dDIF; mD(mD==0) = NaN;


% GET MEAN CV PER-PERSON, FOR LOW-P SNPs, FOR SNPs A PERSON HAS ALT GENOTYPE
tSNP.mM = nanmean( mM ,2);
tSNP.mR = nanmean( mR ,2);
tSNP.mA = nanmean( mA ,2);
tSNP.mD = nanmean( mD ,2);


close all;
fh01 = figure('Units','pixels','Position',[10 35 900 500],'Color','w'); 
ax01 = axes('Position',[.06 .06 .9 .9],'Color','none');
subplot(1,4,1); boxplot(tSNP.mM,tSNP.BRAAK,'PlotStyle','compact','Whisker',10); ylim([-.5 .5]); title('G1 REAL REF');
subplot(1,4,2); boxplot(tSNP.mR,tSNP.BRAAK,'PlotStyle','compact','Whisker',10); ylim([-.5 .5]); title('G2 FAUX REF');
subplot(1,4,3); boxplot(tSNP.mA,tSNP.BRAAK,'PlotStyle','compact','Whisker',10); ylim([-.5 .5]); title('G2 REAL ALT')
subplot(1,4,4); boxplot(tSNP.mD,tSNP.BRAAK,'PlotStyle','compact','Whisker',10); ylim([-.3 .3]); title('G2 DIF')




% muM = nanmean( mM ,2);
% muR = nanmean( mR ,2);
% muA = nanmean( mA ,2);
% muD = nanmean( mD ,2);

TSNP.iSNP = tSNP;

%==========================================================================
%% LINEAR REGRESSION
%==========================================================================
clc; close all; clearvars -except P TSNP RSNP
%----------------
tSNP = TSNP.iSNP;
rSNP = RSNP.iSNP;
tLOC = TSNP.iLOC;
rLOC = RSNP.iLOC;
%----------------




ai = ~isnan(tSNP.BRAAK);
bi = ~isnan(tSNP.mR);
ci = ~isnan(tSNP.mA);
di = ~isnan(tSNP.mD);
T = tSNP(ai&bi&ci&di,:);


% BRAAK = T.BRAAK - nanmean(T.BRAAK);
% AD = T.AD - nanmean(T.AD);
% CVref = T.mR - nanmean(T.mR);
% CValt = T.mA - nanmean(T.mA);
% CVdif = T.mD - nanmean(T.mD);


BRAAK = T.BRAAK;
AD = zscore(T.AD);
zCVref = zscore(T.mR);
zCValt = zscore(T.mA);
zCVdif = zscore(T.mD);


TBL = table(BRAAK,AD,zCVref,zCValt,zCVdif);




clc
disp('-------------------------------------------------------------------')
lm_mR = fitlm(TBL,'BRAAK~zCVref')
disp('-------------------------------------------------------------------')
lm_mA = fitlm(TBL,'BRAAK~zCValt')
disp('-------------------------------------------------------------------')
lm_mD = fitlm(TBL,'BRAAK~zCVdif')
disp('-------------------------------------------------------------------')



% lm = fitlm(tSNP,'BRAAK~mD+AD')




%%

tSNP(tSNP.AGE<60,:) = [];

%==========================================================================
%% GENERATE BOXPLOT AND ERRORBAR PLOTS FOR CV-X-AGE
%==========================================================================


N=8;

[y,Y] = discretize(tSNP.AGE,N)

x = Y(1:N);

%--- REF
R = tSNP.mR .* (y==(1:N)); R(R==0) = NaN;
Rmu = nanmean(R); Rme = nanmedian(R); Rse = nanstd(R) ./ sqrt(numel(nanstd(R)));

%--- ALT
A = tSNP.mA .* (y==(1:N)); A(A==0) = NaN;
Amu = nanmean(A); Ame = nanmedian(A); Ase = nanstd(A) ./ sqrt(numel(nanstd(A)));

%--- DIF
D = tSNP.mD .* (y==(1:N)); D(D==0) = NaN;
Dmu = nanmean(D); Dme = nanmedian(D); Dse = nanstd(D) ./ sqrt(numel(nanstd(D)));


close all; figure('Units','normalized','Position',[.1 .08 .4 .78])
ax1=subplot(2,3,1); boxplot(tSNP.mR,Y(y),'PlotStyle','compact','Whisker',100); axis([0 7 -.5 .5]); title('REF'); 
ax2=subplot(2,3,2); boxplot(tSNP.mA,Y(y),'PlotStyle','compact','Whisker',100); axis([0 7 -.5 .5]); title('ALT');
ax3=subplot(2,3,3); boxplot(tSNP.mD,Y(y),'PlotStyle','compact','Whisker',100); axis([0 7 -.2 .2]); title('DIF');


ax4=subplot(2,3,4); errorbar(x,Rmu,Rse,'.','MarkerSize',40); axis([55 90 -.3 .31]); title('REF'); 
ax5=subplot(2,3,5); errorbar(x,Amu,Ase,'.','MarkerSize',40); axis([55 90 -.3 .31]); title('ALT');
ax6=subplot(2,3,6); errorbar(x,Dmu,Dse,'.','MarkerSize',40); axis([55 90 -.1 .11]); title('DIF');



%==========================================================================
%% HISTOGRAM OF SHIFT VALUES FROM PEOPLE WHO ACTUALLY HAVE ALT ALLELE
%==========================================================================
clc; close all; clearvars -except P TSNP RSNP
%----------------
tSNP = TSNP.iSNP;
rSNP = RSNP.iSNP;
tLOC = TSNP.iLOC;
rLOC = RSNP.iLOC;
%----------------


clc; close all; figure
histogram(tSNP.mD, (-.5:.01:.5) ,'FaceColor',[.9 .1 .1])
hold on
histogram(tSNP.mA, (-.5:.01:.5) ,'FaceColor',[.01 .21 .99])
hold on
histogram(tSNP.mR, (-.5:.01:.5) ,'FaceColor',[.3 .9 .3],'FaceAlpha',.4)
legend({'DIF','ALT','REF'})



%==========================================================================
%% BOXPLOT OF SHIFT VALUES FROM PEOPLE WHO ACTUALLY HAVE ALT ALLELE
%==========================================================================
clc; close all; clearvars -except P TSNP RSNP
%----------------
tSNP = TSNP.iSNP;
rSNP = RSNP.iSNP;
tLOC = TSNP.iLOC;
rLOC = RSNP.iLOC;
%----------------




clc; close all;
fh1 = figure('Units','normalized','OuterPosition',[.01 .06 .4 .8],'Color','w');
ax1 = axes('Position',[.09 .10 .4 .8],'Color','none');
ax2 = axes('Position',[.57 .10 .4 .8],'Color','none');



    axes(ax1); 
ph1=boxplot(tSNP.mA,tSNP.BRAAK,...
    'OutlierSize',1,'Symbol','.','Widths',.6,...
    'FullFactors','on','FactorGap',[10],...
    'Jitter',0,'Notch','marker','ExtremeMode','compress');
    ax1.YLim=[-.5 .5];
    title('BRAAK vs ALT')


    axes(ax2); 
ph2=boxplot(tSNP.mD,tSNP.BRAAK,...
    'OutlierSize',1,'Symbol','.','Widths',.6,...
    'FullFactors','on','FactorGap',[10],...
    'Jitter',0,'Notch','marker','ExtremeMode','compress');
    ax2.YLim=[-.3 .3];
    title('BRAAK vs DIF')









%%
close all; clc; clearvars -except P ADSP TSNP RSNP
%----------------
tSNP = TSNP.iSNP;
tLOC = TSNP.iLOC;
%----------------


i = ~isnan(tSNP.BRAAK);
A = iosr.statistics.tab2box(tSNP.BRAAK(i),tSNP.mA(i));
D = iosr.statistics.tab2box(tSNP.BRAAK(i),tSNP.mD(i));



clc; close all;
fh1 = figure('Units','normalized','OuterPosition',[.01 .06 .4 .8],'Color','w');
ax1 = axes('Position',[.09 .10 .4 .8],'Color','none');
ax2 = axes('Position',[.57 .10 .4 .8],'Color','none');


    axes(ax1)
iosr.statistics.boxPlot(A,'boxColor',[.3 .7 .9],'lineWidth',1.8,...
    'medianColor',[.3 .3 .3],'notchLine',true,'boxAlpha',.9,...
    'outlierSize',50,'symbolMarker','.','symbolColor',[.2 .5 .7],'showViolin',false,...
    'showScatter',true,'scatterAlpha',.6,'scatterMarker','.','scatterColor',[.3 .5 .6])

    xlabel('BRAAK'); ylabel('ALT CV');
    axis([0 8 -.5 .5])
    ax1.XTickLabel=0:6;


    axes(ax2)
iosr.statistics.boxPlot(D,'boxColor',[.3 .7 .9],'lineWidth',1.8,...
    'medianColor',[.3 .3 .3],'notchLine',true,'boxAlpha',.9,...
    'outlierSize',50,'symbolMarker','.','symbolColor',[.2 .5 .7],'showViolin',false,...
    'showScatter',true,'scatterAlpha',.6,'scatterMarker','.','scatterColor',[.3 .5 .6])
    xlabel('BRAAK'); ylabel('DIF CV');
    axis([0 8 -.3 .31])
    ax2.XTickLabel=0:6;







%==========================================================================
%% GRAMM TOOLBOX PLOTS
%==========================================================================
%{
%        addPrctiles = []             % Additional percentiles to plot.
%        addPrctilesColors            % Colors for the additional percentile markers.
%        addPrctilesLabels = {}       % Labels for additional percentiles.
%        addPrctilesMarkers = {}      % Markers for additional percentiles.
%        addPrctilesSize              % Size of the additional percentile markers.
%        addPrctilesTxtSize           % Size of the additional percentile labels text.
%        boxAlpha = 1                 % The transparency of the boxes.
%        boxColor = 'none'            % Fill color of the boxes.
%        boxWidth = 'auto'            % The width of the boxes.
%        groupLabelFontSize = 9       % The font size of the group labels.
%        groupLabelHeight = 'auto'    % The height of the area reserved for group labels.
%        groupLabels = []             % Labels used to label the boxes for each x-group.
%        groupWidth = 0.75            % The proportion of the x-axis interval across which each x-group of boxes should be spread.
%        limit = '1.5IQR'             % Mode indicating the limits that define outliers.
%        lineColor = 'k'              % Color of the box outlines and whiskers.
%        lineStyle = '-'              % Style of the whisker lines.
%        lineWidth = 1                % Width, in points, of the box outline, whisker lines, notch line, and outlier marker edges.
%        meanColor = 'auto'           % Color of the mean marker.
%        meanMarker = '+'             % Marker used for the mean.
%        meanSize = 6                 % Size, in point units, of the mean markers.
%        medianColor = 'auto'         % Color of the median line.
%        method = 'R-8'               % The method used to calculate the quantiles.
%        notch = false                % Whether the box should have a notch.
%        notchDepth = 0.4             % Depth of the notch as a proportion of half the box width.
%        notchLineColor = 'k'         % Color of the notch line.
%        notchLineStyle = ':'         % Line style of the notch line.
%        notchLine = false            % Whether to draw a horizontal line in the box at the extremes of the notch.
%        outlierSize = 6^2            % Size, in square points, of the outlier markers.
%        percentile = [25,75]         % Percentile limits of the box.
%        sampleFontSize = 9           % Specify the font size of the sample size display.
%        sampleSize = false           % Whether to display the sample size in the top-right of each box.
%        scaleWidth = false           % Scale the width of each box according to the square root of the sample size.
%        scatterAlpha = 1             % The transparency of the scatter markers.
%        scatterColor = [.5 .5 .5]    % Scatter marker color for scatter plots of underlying data.
%        scatterLayer = 'top'         % The layer of scatter plots with respect to the boxes.
%        scatterMarker = 'x'          % Marker used for scatter plots of underlying data.
%        scatterSize = 6^2            % Size, in square points, of the scatter markers.
%        showLegend = false           % Draw a legend.
%        showMean = false             % Display the mean of the data for each box.
%        showOutliers = true          % Display outliers.
%        showScatter = false          % Display a scatter plot of the underlying data for each box.
%        showViolin = false           % Display a violin (kernel density) plot for the underlying data.
%        style = 'normal'             % Determine whether to show additional x-axis labels for the data.
%        symbolColor = 'auto'         % Outlier marker color.
%        symbolMarker = 'o'           % Marker used to denote outliers.
%        theme = 'default'            % Control a range of display properties.
%        themeColors = 'auto'         % Colors used when creating the theme.
%        violinBins = 'auto'          % The bins used to calculate the violins.
%        violinBinWidth = 'auto'      % The width of the bins used to calculate the violins.
%        violinAlpha = 1              % Alpha of the violins.
%        violinColor = 'none'         % Color of the violins.
%        violinWidth = 'auto'         % The width of the violins.
%        violinKernel = 'normal'      % The violin kernel used to calculate the kernel density.
%        xSeparator = false           % Add a separator line between x groups.
%        xSpacing = 'x'               % Determine the x-axis spacing of boxes.
%}

% 'scatterAlpha',.4,'scatterMarker','.','scatterSize',10

%{
% g(1,1)=gramm('x',tSNP.BRAAK,'y',tSNP.mA,'color',tSNP.AD);
% g(1,1).set_names('x','BRAAK','y','ALT','color','AD');
g(1,1)=gramm('x',tSNP.BRAAK,'y',tSNP.mA);
g(1,1).set_names('x','BRAAK','y','ALT');
g(1,1).stat_boxplot('width',0.5,'dodge',.8,'notch',true);
g(1,1).set_color_options('map','d3_10');

g(1,2)=gramm('x',tSNP.BRAAK,'y',tSNP.mD);
g(1,2).set_names('x','BRAAK','y','DIF');
g(1,2).stat_boxplot('width',0.5,'dodge',.8,'notch',true);
g(1,2).set_color_options('map','d3_10');

figure('Position',[100 100 800 800]);
g.draw();



C = load('cars.mat');
cars = C.CARS;

g(1,1)=gramm('x',cars.ORG,'y',cars.HP,'color',cars.CYL,'subset',cars.CYL~=3 & cars.CYL~=5);
g(1,1).set_names('x','ORG','y','HP','color','# CYL');
g(1,2)=copy(g(1,1));
g(1,3)=copy(g(1,1));
g(2,1)=copy(g(1,1));
g(2,2)=copy(g(1,1));



g(2,3)=gramm('x',cars.ORG,'y',cars.HP,'color',cars.USA,'subset',cars.CYL~=3 & cars.CYL~=5);
g(2,3).set_names('x','ORG','y','HP','color','USA');
g(2,3).stat_violin('normalization','area','dodge',0,'fill','edge');
g(2,3).stat_boxplot('width',0.15);
g(2,3).set_title('with stat_boxplot()');
g(2,3).set_color_options('map','brewer_dark');

g.set_title('Options for stat_violin()');
figure('Position',[100 100 800 600]);
g.draw();



close all


%Create data
x=repmat(1:10,1,100);
catx=repmat({'A' 'B' 'C' 'F' 'E' 'D' 'G' 'H' 'I' 'J'},1,100);
y=randn(1,1000)*3;
c=repmat([1 1 1 1 1 1 1 1 1 1    1 1 2 2 2 2 2 3 2 2     2 1 1 2 2 2 2 3 2 3     2 1 1 2 2 2 2 2 2 2],1,25);
y=2+y+x+c*0.5;


clear g
g(1,1)=gramm('x',catx,'y',y,'color',c);
g(2,1)=copy(g(1));
g(3,1)=copy(g(1));
g(4,1)=copy(g(1));
g(5,1)=copy(g(1));


g(1,1).stat_boxplot();
g(1,1).geom_vline('xintercept',0.5:1:10.5,'style','k-');
g(1,1).set_title('''width'',0.6,''dodge'',0.7 (Default)');

g(2,1).stat_boxplot('width',0.5,'dodge',0);
g(2,1).geom_vline('xintercept',0.5:1:10.5,'style','k-');
g(2,1).set_title('''width'',0.5,''dodge'',0');

g(3,1).stat_boxplot('width',1,'dodge',1);
g(3,1).geom_vline('xintercept',0.5:1:10.5,'style','k-');
g(3,1).set_title('''width'',1,''dodge'',1');

g(4,1).stat_boxplot('width',0.6,'dodge',0.4);
g(4,1).geom_vline('xintercept',0.5:1:10.5,'style','k-');
g(4,1).set_title('''width'',0.6,''dodge'',0.4');

g(5,1).facet_grid([],c);
g(5,1).stat_boxplot('width',0.5,'dodge',0,'notch',true);
g(5,1).set_title('''width'',0.5,''dodge'',0,''notch'',true');

g.set_title('Dodge and spacing options for stat_boxplot()');

figure('Position',[100 100 800 1000]);
g.draw();
%}











%==========================================================================
%% COMPUTE CV DISTRIBUTIONS FOR EACH APOE SUBTYPE
%==========================================================================
close all; clc; clearvars -except P ADSP TSNP RSNP
%----------------
tSNP = TSNP.iSNP;
tLOC = TSNP.iLOC;
%----------------



EDGES = -.5:.01:.5;

close all
fh01 = figure('Units','pixels','Position',[10 35 900 500],'Color','w'); 
ax01 = axes('Position',[.1 .1 .8 .8],'Color','none');
set(groot,'defaultAxesFontSize',24)

histogram(mean(tSNP.ALT(tSNP.APOE == 24,:)),EDGES)
hold on
histogram(mean(tSNP.ALT(tSNP.APOE == 33,:)),EDGES)
hold on
histogram(mean(tSNP.ALT(tSNP.APOE == 34,:)),EDGES)
hold on
histogram(mean(tSNP.ALT(tSNP.APOE == 23,:)),EDGES)
hold on
histogram(mean(tSNP.ALT(tSNP.APOE == 44,:)),EDGES)
hold on
histogram(mean(tSNP.ALT(tSNP.APOE == 22,:)),EDGES)
% legend({'24','33','34','23','44','22',})
xlim([-.5 .5]); box off










%==========================================================================
%% GET MEAN BRAAK SCORE FOR NAT AFTER DETERMINING abs(DIF) > .13
%==========================================================================
close all; clc; clearvars -except P ADSP SNP
%----------------
iSNP = SNP.iSNP;
iLOC = SNP.iLOC;
%----------------


% GET AVERAGE DIF SHIFT PER SNP
muDIF = mean(iSNP.DIF);


% GET INDEX FOR SIGNIFICANT SNPS
hi = muDIF >  .13;
lo = muDIF < -.13;


% COUNT NUMBER OF SNPS THAT PASSED SIGNIFICANCE
nHI = sum(hi)
nLO = sum(lo)



% CREATE TABLE THAT INCLUDES ONLY COLUMNS FROM SIG SNPS 
hiSNP = iSNP(:,1:23);
hiSNP.PVX = iSNP.PVX(:,hi) > 0;
hiSNP.REF = iSNP.REF(:,hi);
hiSNP.ALT = iSNP.ALT(:,hi);
hiSNP.DIF = iSNP.DIF(:,hi);


loSNP = iSNP(:,1:23);
loSNP.PVX = iSNP.PVX(:,lo) > 0;
loSNP.REF = iSNP.REF(:,lo);
loSNP.ALT = iSNP.ALT(:,lo);
loSNP.DIF = iSNP.DIF(:,lo);



hiBRAAK = repmat(hiSNP.BRAAK,1,nHI);
loBRAAK = repmat(loSNP.BRAAK,1,nLO);

shiBRAAK = hiBRAAK(hiSNP.PVX); 
sloBRAAK = loBRAAK(loSNP.PVX);



nanmean(shiBRAAK)

nanmean(sloBRAAK)

BG = [shiBRAAK   ones( size(shiBRAAK,1),1);...
     sloBRAAK  zeros( size(sloBRAAK,1),1)];


close all;
fh1 = figure('Units','normalized','OuterPosition',[.01 .06 .35 .85],'Color','w');
ax1 = axes('Position',[.09 .09 .8 .8],'Color','none');

    axes(ax1)
ph1=boxplot(BG(:,1),BG(:,2),...
    'OutlierSize',1,'Symbol','.','Widths',.6,...
    'FullFactors','on','FactorGap',[10],...
    'Jitter',0,'Notch','marker','ExtremeMode','compress');
    %ax1.YLim=[-.5 .5];




%==========================================================================
%% EXAMINE CV
%==========================================================================
close all; clc; clearvars -except P ADSP SNP
%----------------
iSNP = SNP.iSNP;
iLOC = SNP.iLOC;
%----------------


close all;
plot(iSNP.REF(1,:))
xlim([-.5 .5])




CASEiSNP = iSNP(iSNP.AD==1,:);
CTRLiSNP = iSNP(iSNP.AD~=1,:);



close all;
subplot(1,2,1); 
plot(mean(CASEiSNP.REF))
hold on
plot(mean(CTRLiSNP.REF))
ylim([-.5 .5])

subplot(1,2,2); 
histogram(mean(CASEiSNP.REF))
hold on;
histogram(mean(CTRLiSNP.REF))






%==========================================================================
%% PLOT A HISTOGRAM OF ALL VALUES
%==========================================================================
close all; clc; clearvars -except P ADSP SNP
%----------------
iSNP = SNP.iSNP;
iLOC = SNP.iLOC;
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
xlim([-.5 .5]); title('REF')

axes(ax02)
h=histogram(iSNP.ALT,'LineWidth',.1,'Normalization','pdf','LineWidth',.01);
xlim([-.5 .5]); title('ALT')

axes(ax03)
h=histogram(iSNP.DIF,'LineWidth',.1,'Normalization','pdf','LineWidth',.01);
xlim([-.5 .5]); title('DIF')

axes(ax04)
h=histogram(iSNP.NAT,'LineWidth',.1,'Normalization','pdf','LineWidth',.01);
xlim([-.5 .5]); title('NAT')




%==========================================================================
%% LOAD RANDOM TARGETS & REAL TARGETS, COMPARE ECDF
%==========================================================================
close all; clc; clearvars -except P ADSP
A=[P.home P.f 'genos_data' P.f 'SERIAL' P.f 'RUN6_INCLUDE_APOE' P.f 'DAT' P.f 'SNP.mat'];
B=[P.home P.f 'genos_data' P.f 'SERIAL' P.f 'RUN5_RANDTARGS'    P.f 'DAT' P.f 'SNP.mat'];
SNP6 = load(A,'iSNP','iLOC');
SNP5 = load(B,'iSNP','iLOC');
%%
close all; clc; clearvars -except P ADSP SNP5 SNP6
%----------------
tSNP = SNP6.iSNP;
rSNP = SNP5.iSNP;
tLOC = SNP6.iLOC;
rLOC = SNP5.iLOC;
%----------------



tDIF = mean(tSNP.DIF,1);
rDIF = mean(rSNP.DIF,1);
%distributionFitter(tDIF)
tREF = mean(tSNP.REF,1);
rREF = mean(rSNP.REF,1);
%distributionFitter(tREF)
tALT = mean(tSNP.ALT,1);
rALT = mean(rSNP.ALT,1);
%distributionFitter(tALT)



% DETERMINE EMPIRICAL CDF
[tDIFf,tDIFx] = ecdf(tDIF); tDIFcdf = [tDIFf,tDIFx];
[rDIFf,rDIFx] = ecdf(rDIF); rDIFcdf = [rDIFf,rDIFx];



% GET THE UPPER AND LOWER BOUND P=.05 CUTOFFS
rDIFcut = [max(rDIFcdf(rDIFf<.05,2))  min(rDIFcdf(rDIFf>.95,2))];



% DETERMINE WHICH SNPS ARE MORE EXTREME THAN THE CDF CUTOFF VALUES
lo = tDIF<=rDIFcut(1);
hi = tDIF>=rDIFcut(2);
sum(lo); sum(hi);
% hi = hi(1:50);
% lo = lo(1:50);



% GET cvhi & cvlo -- THE INDEX SNPS THAT PASS THRESHOLD(NOT PEOPLE) 
cvlo = tDIF(tDIF<=rDIFcut(1));
cvhi = tDIF(tDIF>=rDIFcut(2));
close all; histogram([cvlo cvhi], -.5:.01:.5)




% TAG THE SNPS IN THE LOCI TABLE THAT PASS THRESH
tLOC.sigHi = hi';
tLOC.sigLo = lo';
 
 

% GET ACTUAL NAT GENOTYPES FOR PEOPLE WITH SIGNIFICANT SNPS
tPVXhi = tSNP.PVX(:,hi) > 0;
tPVXlo = tSNP.PVX(:,lo) > 0;
 
 

% NUMBER OF SIGNIFICANT SNPS EACH PERSON HAS
nHi = sum(tPVXhi,2);   % RISK SNPS
nLo = sum(tPVXlo,2);   % PROT SNPS




% PLOT A HISTOGRAM SHOWING THE NUMBER OF PRO & RISK SNPS PER PERSON
close all
subplot(1,2,1); histogram(nLo); title('Pro SNP count per person')
subplot(1,2,2); histogram(nHi); title('Rsk SNP count per person')




% PLOT A HISTOGRAM SHOWING NUMBER SNPS X PRO|RISK CASE|CTRL
%---------------------------------------------------------------

% AVERAGE NUMBER OF RISK SNPS PER CASE
VCaRsk = nHi(tSNP.AD==1);

% AVERAGE NUMBER OF RISK SNPS PER CTRL
VCoRsk = nHi(tSNP.AD~=1);
 
% AVERAGE NUMBER OF PRO SNPS PER CASE
VCaPro = nLo(tSNP.AD==1);

% AVERAGE NUMBER OF PRO SNPS PER CTRL
VCoPro = nLo(tSNP.AD~=1);

close all
subplot(1,2,1); histogram(VCaRsk); title('DIF:Rsk allele counts per CASE|CTRL')
hold on;        histogram(VCoRsk); legend({'CASE','CTRL'})
subplot(1,2,2); histogram(VCaPro); title('DIF:Pro allele counts per CASE|CTRL')
hold on;        histogram(VCoPro);  legend({'CASE','CTRL'})





%%



Y = [PVX_muCaRsk PVX_muCoRsk PVX_muCaPro PVX_muCoPro];
X = [ 1 2 3 4];

close all
ph1 = superbar(X , Y, 'BaseValue', 0, ... 
    'BarFaceColor', parula,'BarEdgeColor', 'k'); 
 


title('This plot requires SUPERBARS file exchange'); 
hax1.XTick = [1 2 3 5 6 7]; 
hax1.XTickLabel = {'A' , 'B' , 'C' , 'D' , 'E' , 'E', 'F'}; 




%%
% GET REF|ALT|DIF VALUES FOR SIG SNPS
tPVXa = tSNP.PVX(:,hi);
tPVXb = tSNP.PVX(:,lo);





%%
tPVXa = tSNP.PVX(:,lo);
tPVXb = tSNP.PVX(:,hi);

CaPVXa = tPVXa(tPVXa.AD==1);
CaPVXb = tPVXb(tPVXb.AD==1);

CoPVXb = tPVXa(tPVXa.AD~=1);
CoPVXa = tPVXb(tPVXb.AD~=1);





close all
fh01 = figure('Units','normalized','OuterPosition',[.1 .1 .55 .7],'Color','w');
ax01 = axes('Position',[.1 .1 .8 .8],'Color','none','FontSize',16); hold on;

pn01 = probplot([rDIF(1:4000), tDIF(1:4000)]);
legend('Random SNPs','Target SNPs','Location','best')






close all;
qqplot(abs(DIF5),abs(DIF6))



% cdfplot    Empirical cumulative distribution function (eCDF) plot
% ecdf       Empirical cumulative distribution function
% ecdfhist   Histogram based on eCDF
% ksdensity  Kernel smoothing function for uni/bivariate data


[D3f,D3x,D3ci,D3cj] = ecdf(DIF6);
[D5f,D5x,D5ci,D5cj] = ecdf(DIF5);

 ecdf(DIF6)
hold on
 ecdf(DIF5)


close all
[h,stats] = cdfplot(DIF5);


close all
[h,stats] = ecdf(DIF5);

clc
[h,p,ks2stat] = kstest2(DIF6,DIF5,'Alpha',0.01)



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



%==========================================================================
%% PCA
%==========================================================================
close all; clc; clearvars -except P ADSP SNP iSNP iLOC LOOP
%----------------
iSNP = SNP.iSNP;
iLOC = SNP.iLOC;
%----------------

iSNP.AGE(iSNP.AGE<60) = 60;



% PCA
%--------------------------------
ss = statset('pca');
ss.Display = 'none';
ss.MaxIter = 100;
ss.TolFun = 1e4;
ss.TolX = 1e4;
ss.UseParallel = false;
%--------------------------------
[PCC,PCS] = pca(iSNP.DIF,'Options',ss);




close all; 
fh1 = figure('Units','normalized','Position',[.1 .1 .55 .7],'Color','w');
ax1 = axes('Position',[.1 .1 .8 .8],'Color','none','FontSize',16); hold on;

scatter(PCS(:,1),PCS(:,2),100,iSNP.AGE,'.')
ax1.FontSize = 24;

colormap((cmocean('phase'))); colorbar
colormap((AdvancedColormap('hsv'))); colorbar
colormap((AdvancedColormap('jet2'))); colorbar


%BEST COLORMAPS
%-------------------------------
%AdvancedColormap('thermal');
%AdvancedColormap('jet4');
%AdvancedColormap('bled');
%AdvancedColormap('hicontrast');
%AdvancedColormap('hsv2');
%AdvancedColormap('pastel');
%AdvancedColormap('temp');
%AdvancedColormap('fireice');
%AdvancedColormap('fireicedark');
%AdvancedColormap('rgb');
%AdvancedColormap('cmy');
%AdvancedColormap('krpoglabsvcmyw');
%cmocean('thermal')
%cmocean('phase')
%cmocean('balance')
%cmocean('tarn')
%cmocean('curl')
%cmocean('matter')
%cmocean('deep')
%cmocean('dense')
%-------------------------------


colormap((cmocean('thermo'))); colorbar
colormap(flipud(cmocean('matter'))); colorbar
AdvancedColormap('bled'); colorbar
%--------------------------------




%--------------------------------
[PCC,PCS] = pca(  iSNP.DIF' , 'Options',ss);
% close all; 
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
close all; clc; clearvars -except P ADSP SNP iSNP iLOC
%----------------
iSNP = SNP.iSNP;
iLOC = SNP.iLOC;
%----------------







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
close all; clc; clearvars -except P ADSP SNP iSNP iLOC
%----------------
iSNP = SNP.iSNP;
iLOC = SNP.iLOC;
%----------------


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
%%
close all; clc; clearvars -except P ADSP SNP
%----------------
iSNP = SNP.iSNP;
iLOC = SNP.iLOC;
%----------------

CV=mean(iSNP.REF,2);
[Y,E] = discretize(CV, (-.5:.02:.5) ); E=E';
CVmu = zeros(size(E));
for i = 1:numel(unique(E))
    CVmu(i) = mean(iSNP.AD(Y==i));
end
scatter(E,CVmu)
axis([-.5 .5 0 1])


% % make person index by CV bins
% clear CVbin CVbinN pct_true N
% for n=1:50
% CVbin{n}=logical(((CV>((n-1)*.02-.5)).*(CV<(n*.02-.5))));
% CVbinN(n)=sum(CVbin{n});
% end
% % calculate %case for each bin
% for n=1:50
%     if sum(CVbin{z}>0);N(n)=1;
%         pct_true(n)=sum((CVbin{n}).*(iSNP.AD==1))/sum(CVbin{n});
%     end
% end
% NN=[1:numel(pct_true)];
% fitlm(NN(N==1),pct_true(N==1))
% figure; plot(pct_true,'o');
% title('RUN3')





