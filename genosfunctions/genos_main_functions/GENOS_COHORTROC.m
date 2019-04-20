% GENETIC VARIANT PROCESSING PIPELINE TO PREP FOR NEURAL NET TRAINING
% DATASET: ADSP WES VCF (LATEST RELEASE)
%
% Overview of goals, issues, and processing steps:
%{
The neural network classifier requires a 2D rectangular matrix for 
training. This PERSON x LOCI matrix, will indicate for each person 
at each locus (chr:pos) whether they are homozygous reference ('0'),
heterozygous ('1'), or homozygous alternate ('2').


The pipeline will prepare to matrices, each composed of a unique group
of people. This will provide one dataset for neural net training, and
another for testing the accuracy of the classifier.


Instead of listing every processing detail here, those will be explained
in each block of code along the way. Though, I will explain a data
stratification issue we have previously identified so it's clear where
things are going out-of-the-gate. 


ADSP has combined sequencing data from 25 different consortium studies
(wherever you see the code mention 'STUDY' or 'COHORT', assume they are 
all synonymously referring to these 25 data-sources). Because of this, 
major stratification issues lurk in this dataset. Further, many of these 
PHEN provide very unequal numbers of Alzheimer's patients (CASE) and 
controls (CTRL).


Here is a summary table of those 25 STUDY-PHEN, and their 
share of participant contributions...

%% COID CONSORT	STUDY_NAME  STUDY   FAM PEOPLE  ENR     CASE	CTRL	TOT    %CA %CO  EQL	CULL	NOTES
01  DGC     Adult Chng  ACT     0	1277	0       323     945     1268	25	75	323	-622	KEEP
02  ADGC	AD Centers  ADC     0	3263	0       2438	817     3255	75	25	817	1621	KEEP
03  CHARGE	Athrosclro  ARIC	0	57      0       39      18      57      68	32	18	21      REMOVE
04  CHARGE	Aus Stroke	ASPS	0	126     0       121     5       126     96	4	5	116     REMOVE
05  ADGC	Chic Aging  CHAP	0	231     0       27      204     231     12	88	27	-177	REMOVE
06  CHARGE	CardioHlth  CHS     0	840     0       250     583     833     30	70	250	-333	KEEP
07  ADGC	Hispanic S  CUHS	361	0       331     160     171     331     48	52	160	-11     HISPANIC
08  CHARGE	Erasmus Er  ERF     30	45      0       45      0       45      100	0	0	45      REMOVE
09  CHARGE	Framingham  FHS     0	581     0       157     424     581     27	73	157	-267	KEEP
10	ADGC	Gene Diffs  GDF     0	207     0       111     96      207     54	46	96	15      KEEP
11	ADGC	NIA - LOAD  LOAD	80  113     363     367     109     476     77	23	109	258     KEEP
12	ADGC	Aging Proj  MAP     0	415     0       138     277     415     33	67	138	-139	KEEP
13	ADGC	Mayo Clini  MAYO	0	349     0       250     99      349     72	28	99	151     KEEP
14	ADGC	Miami Univ  MIA     61	181     19      186     14      200     93	7	14	172     REMOVE
15	ADGC	AD Genetic  MIR     0	284     47      316     15      331     95	5	15	301     REMOVE
16	ADGC	Mayo cl PD  MPD     0	20      0       0       20      20      0	100	0	-20     REMOVE
17	ADGC	NationC AD  NCRD	18	108     52      160     0       160     100	0	0	160     REMOVE
18	ADGC	Wash Unive  RAS     0	0       46      46      0       46      100	0	0	46      REMOVE
19	ADGC	Relig Ordr  ROS     0	351     0       154     197     351     44	56	154	-43     KEEP
20	CHARGE	RotterdamS  RS      0	1089	0       276     813     1089	25	75	276	-537	KEEP
21	ADGC	Texas AD S  TARC	0	144     0       132     12      144     92	8	12	120     REMOVE
22	ADGC	Un Toronto  TOR     0	0       9       9       0       9       100	0	0	9       REMOVE
23	ADGC	Vanderbilt  VAN     6	235     1       210     26      236     89	11	26	184     REMOVE
24	ADGC	WashNY Age  WCAP	0	147     3       34      116     150     23	77	34	-82     REMOVE
25	ADGC 	Univ Penns  UPN     40 	0       3       0       0       0       0	0	0	0       REMOVE



Along the way, the code will take steps to mitigate stratification
issues due to COHORT-specific effects. 

...I'm trying to think if there is anything else you should know...
I can't, sooo, I guess we should just go-ahead and start with step #1.
(aka the first chance for something to go wrong)...

%}

%% STEP-1: LOAD THE DATASET

close all; clear; clc; rng('shuffle');
genosdir = fileparts(which('GENOS.m'));
cd(genosdir);


subfuncpath = [genosdir '/genosfunctions'];
datasetpath = [genosdir '/genosdata'];
gpath = [genosdir ':' subfuncpath ':' datasetpath];
addpath(gpath)


% which('ADSP.mat')
% load('ADSP.mat')
which('ADSP_MINI.mat')
load('ADSP_MINI.mat')


clc; rng('shuffle'); 
disp('dataset loaded')
clearvars -except ADSP



%% CARBON COPY MAIN VARIABLES FROM ADSP.STRUCT


LOCI = ADSP.LOCI;
CASE = ADSP.CASE;
CTRL = ADSP.CTRL;
PHEN = ADSP.PHEN;
USNP = ADSP.USNP;


clearvars -except ADSP LOCI CASE CTRL PHEN USNP
clc; disp(LOCI(1:9,:))




%% IDENTIFY KEEPER PHEN AND MAKE CASE:CTRL EQUAL SIZE


PHECASE = PHEN(PHEN.AD==1,:);
PHECTRL = PHEN(PHEN.AD==0,:);

nCA = 1:size(PHECASE,1);
nCO = 1:size(PHECTRL,1);

iCA = randperm( size(PHECASE,1) , round(size(PHECASE,1)*0.7) );
iCO = randperm( size(PHECTRL,1) , round(size(PHECTRL,1)*0.7) );


[aCA,bCA] = ismember(nCA,iCA);
[aCO,bCO] = ismember(nCO,iCO);



TRCASE = PHECASE(aCA,:);
TRCTRL = PHECTRL(aCO,:);

TECASE = PHECASE(~aCA,:);
TECTRL = PHECTRL(~aCO,:);





clc; close all;
disp(TRCASE(1:9,:))
subplot(2,2,1), histogram(TRCASE.TOTvars); title('TRAIN CASE')
subplot(2,2,2), histogram(TRCTRL.TOTvars); title('TRAIN CTRL')
subplot(2,2,3), histogram(TECASE.TOTvars); title('TEST  CASE')
subplot(2,2,4), histogram(TECTRL.TOTvars); title('TEST  CTRL')
disp('------------------------------------')
disp('TOTAL VARIANTS PER-PERSON')
disp('                MIN    MAX')
fprintf('TRAINING CASE: %.0f  %.0f \n',min(TRCASE.TOTvars),  max(TRCASE.TOTvars))
fprintf('TRAINING CTRL: %.0f  %.0f \n',min(TRCTRL.TOTvars),  max(TRCTRL.TOTvars))
fprintf('TESTING  CASE: %.0f  %.0f \n',min(TECASE.TOTvars),  max(TECASE.TOTvars))
fprintf('TESTING  CTRL: %.0f  %.0f \n',min(TECTRL.TOTvars),  max(TECTRL.TOTvars))
disp('------------------------------------')
pause(1)
clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
TRCASE TRCTRL TECASE TECTRL











%###############################################################
%%  USE PER-PERSON VARIANT COUNTS TO IDENTIFY OUTLIERS
%###############################################################

% % MAKE SURE EVERYONE IN CTRL GROUP IS OVER 65 YEARS OLD
% PHETRCTRL(PHETRCTRL.AGE < 65 , :) = [];
% PHETECTRL(PHETECTRL.AGE < 65 , :) = [];
% 
% 
% Qlo = quantile(PHETRCASE.TOTvars,.01);
% Qhi = quantile(PHETRCASE.TOTvars,.99);
% TRCASE  = PHETRCASE(((PHETRCASE.TOTvars > Qlo) & (PHETRCASE.TOTvars < Qhi)),:);
% 
% 
% Qlo = quantile(PHETRCTRL.TOTvars,.01);
% Qhi = quantile(PHETRCTRL.TOTvars,.99);
% TRCTRL  = PHETRCTRL(((PHETRCTRL.TOTvars > Qlo) & (PHETRCTRL.TOTvars < Qhi)),:);
% 
% 
% Qlo = quantile(PHETECASE.TOTvars,.01);
% Qhi = quantile(PHETECASE.TOTvars,.99);
% TECASE  = PHETECASE(((PHETECASE.TOTvars > Qlo) & (PHETECASE.TOTvars < Qhi)),:);
% 
% 
% Qlo = quantile(PHETECTRL.TOTvars,.01);
% Qhi = quantile(PHETECTRL.TOTvars,.99);
% TECTRL  = PHETECTRL(((PHETECTRL.TOTvars > Qlo) & (PHETECTRL.TOTvars < Qhi)),:);
% 


% clc; close all;
% subplot(2,2,1), histogram(TRCASE.TOTvars); title('TRAIN CASE')
% subplot(2,2,2), histogram(TRCTRL.TOTvars); title('TRAIN CTRL')
% subplot(2,2,3), histogram(TECASE.TOTvars); title('TEST  CASE')
% subplot(2,2,4), histogram(TECTRL.TOTvars); title('TEST  CTRL')
% disp('------------------------------------')
% disp('TOTAL VARIANTS PER-PERSON')
% disp('                MIN    MAX')
% fprintf('TRAINING CASE: %.0f  %.0f \n',min(TRCASE.TOTvars),  max(TRCASE.TOTvars))
% fprintf('TRAINING CTRL: %.0f  %.0f \n',min(TRCTRL.TOTvars),  max(TRCTRL.TOTvars))
% fprintf('TESTING  CASE: %.0f  %.0f \n',min(TECASE.TOTvars),  max(TECASE.TOTvars))
% fprintf('TESTING  CTRL: %.0f  %.0f \n',min(TECTRL.TOTvars),  max(TECTRL.TOTvars))
% disp('------------------------------------')
% pause(1)
% clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
% TRCASE TRCTRL TECASE TECTRL








%% DISCARD VARIANT LOCI WITH KNOWN ULTRA-RARE ALT ALLELE COUNTS

PASS = (LOCI.CASEALT > 10) | (LOCI.CTRLALT > 10);

LOCI  = LOCI(PASS,:);
CASE  = CASE(PASS);
CTRL  = CTRL(PASS);
USNP  = USNP(PASS);





%###############################################################
%%          COUNT NUMBER OF VARIANTS PER LOCI
%###############################################################


[TRCASEN,TRCTRLN,TECASEN,TECTRLN,TRCASEUN,TRCTRLUN,TECASEUN,TECTRLUN] = ...
    snpsum(CASE,CTRL,USNP,TRCASE.SRR,TRCTRL.SRR,TECASE.SRR,TECTRL.SRR);


% SAVE COUNTS AS NEW TABLE COLUMNS
LOCI.TRCASEREF = (numel(TRCASE.SRR)*2) - (TRCASEUN.*2) - TRCASEN;
LOCI.TRCTRLREF = (numel(TRCTRL.SRR)*2) - (TRCTRLUN.*2) - TRCTRLN;
LOCI.TRCASEALT = TRCASEN;
LOCI.TRCTRLALT = TRCTRLN;


LOCI.TECASEREF = (numel(TECASE.SRR)*2) - (TECASEUN.*2) - TECASEN;
LOCI.TECTRLREF = (numel(TECTRL.SRR)*2) - (TECTRLUN.*2) - TECTRLN;
LOCI.TECASEALT = TECASEN;
LOCI.TECTRLALT = TECTRLN;




LOCI.CASEREF = LOCI.TRCASEREF;
LOCI.CTRLREF = LOCI.TRCTRLREF;
LOCI.CASEALT = LOCI.TRCASEALT;
LOCI.CTRLALT = LOCI.TRCTRLALT;






close all
fh01 = figure('Units','normalized','OuterPosition',[.01 .25 .48 .62],'Color','w');
ax01 = axes('Position',[.06 .06 .9 .9],'Color','none');
ph01 = histogram(TRCASEN(TRCASEN>5), 40,'DisplayStyle','stairs',...
    'LineWidth',3,'Normalization','count'); hold on;
ph02 = histogram(TRCTRLN(TRCTRLN>5), 40,'DisplayStyle','stairs',...
    'LineWidth',3,'Normalization','count'); hold on;
ph03 = histogram(TECASEN(TECASEN>5), 40,'DisplayStyle','stairs',...
    'LineWidth',3,'Normalization','count'); hold on;
ph04 = histogram(TECTRLN(TECTRLN>5), 40,'DisplayStyle','stairs',...
    'LineWidth',3,'Normalization','count'); hold on;
title('DISTRIBUTION OF ALT ALLELES')

fh21 = figure('Units','normalized','OuterPosition',[.5 .25 .48 .62],'Color','w');
ax21 = axes('Position',[.06 .06 .9 .9],'Color','none');
ph21 = histogram(LOCI.TRCASEREF(LOCI.TRCASEREF<max(LOCI.TRCASEREF)-10),...
    40,'DisplayStyle','stairs','LineWidth',3,'Normalization','count'); hold on;
ph22 = histogram(LOCI.TRCTRLREF(LOCI.TRCTRLREF<max(LOCI.TRCTRLREF)-10),...
    40,'DisplayStyle','stairs','LineWidth',3,'Normalization','count'); hold on;
ph23 = histogram(LOCI.TECASEREF(LOCI.TECASEREF<max(LOCI.TECASEREF)-10),...
    40,'DisplayStyle','stairs','LineWidth',3,'Normalization','count'); hold on;
ph24 = histogram(LOCI.TECTRLREF(LOCI.TECTRLREF<max(LOCI.TECTRLREF)-10),...
    40,'DisplayStyle','stairs','LineWidth',3,'Normalization','count'); hold on;
title('DISTRIBUTION OF REF ALLELES')


pause(2); close all;

clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
TRCASE TRCTRL TECASE TECTRL






%###############################################################
%%               COMPUTE FISHER'S P-VALUE
%###############################################################


% COMPUTE FISHERS STATISTICS FOR THE TRAINING GROUP
[FISHP, FISHOR] = fishp_mex(LOCI.TRCASEREF,LOCI.TRCASEALT,...
                            LOCI.TRCTRLREF,LOCI.TRCTRLALT);

LOCI.TRFISHP  = FISHP;
LOCI.TRFISHOR = FISHOR;




% COMPUTE FISHERS STATISTICS FOR THE TRAINING GROUP
[FISHP, FISHOR] = fishp_mex(LOCI.TECASEREF, LOCI.TECASEALT,...
                            LOCI.TECTRLREF, LOCI.TECTRLALT);

LOCI.TEFISHP  = FISHP;
LOCI.TEFISHOR = FISHOR;




close all
x = -log(LOCI.TRFISHP);
y = -log(LOCI.TEFISHP);
histogram(x(x>6.9)); hold on
histogram(y(y>6.9));
pause(3)


clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
TRCASE TRCTRL TECASE TECTRL




%% COMPUTE CHI SQUARE VALUE
%{
% % COMPUTE CHI SQUARED STATISTICS FOR THE TRAINING GROUP
% [CHIP, CHIOR] = chisq(LOCI.TRCASEREF,LOCI.TRCASEALT,...
%                       LOCI.TRCTRLREF,LOCI.TRCTRLALT);
% 
% LOCI.TRCHIP  = CHIP;
% LOCI.TRCHIOR = CHIOR;
% 
% 
% 
% % COMPUTE CHI SQUARED STATISTICS FOR THE TRAINING GROUP
% [CHIP, CHIOR] = chisq(LOCI.TECASEREF, LOCI.TECASEALT,...
%                       LOCI.TECTRLREF, LOCI.TECTRLALT);
% 
% LOCI.TECHIP  = CHIP;
% LOCI.TECHIOR = CHIOR;
%
% clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
% TRCASE TRCTRL TECASE TECTRL
%}










%###############################################################
%###############################################################
%###############################################################
%%
%
% PREPARE DATA FOR NEURAL NET CLASSIFIER SUPERVISED LEARNING
%
%###############################################################
%###############################################################
%###############################################################
clc;




%% STORE VARIABLES FOR NEURAL NETWORK TRAINING AS 'VLOCI'

VLOCI     = LOCI;
VCASE     = CASE;
VCTRL     = CTRL;
VUSNP     = USNP;
VTRCASE   = TRCASE;
VTRCTRL   = TRCTRL;
VTECASE   = TECASE;
VTECTRL   = TECTRL;



clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL






%% REMOVE KNOWN AD ALZ GENES BY NAME
%{

% ALZGENES = string({"APOE";"BIN1";"CLU";"ABCA7";"CR1";...
%                    "PICALM";"MS4A6A";"CD33";"MS4A4E";"CD2AP"});
% 
% for nn = 1:numel(ALZGENES)
% 
%     x = strcmp(VLOCI.GENE,ALZGENES{nn});
%     sum(x)
% 
%     VLOCI(x,:) = [];
%     VCASE(x) = [];
%     VCTRL(x) = [];
% 
% end
% VLOCI.VID  = (1:size(VLOCI,1))';



% ALZGENES = string({"APOE"});
% 
% for nn = 1:numel(ALZGENES)
% 
%     x = strcmp(VLOCI.GENE,ALZGENES{nn});
%     sum(x)
% 
%     VLOCI(x,:) = [];
%     VCASE(x) = [];
%     VCTRL(x) = [];
% 
% end
% VLOCI.VID  = (1:size(VLOCI,1))';
% 
% clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
% VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL
%}






%% TAKE THE TOP N GENES FOR NEURAL NET CLASSIFIER TRAINING



% FIRST SORT BY TRAINING GROUP FISHP
[~,i]  = sort(VLOCI.TRFISHP);
VLOCI  = VLOCI(i,:);
VCASE  = VCASE(i);
VCTRL  = VCTRL(i);
VUSNP  = VUSNP(i);





% EXTRACT TOP-N NUMBER OF VARIANTS
N = 800;

VLOCI  = VLOCI(1:N,:);
VCASE  = VCASE(1:N);
VCTRL  = VCTRL(1:N);
VUSNP  = VUSNP(1:N);






% SET MAIN FISHP TO TRAINING GROUP FISHP
VLOCI.FISHP      = VLOCI.TRFISHP;
VLOCI.FISHOR     = VLOCI.TRFISHOR;
VLOCI.CASEREF    = VLOCI.TRCASEREF;
VLOCI.CASEALT    = VLOCI.TRCASEALT;
VLOCI.CTRLREF    = VLOCI.TRCTRLREF;
VLOCI.CTRLALT    = VLOCI.TRCTRLALT;





clc; disp(VLOCI(1:9,:))
fprintf('\n LOCI COUNT FOR FEATURE MATRIX: %.0f \n\n',size(VLOCI,1))

clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL





%% VISUALIZE CORRELATION BETWEEN TRAINING AND TESTING FISHER'S P-VALUES


clc; close all;
fh1 = figure('Units','normalized','OuterPosition',[.1 .05 .7 .92],...
             'Color','w','MenuBar','none');
ax1 = axes('Position',[.06 .06 .9 .9],'Color','none');


ph1 = loglog(-log(VLOCI.TRFISHP(1:end))+1 , -log(VLOCI.TEFISHP(1:end))+1,'-mo',...
    'LineStyle','none',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.49 1 .63],...
    'MarkerSize',10);
grid on; axis tight



[RHO,PVAL] = corr(-log(VLOCI.TRFISHP),-log(VLOCI.TEFISHP));





clc; disp(VLOCI(1:50,:))
fprintf('\n CORRELATION: %.2f (P ~ %.1g) \n\n',RHO,PVAL)
pause(2); close all; disp(' ')

clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL





%###################################################################
%%      MAKE  RECTANGLE  NEURAL NET  VARIANT MATRIX
%###################################################################



TRPHE = [VTRCASE; VTRCTRL];
TEPHE = [VTECASE; VTECTRL];



% SCRAMBLE TRAINING PHENOTYPE ORDER
N      = size(TRPHE,1);  % Total number of people
k      = randperm(N)';   % Get N random ints in range 1:N
TRPHE  = TRPHE(k,:);     % Scramble Phenotype table

% SCRAMBLE TESTING PHENOTYPE ORDER
N      = size(TEPHE,1);  % Total number of people
k      = randperm(N)';   % Get N random ints in range 1:N
TEPHE  = TEPHE(k,:);     % Scramble Phenotype table




% MAKE THE NEURAL NET TRAINING & TESTING MATRICES

[FULLTRX, TRX, TRL] = makennet(VLOCI,VCASE,VCTRL,VUSNP,TRPHE,0);

[FULLTEX, TEX, TEL] = makennet(VLOCI,VCASE,VCTRL,VUSNP,TEPHE,0);



clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE...
FULLTRX FULLTEX TRX TRL TEX TEL




%% EXPORT TRAINING AND TESTING MATRICES FOR USE IN TENSORFLOW
%
% writetable(VLOCI,'VARIANTLOCI.csv')
% 
% writetable(TRPHE,'TRAINPHENOTYPES.csv')
%
% writetable(TEPHE,'TESTPHENOTYPES.csv')
% 
% writetable(array2table(TRX),'TRAINMATRIX.csv')
% 
% writetable(array2table(TRL),'TRAINLABELS.csv')
% 
% writetable(array2table(TEX),'TESTMATRIX.csv')
% 
% writetable(array2table(TEL),'TESTLABELS.csv')


return
%##########################################################################
%
%%          MACHINE LEARNING USING ARTIFICIAL NEURAL NETWORKS
%
%##########################################################################


% TRL_COH = TRPHE.COHORTNUM;
% TEL_COH = TEPHE.COHORTNUM;
TRL_COH = FULLTRX(2:end,4);
TEL_COH = FULLTEX(2:end,4);

TRC = labelmxr(TRL_COH)';
TEC = labelmxr(TEL_COH)';



NN = patternnet([300 200 100]);

net = train(NN,TRX,TRC);

[NNPerf] = netstats(net,TRX,TRC,TEX,TEC,.05,.85,.15,1);






% NN = patternnet([300 200 100]);
% 
% net = train(NN,TRX,TRL);
% 
% [NNPerf] = netstats(net,TRX,TRL,TEX,TEL,.05,.85,.15,1);




clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE...
FULLTRX FULLTEX TRX TRL TEX TEL net netz




%##########################################################################
%%        ASSESS ABILITY OF NEURAL NET TO GUESS PATIENT COHORT
%##########################################################################



ACTIVITY = net(TEX);
CLASSES    = vec2ind(ACTIVITY);
ACTIVATION = ACTIVITY(1,:);
GUESSES    = round(ACTIVITY);

iCORRECT = TEL(1,:) == GUESSES;

iHICASE  = ACTIVATION>.8;
iHICTRL  = ACTIVATION<.2;
iHI      = (ACTIVATION<.20) | (ACTIVATION>.80);

mean(iCORRECT(iHI))




























%% RUN THE TRAINING A BUNCH OF TIMES (DUAL HIDDEN LAYER)


nLoops = 10;

netz={};
NETi = zeros(nLoops,7);
for nn = 1:nLoops
    netz{nn} = patternnet([300 100 100]);
    netz{nn} = train(netz{nn},TRX,TRL);

    [NNPerf] = netstats(netz{nn},TRX,TRL,TEX,TEL,.05,.85,.15,0);
    NETi(nn,:) = [nn NNPerf];


disp(' ')
disp(['        NETn    All TRAIN     ALL TEST        '...
      'MID CORRECT    MID POP   HIGH CORRECT   HIGH POP'])
disp(NETi)
disp(' ')
end



clc
disp('Dual-hidden layer L1:300n L2:100n L3:100')
disp(' ')
disp(['        NETn    All TRAIN     ALL TEST        '...
      'MID CORRECT    MID POP   HIGH CORRECT   HIGH POP'])
disp(NETi)
disp(' ')
disp('Mean:'); disp(mean(NETi))
disp('ALL VALUES ABOVE ARE PERCENTAGES')




clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE TR TE...
TRX TRL TEX TEL net netz NETi











%##########################################################################
%
%%     ASSESS NEURAL NET PERFORMANCE BASED ON CLASSIFIER CONFIDENCE
%
%##########################################################################




%% SET THE CURRENT NEURAL NET TO THE BEST TRAINED NET FROM ABOVE

NS = NETi ./ size(netz,2);
BMOscore = NS(:,6) .* (NS(:,7) ./3) + NS(:,4) + NS(:,3) + NS(:,6);
[~,topBMO] = max(BMOscore);
[~,i] = sort(BMOscore,'descend');

netz = netz(i);
NETi = NETi(i,:);

disp(['  Rank TrAll  ACOR  MCOR  MPOP  HCOR  HPOP'])
disp(round(NETi,0))

net = netz{1};





clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE TR TE...
TRX TRL TEX TEL net netz NETi















%##########################################################################
%
%%       ASSESS NEURAL NET ABILITY TO DETERMINE COHORT ORIGIN
%
%##########################################################################



ACTIVATION = net(TEX);
CLASSES    = vec2ind(ACTIVATION);
ACTIVATION = ACTIVATION(1,:);
GUESSES    = round(ACTIVATION);

iCORRECT = TEL(1,:) == GUESSES;

iHICASE  = ACTIVATION>.8;
iHICTRL  = ACTIVATION<.2;
iHI      = (ACTIVATION<.20) | (ACTIVATION>.80);

mean(iCORRECT(iHI))


















































%% GET INDEX OF HIGH CONFIDENCE CORRECT DECISIONS ON **TEST** SET

ACTIVATION = net(TEX);
CLASSES    = vec2ind(ACTIVATION);
ACTIVATION = ACTIVATION(1,:);
GUESSES    = round(ACTIVATION);

iCORRECT = TEL(1,:) == GUESSES;

iHICASE  = ACTIVATION>.8;
iHICTRL  = ACTIVATION<.2;
iHI      = (ACTIVATION<.20) | (ACTIVATION>.80);

mean(iCORRECT(iHI))



clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE TR TE...
TRX TRL TEX TEL net netz NETi






%##########################################################################
%%       MAKE UNEQUAL SIZE TEST SETS WITH X% CASES Y% CTRL
%##########################################################################
clc

TESL = TEL;
TESX = TEX;

[~,i,~] = find(TESL(1,:));  % index of cases in label & test matrix
[~,j,~] = find(TESL(2,:));  % index of ctrls in label & test matrix

nTOTA = size(TESL,2);       % number of people in test set
nCASE = numel(i);           % number of cases in test set
nCTRL = numel(j);           % number of ctrls in test set

MAXN = min([nCASE nCTRL]);
n = round(MAXN / .7);


% BUILD SET WHERE CASE PCT > CTRL PCT

GetNca = MAXN;
GetNco = n - MAXN;

CAi = i(1:GetNca);
COi = j(1:GetNco);

TEX_70Ca_30Co = TESX(:,[CAi COi]);
TEL_70Ca_30Co = TESL(:,[CAi COi]);


fprintf('\nN case: %.0f (P %0.2g)\n',GetNca,round(GetNca/(GetNca+GetNco),2))
fprintf('\nN ctrl: %.0f (P %0.2g)\n',GetNco,round(GetNco/(GetNca+GetNco),2))



% BUILD SET WHERE CTRL PCT > CASE PCT

GetNca = n - MAXN;
GetNco = MAXN;

CAi = i(1:GetNca);
COi = j(1:GetNco);

TEX_30Ca_70Co = TESX(:,[CAi COi]);
TEL_30Ca_70Co = TESL(:,[CAi COi]);


fprintf('\nN case: %.0f (P %.2g)\n',GetNca,round(GetNca/(GetNca+GetNco),2))
fprintf('\nN ctrl: %.0f (P %.2g)\n',GetNco,round(GetNco/(GetNca+GetNco),2))



clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE TR TE...
TRX TRL TEX TEL net netz NETi ...
TEX_70Ca_30Co TEL_70Ca_30Co TEX_30Ca_70Co TEL_30Ca_70Co





%% EVALUATE PERFORMANCE OF UNEQUAL SIZE SETS CREATED ABOVE


%---- RUN 70:30 CASE:CTRL MATRIX THROUGH CLASSIFIER
NNOUT     = net(TEX_70Ca_30Co);
ACTIVITY  = NNOUT(1,:);
CLASSES   = vec2ind(NNOUT);
GUESSES   = round(ACTIVITY);

iHI = (ACTIVITY<.20) | (ACTIVITY>.80);

iCORRECT = TEL_70Ca_30Co(1,:) == GUESSES;

PCOR_70Ca30Co = mean(iCORRECT(iHI));



%---- RUN 30:70 CASE:CTRL MATRIX THROUGH CLASSIFIER
NNOUT     = net(TEX_30Ca_70Co);
ACTIVITY  = NNOUT(1,:);
CLASSES   = vec2ind(NNOUT);
GUESSES   = round(ACTIVITY);

iHI = (ACTIVITY<.20) | (ACTIVITY>.80);

iCORRECT = TEL_30Ca_70Co(1,:) == GUESSES;

PCOR_30Ca70Co = mean(iCORRECT(iHI));



%---- REPORT PERFORMANCE OF 70:30 TEST SETS
clc;
fprintf('\nPCT correct when 70:30 case:ctrl %.3f \n',PCOR_70Ca30Co)
fprintf('\nPCT correct when 30:70 case:ctrl %.3f \n',PCOR_30Ca70Co)



clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE TR TE...
TRX TRL TEX TEL net netz NETi ...
TEX_70Ca_30Co TEL_70Ca_30Co TEX_30Ca_70Co TEL_30Ca_70Co










%################################################################
%%  ASSESS VARIANT PROFILE OF HIGH-CONFIDENCE POPULATION
%################################################################


NN = patternnet([200 100 100]);

NN = configure(NN,TRAINMX,TRAINLAB);

net = train(NN,TRAINMX,TRAINLAB);

% perf = perform(net,TESTLAB,ACTIVATION);

HC = [.80 .20]

[NNPerf] = netstats(net,TRAINMX,TRAINLAB,TESTMX,TESTLAB,.05,HC(1),HC(2),1);









%% EVAL USING TRAINING DATA

ACTIVATION = net(TRAINMX);
GUESSES    = round(ACTIVATION);
CLASSES    = vec2ind(ACTIVATION);


CORRECT = sum(TRAINLAB == GUESSES)>0;

HICONF  = sum(ACTIVATION>HC(1) | ACTIVATION<HC(2))>0;

HICONFCASE  = ACTIVATION(1,:)>HC(1);
HICONFCTRL  = ACTIVATION(2,:)>HC(1);

LOCONF  = sum((ACTIVATION<.6)&(ACTIVATION>.4))>0;


HICASE = TRAINMX(:,HICONFCASE);
HICTRL = TRAINMX(:,HICONFCTRL);
HICASE = HICASE(:,1:100);
HICTRL = HICTRL(:,1:100);
% HICASE = mean(HICASE,2);
% HICTRL = mean(HICTRL,2);



close all
fh1 = figure('Units','normalized','OuterPosition',[.01 .05 .95 .92],...
             'Color','w','MenuBar','none');
ax1 = axes('Position',[.06 .05 .43 .90],'Color','none');
ax2 = axes('Position',[.55 .05 .43 .90],'Color','none');


axes(ax1)
imagesc(HICASE)
title('\fontsize{16} High Confidence Cases')
ax1.YTick = 1:size(HICASE,1);
YLabs = cellstr([repmat('\fontsize{8} ',length(VLOCI.GENE),1), char(VLOCI.GENE)]);
ax1.YTickLabel = cellstr(YLabs);
ax1.FontName = 'Courier';



axes(ax2)
imagesc(HICTRL)
title('\fontsize{16} High Confidence Controls')
ax2.YTick = 1:size(HICTRL,1);
ax2.YTickLabel = cellstr(YLabs);
ax2.FontName = 'Courier';


map =  [.1  .1  .6 ;
        .2  .2  .8 ;
        .95  .4  .2 ;
        1   1   .6];
colormap(map)

% C = jet();
% colormap(C(1:end-6,:))
% colorbar





%% AVERAGE VARIANT PROFILE FOR HIGH CONFIDENCE PREDICTIONS

HICASE = TRAINMX(:,HICONFCASE);
HICTRL = TRAINMX(:,HICONFCTRL);

HICASE = mean(HICASE,2);
HICTRL = mean(HICTRL,2);



close all
fh1 = figure('Units','normalized','OuterPosition',[.01 .05 .95 .92],...
             'Color','w','MenuBar','none');
ax1 = axes('Position',[.06 .05 .25 .90],'Color','none');
ax2 = axes('Position',[.33 .05 .27 .90],'Color','none');
ax3 = axes('Position',[.66 .05 .25 .90],'Color','none');


axes(ax1)
imagesc(HICASE)
title('\fontsize{16} Average High Confidence Case')
ax1.YTick = 1:size(HICASE,1);
YLabs = cellstr([repmat('\fontsize{9} ',length(VLOCI.GENE),1), char(VLOCI.GENE)]);
ax1.YTickLabel = cellstr(YLabs);
ax1.XTickLabel = [];
ax1.FontName = 'Courier';

axes(ax2)
imagesc(HICTRL)
title('\fontsize{16} Average High Confidence Control')
ax2.YTick = 1:size(HICTRL,1);
ax2.YTickLabel = [];
ax2.XTickLabel = [];
ax2.FontName = 'Courier';

C = jet();
colormap(C(1:end-6,:))
colorbar


axes(ax3)
imagesc( abs(HICASE-HICTRL))
title('\fontsize{16} abs |CASE mean - CTRL mean|')
ax3.YTick = 1:size(HICTRL,1);
ax3.YTickLabel = cellstr(YLabs);
ax3.XTickLabel = [];
ax3.FontName = 'Courier';


C = jet();
colormap(ax3, C(12:end-4,:))
colorbar



%%
clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL...
TRNN TENN TRAINMX TRAINLAB TESTMX TESTLAB TRAINX TESTX...
NN net NNPerf HC







%% EVAL USING HOLD-OUT TEST DATA

ACTIVATION = net(TESTMX);
GUESSES    = round(ACTIVATION);
CLASSES    = vec2ind(ACTIVATION);


CORRECT = sum(TESTLAB == GUESSES)>0;

HICONF  = sum(ACTIVATION>HC(1) | ACTIVATION<HC(2))>0;

HICONFCASE  = ACTIVATION(1,:)>HC(1);
HICONFCTRL  = ACTIVATION(2,:)>HC(1);

LOCONF  = sum((ACTIVATION<.6)&(ACTIVATION>.4))>0;


HICASE = TESTMX(:,HICONFCASE);
HICTRL = TESTMX(:,HICONFCTRL);
HICASE = HICASE(:,1:50);
HICTRL = HICTRL(:,1:50);



close all
fh1 = figure('Units','normalized','OuterPosition',[.01 .06 .9 .7],...
             'Color','w','MenuBar','none');
ax1 = axes('Position',[.05 .09 .42 .83],'Color','none');
ax2 = axes('Position',[.55 .09 .42 .83],'Color','none');


axes(ax1)
imagesc(HICASE)
title('High Confidence Cases')

axes(ax2)
imagesc(HICTRL)
title('High Confidence Controls')


map =  [.2  .2  .7
        .1  .1  .6
        .9  .4  .4
        1   1   .6];
colormap(map)
colorbar








clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL...
TRNN TENN TRAINMX TRAINLAB TESTMX TESTLAB TRAINX TESTX...
NN net NNPerf HC











return
%################################################################
%%  PERFORM MACHINE LEARNING USING HOMEBREW METHODS
%################################################################
% LOGISTIC REGRESSION NEURAL NET WITH GRADIENT DECENT

Niters = 20;
RES = zeros(Niters,5);

% szTx = size(TRAINX,1);
% nP = fliplr(round(linspace(100,szTx,Niters)))';



for nn = 1:Niters


clearvars -except ADSP GENB LOCI CASE CTRL PHEN RES nn...
TRCASE TRCTRL TECASE TECTRL TRCASEN TRCTRLN TECASEN TECTRLN...
VLOCI VCASE VCTRL VTRCASE VTRCTRL VTECASE VTECTRL...
TRAINMX TRAINLAB TESTMX TESTLAB TRNN TENN nP


% The *TRNN* matrix is formatted as follows
%
% ROW-1: VLOCI.VID of variant for TRNN(1,4:end)
% COL-1: PHEN.SRR participant ID
% COL-2: PHEN.AD participant AD CASE\CTRL status
% COL-3: bias input node 'ones'
%
% TRNN(2:end,3:end) ... the bias + the variant data
% TRNN(2:end,4:end) ... just the variant data



% DEAL OUT VARIANTS AND LABELS FOR TRAINING DATA
ADNN    =  TRNN;
TRAINX  =  ADNN(2:end,4:end);  
TRAINL  =  ADNN(2:end,2) + 1;


% DEAL OUT VARIANTS AND LABELS FOR TESTING DATA
TESTX   =  TENN(2:end,4:end);
TESTL   =  TENN(2:end,2) + 1;



% RANDOMIZE ROWS
randp        = randperm(size(TRAINX,1));
TRAINX       = TRAINX(randp,:);
TRAINL       = TRAINL(randp,:);


% % UNCOMMENT FOR SERIAL REDUCTION
% TRAINX = TRAINX(  1:nP(nn) ,:);
% TRAINL = TRAINL(  1:nP(nn) ,:);



% SET NUMBER OF L1 NEURONS EQUAL TO NUMBER OF VARIANTS
L1neurons  = size(TRAINX,2);
nLabels = 2;



% NEURAL NET PARAMETERS
%----------------------------------------------------------------------
lambda = 0.005;                 % .001 - .01 is a good operational window
epsInit = 0.22;                 % random initial theta weights
maxIters  = 50;                 % 20-50 iterations should be sufficient
L2neurons = 35;                 % 10 - 50 neurons should be sufficient
%L2neurons = 20 + randi(20);    % 10 - 50 neurons should be sufficient
%maxIters  = 50 + randi(20);    % 20-50 iterations should be sufficient
%----------------------------------------------------------------------



% ESTABLISH NEURON ACTIVATION FUNCTION
sig    = @(z) 1 ./ (1 + exp(-z) );     % sigmoid activation function
sigrad = @(g) sig(g) .* (1-sig(g));    % sigmoid gradient function
sstf = @(n) 2 ./ (1 + exp(-2*n)) - 1;  % sigmoid symetric transfer func
% subplot(2,1,1), plot(sstf(-5:.1:5))
% subplot(2,1,2), plot(sig(-5:.1:5))



clearvars -except ADSP GENB LOCI CASE CTRL PHEN RES nn...
TRCASE TRCTRL TECASE TECTRL TRCASEN TRCTRLN TECASEN TECTRLN...
VLOCI VCASE VCTRL VTRCASE VTRCTRL VTECASE VTECTRL...
TRAINMX TRAINLAB TESTMX TESTLAB TRNN TENN ADNN TRAINX TRAINL TESTX TESTL...
nLabels maxIters lambda L1neurons L2neurons epsInit sig sigrad sstf nP





% PREP NEURAL NET SYNAPTIC WEIGHTS (THETAS) & GRADIENT DECENT COST FUNC
%----------------------------------------------------------------------

% INITIALIZE THETA WEIGHTS
initTheta1 = rand(L2neurons, L1neurons+1) * 2 * epsInit - epsInit;
initTheta2 = rand(nLabels, L2neurons+1) * 2 * epsInit - epsInit;


% UNROLL THETA WEIGHTS 
initial_Thetas = [initTheta1(:) ; initTheta2(:)];
nn_Thetas = initial_Thetas;


% ESTABLISH COST FUNCTION
J = NNCostF(nn_Thetas, L1neurons, L2neurons, nLabels, TRAINX, TRAINL, lambda);








% TRAIN NEURAL NETWORK USING TRAINING DATASET
%----------------------------------------------------------------------
options = optimset('MaxIter', maxIters);

GcostFun = @(T) NNCostF(T, L1neurons, L2neurons, nLabels, TRAINX, TRAINL, lambda);

[nn_Thetas, cost] = NNfmincg(GcostFun, initial_Thetas, options);
%----------------------------------------------------------------------




% REROLL THETA WEIGHTS
Theta1 = reshape(nn_Thetas(1:L2neurons * (L1neurons + 1)), L2neurons, (L1neurons + 1));

Theta2 = reshape(nn_Thetas((1 + (L2neurons * (L1neurons + 1))):end), nLabels, (L2neurons + 1));




clearvars -except ADSP GENB LOCI CASE CTRL PHEN RES nn...
TRCASE TRCTRL TECASE TECTRL TRCASEN TRCTRLN TECASEN TECTRLN...
VLOCI VCASE VCTRL VTRCASE VTRCTRL VTECASE VTECTRL...
TRAINMX TRAINLAB TESTMX TESTLAB TRNN TENN ADNN TRAINX TRAINL TESTX TESTL...
nLabels maxIters lambda L1neurons L2neurons epsInit sig sigrad sstf nP...
Theta1 Theta2




% EVALUATE PERFORMANCE ON **TRAINING** SET
%----------------------------------------------------------------------

[p , a , h] = NNpredict(Theta1, Theta2, TRAINX);

% p : prediction label {1,2}
% a : confidence level of p
% h : confidence level of p and min label



TRAINPCTCORRECT = mean(p == TRAINL);

disp('Percent accuracy on training data:')
disp( TRAINPCTCORRECT )





% EVALUATE NEURAL NETWORK PERFORMANCE ON **TEST** SET
%----------------------------------------------------------------------

[p , a , h] = NNpredict(Theta1, Theta2, TESTX);


TESTPCTCORRECT = mean(p == TESTL);
disp('Percent accuracy on testing data:')
disp( TESTPCTCORRECT )






% EXAMINE & REPORT ACCURACY FOR HIGH CONFIDENCE DATA
%----------------------------------------------------------------------
MidConfThresh = .70;
HiConfThresh  = .85;


y = TESTL;

numCorrect = p == y;

MidConf      = a >= MidConfThresh & a < HiConfThresh;
MidConfPCT   = (mean(p(MidConf) == y(MidConf)))*100;
MidConfNp    = numel(p(MidConf));
MidConfPCTp  = (numel(p(MidConf)) / numel(p))*100;

MidHiConf      = a >= MidConfThresh;
MidHiConfPCT   = (mean(p(MidHiConf) == y(MidHiConf)))*100;
MidHiConfNp    = numel(p(MidHiConf));
MidHiConfPCTp  = (numel(p(MidHiConf)) / numel(p))*100;

HiConf       = a >= HiConfThresh;                   % Logical index of high confidence predictions
HiConfPCT    = (mean(p(HiConf) == y(HiConf)))*100;  % Percent correct hi conf predictions
HiConfNp     = numel(p(HiConf));                    % Total number of hi conf predictions
HiConfPCTp   = (numel(p(HiConf)) / numel(p))*100;   % Percent of hi conf predictions
HiConfNCorr  = HiConfNp * (HiConfPCT / 100);        % Total number of correct hi conf predictions


fprintf('\nPercent accuracy on MID conf (%.2f < a < %.2f) testing data: >>> %.1f%% <<<\n',...
    MidConfThresh, HiConfThresh, MidConfPCT)
fprintf('    (includes %.0f of %.0f [%.1f%%] Subset:Total patients) \n\n',...
    MidConfNp , numel(p) , MidConfPCTp )


fprintf('\nPercent accuracy on MID+HI conf (a > %.2f) testing data: >>> %.1f%% <<<\n',...
    MidConfThresh, MidHiConfPCT)
fprintf('    (includes %.0f of %.0f [%.1f%%] Subset:Total patients) \n\n',...
    MidHiConfNp , numel(p) , MidHiConfPCTp )


fprintf('\nPercent accuracy on HIGH conf (a > %.2f) testing data: >>> %.1f%% <<<\n',...
    HiConfThresh, HiConfPCT)
fprintf('    (includes %.0f of %.0f [%.1f%%] Subset:Total patients) \n\n',...
    HiConfNp , numel(p) , HiConfPCTp )





RES(nn,:) = [TESTPCTCORRECT MidHiConfPCT MidHiConfPCTp HiConfPCT HiConfPCTp];

clearvars -except ADSP GENB LOCI CASE CTRL PHEN RES nn...
TRCASE TRCTRL TECASE TECTRL TRCASEN TRCTRLN TECASEN TECTRLN...
VLOCI VCASE VCTRL VTRCASE VTRCTRL VTECASE VTECTRL...
TRAINMX TRAINLAB TESTMX TESTLAB TRNN TENN ADNN TRAINX TRAINL TESTX TESTL...
nLabels maxIters lambda L1neurons L2neurons epsInit sig sigrad sstf nP...
Theta1 Theta2

end

clc;
RES(:,1) = RES(:,1) .*100;
disp(RES)
disp(RES(:,2:end))





%% HOW DOES NN SCORE FAKE PERSON WITH 1 SINGLE VARIANT... IN THE TOP GENE
% THE TOP GENE IS USUALLY APOE


sig    = @(z) 1 ./ (1 + exp(-z) );      % sigmoid activation function
sigrad = @(g) sig(g) .* (1-sig(g));     % sigmoid gradient function
sstf   = @(n) 2 ./ (1 + exp(-2*n)) - 1; % sigmoid symetric transfer func

m = size(TRAINX, 1);
h1 = sig([ones(m, 1) TRAINX] * Theta1');
h2 = sig([ones(m, 1) h1] * Theta2');
[a, p] = max(h2, [], 2);

TRAINX(1,:)
h1(1,:)
h2(1,:)

FAKEPEOPLE = eye(size(TRAINX,2));

m = size(FAKEPEOPLE, 1);
h1 = sig([ones(m, 1) FAKEPEOPLE] * Theta1');
h2 = sig([ones(m, 1) h1] * Theta2');
[a, p] = max(h2, [], 2);

VLOCI.NNactiv  = h2(:,2);
VLOCI.NNguess  = p-1;
NNEYE = VLOCI(:,[7 26 27]);
disp(NNEYE(1:10,:))

