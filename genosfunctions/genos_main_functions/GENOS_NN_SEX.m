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

% close all; clear; clc; rng('shuffle');
% genosdir = fileparts(which('GENOS.m'));
% cd(genosdir);
% 
% 
% subfuncpath = [genosdir '/genosfunctions'];
% datasetpath = [genosdir '/genosdata'];
% gpath = [genosdir ':' subfuncpath ':' datasetpath];
% addpath(gpath)
% 
% 
% which('ADSP.mat')
% load('ADSP.mat')
% 
% 
% clc; rng('shuffle');
% clearvars -except ADSP



%% CARBON COPY MAIN VARIABLES FROM ADSP.STRUCT


LOCI = ADSP.LOCI;
CASE = ADSP.CASE;
CTRL = ADSP.CTRL;
USNP = ADSP.UNSNP;
PHEN = ADSP.PHEN;


clearvars -except ADSP LOCI CASE CTRL PHEN USNP
clc; disp(LOCI(1:9,:))




%% IDENTIFY KEEPER PHEN AND MAKE CASE:CTRL EQUAL SIZE
% 1. Identify cohorts in PHEN with a reasonable number of CASE and CTRL
% 2. Take a random selection of people from those cohorts so that
%    there are the same number of CASE and CTRL.
% 3. The rest of the people from the keeper cohorts can be used to
%    increase the size of the NN test group.




% ONLY KEEP PEOPLE WITH A REASONABLE NUMBER OF VARIANTS
% PHE = PHEN(PHEN.TOTvars>14000,:);

PHE = PHEN;


% ADD 11 YEARS TO ALL CTRL AGES
PHE.AGE(PHE.AD==0) = PHE.AGE(PHE.AD==0) + 11;




% ONLY KEEP CTRLS ABOVE 98 YEARS OLD
PHE = PHE( PHE.AD==1  |  (PHE.AD==0 & PHE.AGE>98) , :);





% DESCRETIZE AGE
PHE.AGE = round(PHE.AGE);


PHE.AGE(PHE.AGE<55) = 55;

[Y,E] = discretize(PHE.AGE,[(min(PHE.AGE)-1) 72 80 91 (max(PHE.AGE)+1)]);
for nn = 1:(numel(E)-1)
PHE.AGE(Y==nn) = round((E(nn)+E(nn+1))/2);
end

disp('Number of Cases over 90 years old:')
sum(PHE.AGE>90 & PHE.AD==1)

disp('Number of Controls less than 91 years old:')
sum(PHE.AGE<91 & PHE.AD==0)



close all
histogram(PHE.AGE)
disp(PHE(1:9,:)); pause(1)

clearvars -except ADSP LOCI CASE CTRL PHEN USNP PHE






%###############################################################
%%  USE PER-PERSON VARIANT COUNTS TO IDENTIFY & REMOVE OUTLIERS
%###############################################################

Qlo = quantile(PHE.TOTvars,.01);
Qhi = quantile(PHE.TOTvars,.99);
PHE = PHE(((PHE.TOTvars > Qlo) & (PHE.TOTvars < Qhi)),:);


disp('------------------------------------')
disp('TOTAL VARIANTS PER-PERSON')
disp('      MIN    MAX')
fprintf('     %.0f  %.0f \n',min(PHE.TOTvars),  max(PHE.TOTvars))
disp('------------------------------------')



u = unique(PHE.AGE);


% close all;
% figure('Units','normalized','Position',[.05 .04 .9 .9],'Color','w','MenuBar','none');
% a1=subplot(2,3,1); histogram(PHE.TOTvars(PHE.AGE<=u(1)));                 title(num2str(u(1)))
% a2=subplot(2,3,2); histogram(PHE.TOTvars(PHE.AGE> u(1) & PHE.AGE<=u(2))); title(num2str(u(2)))
% a3=subplot(2,3,3); histogram(PHE.TOTvars(PHE.AGE> u(2) & PHE.AGE<=u(3))); title(num2str(u(3)))
% a4=subplot(2,3,4); histogram(PHE.TOTvars(PHE.AGE> u(3) & PHE.AGE<=u(4))); title(num2str(u(4)))
% a5=subplot(2,3,5); histogram(PHE.TOTvars(PHE.AGE> u(4)));                 title(num2str(u(5)))
% a=min([a1.XLim(1)  a2.XLim(1)  a3.XLim(1)  a4.XLim(1)  a5.XLim(1)]);
% b=max([a1.XLim(2)  a2.XLim(2)  a3.XLim(2)  a4.XLim(2)  a5.XLim(2)]);
% m=[a b];a1.XLim=m; a2.XLim=m;  a3.XLim=m;  a4.XLim=m;  a5.XLim=m;

close all;
figure('Units','normalized','Position',[.05 .04 .9 .9],'Color','w','MenuBar','none');
a1=subplot(2,2,1); histogram(PHE.TOTvars(PHE.AGE<=u(1)));                 title(num2str(u(1)))
a2=subplot(2,2,2); histogram(PHE.TOTvars(PHE.AGE> u(1) & PHE.AGE<=u(2))); title(num2str(u(2)))
a3=subplot(2,2,3); histogram(PHE.TOTvars(PHE.AGE> u(2) & PHE.AGE<=u(3))); title(num2str(u(3)))
a4=subplot(2,2,4); histogram(PHE.TOTvars(PHE.AGE> u(3)));                 title(num2str(u(4)))
a=min([a1.XLim(1)  a2.XLim(1)  a3.XLim(1)  a4.XLim(1)]);
b=max([a1.XLim(2)  a2.XLim(2)  a3.XLim(2)  a4.XLim(2)]);
m=[a b];a1.XLim=m; a2.XLim=m;  a3.XLim=m;  a4.XLim=m;







clearvars -except ADSP LOCI CASE CTRL PHEN USNP PHE















%% STORE VARIABLES FOR NEURAL NETWORK TRAINING AS 'VLOCI'
clc;

VLOCI     = LOCI;
VCASE     = CASE;
VCTRL     = CTRL;
VUSNP     = USNP;
VPHE      = PHE;

clearvars -except ADSP LOCI CASE CTRL PHEN USNP PHE ...
VLOCI VCASE VCTRL VUSNP VPHE



% SORT VARIANT TABLE BY FISHP
[~,i]  = sort(VLOCI.FISHP);
VLOCI  = VLOCI(i,:);
VCASE  = VCASE(i);
VCTRL  = VCTRL(i);
VUSNP  = VUSNP(i);
VLOCI.VID  = (1:size(VLOCI,1))';



% TAKE TOP N VARIANTS FOR NEURAL NET CLASSIFIER TRAINING
N = 1000;

VLOCI  = VLOCI(1:N,:);
VCASE  = VCASE(1:N);
VCTRL  = VCTRL(1:N);
VUSNP  = VUSNP(1:N);
VLOCI.VID  = (1:size(VLOCI,1))';



clc; disp(VLOCI(1:19,:))
fprintf('\n LOCI COUNT FOR FEATURE MATRIX: %.0f \n\n',size(VLOCI,1))

clearvars -except ADSP LOCI CASE CTRL PHEN USNP PHE ...
VLOCI VCASE VCTRL VUSNP VPHE





%###################################################################
%%      MAKE  RECTANGLE  NEURAL NET  VARIANT MATRIX
%###################################################################



% SCRAMBLE PHENOTYPE ORDER
P  = VPHE;
N  = size(P,1);      % Total number of people
k  = randperm(N)';   % Get N random ints in range 1:N
P  = P(k,:);         % Scramble Phenotype table
VPHE = P;



% MAKE THE NEURAL NET TRAINING MATRIX
[NNMX, vMX] = makeAgeNN(VLOCI,VCASE,VCTRL,VUSNP,VPHE,0);
disp(NNMX(1:10,1:10))


TRAIN = NNMX(2:end,:);
TRCHR = NNMX(1,5:end);


% HOLD-OUT A RANDOM SUBSET FOR VALIDATION
idx = randperm(size(TRAIN,1),2000);
TEST = TRAIN(idx,:);
TRAIN(idx,:) = [];



% CREATE FEATURE MATRIX AND LABEL MATRIX
TRAINX   = TRAIN(:,5:end)';
TRAINL   = TRAIN(:,2)';
TRAINALZ = TRAIN(:,3)';
TRAINCOH = TRAIN(:,4)';

TESTX   = TEST(:,5:end)';
TESTL   = TEST(:,2)';
TESTALZ = TEST(:,3)';
TESTCOH = TEST(:,4)';





disp(TRAINX(1:7,1:7))
disp(' ')
disp(TRAINL(1,1:7))

clearvars -except ADSP LOCI CASE CTRL PHEN USNP PHE ...
VLOCI VCASE VCTRL VUSNP VPHE NNMX...
TRAINX TRAINL TRAINAD TESTX TESTL TRAINALZ TESTALZ TRAINCOH TESTCOH





%% MAKE LABEL MATRICES HAVE N COLUMNS EQUAL TO N CLASSES


u = unique(TRAINL);
TRLABS = zeros(numel(u),size(TRAINL,2));
for i=1:size(TRAINL,2)
    TRLABS(:,i) = (u==TRAINL(i))';
end
TRLABS = TRLABS>0;


TELABS = zeros(numel(u),size(TESTL,2));
for i=1:size(TESTL,2)
    TELABS(:,i) = (u==TESTL(i))';
end
TELABS = TELABS>0;




clc; close all; clearvars -except ADSP LOCI CASE CTRL PHEN USNP PHE...
VLOCI VCASE VCTRL VUSNP VPHE NNMX TRAINX TRAINL TRAINAD TESTX TESTL...
TESTAD TRAINALZ TESTALZ TRAINCOH TESTCOH TRLABS TELABS


disp(' '); disp(TRLABS(:,1:7))
disp(' '); disp(TELABS(:,1:7))
% TEclass = vec2ind(TESTLOGIT);
% TEGUESS = round(TESTLOGIT);



% #######################################################################
%%         MACHINE LEARNING USING ARTIFICIAL NEURAL NETWORKS
% #######################################################################



imagesc(TRAINX)

NN = patternnet([300 100 50]);

net = train(NN,TRAINX,TRLABS);


% [NNPerf] = netstats(net,TRAINX,TRAINL,TESTX,TESTL,.05,.85,.15,1);


%% CHECK ROC CURVE PERFORMANCE
TROUT = sim(net,TRAINX);
TEOUT = sim(net,TESTX);

close all; 
fh1=figure('Units','normalized','Position',[.05 .05 .9 .9],'Color','w','MenuBar','none');
ax1=axes('Position',[.05 .09 .42 .83],'Color','none');
ax2=axes('Position',[.55 .09 .42 .83],'Color','none');
axes(ax1); plot((cell2mat((roc(TRLABS,TROUT))'))'); axis tight
line([0 size(TRLABS,2)],[0 1],'Color','k','LineStyle','--')
axes(ax2); plot((cell2mat((roc(TELABS,TEOUT))'))'); axis tight
line([0 size(TELABS,2)],[0 1],'Color','k','LineStyle','--')
pause(2)


clc; clearvars -except ADSP LOCI CASE CTRL PHEN USNP PHE...
VLOCI VCASE VCTRL VUSNP VPHE NNMX TRAINX TRAINL TRAINAD TESTX TESTL...
TESTAD TRAINALZ TESTALZ TRAINCOH TESTCOH TRLABS TELABS NN net NNPerf




%% CHECK PERFORMANCE OF THE TRAINED PATTERN-NET ON TEST DATASET

TRAINout = net(TRAINX);

TRAINLG = vec2ind(TRAINout);



TRAING = TRAINLG;
u = unique(TESTL);
for i = 1:numel(u)
    TRAING(TRAING==i) = u(i);
end

% GUESS = round(outputs);
% errors = gsubtract(TESTL,outputs);
% performance = perform(net,TESTL,outputs)


AGEjit = TRAINL + ((rand(numel(TRAINL),1)-.5).*5)';
GESjit = TRAING + ((rand(numel(TRAING),1)-.5).*5)';




clc; close all
fh1 = figure('Units','normalized','Position',[.05 .04 .6 .8],'Color','w');
ax1 = axes('Position',[.13 .13 .80 .80],'Color','none');


% scatter(GESjit, AGEjit, 50, TRAINL, '.');
scatter(GESjit, AGEjit, 100, TRAINCOH, '.');
colormap(jet)
xlabel({'\fontsize{22} AGE ', '\fontsize{16} (GUESS BIN)'});
ylabel({'\fontsize{22} AGE ', '\fontsize{16} (ACTUAL BIN)'});
ax1.XTick = unique(TRAINL);
ax1.YTick = unique(TRAINL);
ax1.FontSize = 16;




% GET AVERAGE GUESS
u = unique(TRAINL);
MUTRAIN = zeros(size(u));
for i = 1:numel(u)
    MUTRAIN(i) = mean(TRAING(TRAINL==u(i)));
end

hold on;
ph2=plot(MUTRAIN,u,'o-.');
ph2.MarkerSize=35;
ph2.LineWidth=4;
ph2.Color=[.01 .5 .5];
ph2.MarkerEdgeColor='none';
ph2.MarkerFaceColor=[.01 .5 .5];


















%% CHECK PERFORMANCE OF THE TRAINED PATTERN-NET ON TEST DATASET

TESTout = net(TESTX);

TESTLG = vec2ind(TESTout);



TESTG = TESTLG;
u = unique(TESTL);
for i = 1:numel(u)
    TESTG(TESTG==i) = u(i);
end

% GUESS = round(outputs);
% errors = gsubtract(TESTL,outputs);
% performance = perform(net,TESTL,outputs)


AGEjit = TESTL + ((rand(numel(TESTL),1)-.5).*5)';
GESjit = TESTG + ((rand(numel(TESTG),1)-.5).*5)';




clc; close all
fh1 = figure('Units','normalized','Position',[.05 .04 .6 .8],'Color','w');
ax1 = axes('Position',[.13 .13 .80 .80],'Color','none');

scatter(GESjit, AGEjit, 150, TESTL, '.');
colormap(jet)
xlabel({'\fontsize{22} AGE ', '\fontsize{16} (GUESS BIN)'});
ylabel({'\fontsize{22} AGE ', '\fontsize{16} (ACTUAL BIN)'});
ax1.XTick = unique(TESTL);
ax1.YTick = unique(TESTL);
ax1.FontSize = 16;




% GET AVERAGE GUESS
u = unique(TESTL);
MUTEST = zeros(size(u));
for i = 1:numel(u)
    MUTEST(i) = mean(TESTG(TESTL==u(i)));
end

hold on;
ph2=plot(MUTEST,u,'o-.');
ph2.MarkerSize=35;
ph2.LineWidth=4;
ph2.Color=[.01 .5 .5];
ph2.MarkerEdgeColor='none';
ph2.MarkerFaceColor=[.01 .5 .5];



return