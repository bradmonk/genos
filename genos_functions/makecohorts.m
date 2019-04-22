function [TRCASE,TRCTRL,TECASE,TECTRL] = makecohorts(PHEN)
%%          PREPARE TRAINING AND TESTING SETS

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

ID CONSORT	STUDY_NAME  STUDY   FAM PEOPLE  ENR     CASE	CTRL	TOT    %CA %CO  EQL	CULL	NOTES   BRAAK
01  DGC     Adult Chng  ACT     0	1277	0       323     945     1268	25	75	323	-622	01 GOOD    1
02  ADGC	AD Centers  ADC     0	3263	0       2438	817     3255	75	25	817	1621	02 GOOD    0
03  CHARGE	Athrosclro  ARIC	0	57      0       39      18      57      68	32	18	21      03 MEH     0
04  CHARGE	Aus Stroke	ASPS	0	126     0       121     5       126     96	4	5	116     04 BAD     0
05  ADGC	Chic Aging  CHAP	0	231     0       27      204     231     12	88	27	-177	05 OKAY    0
06  CHARGE	CardioHlth  CHS     0	840     0       250     583     833     30	70	250	-333	06 GOOD    1
07  ADGC	Hispanic S  CUHS	361	0       331     160     171     331     48	52	160	-11     07 BAD     0
08  CHARGE	Erasmus Er  ERF     30	45      0       45      0       45      100	0	0	45      08 BAD     0
09  CHARGE	Framingham  FHS     0	581     0       157     424     581     27	73	157	-267	09 GOOD    1
10	ADGC	Gene Diffs  GDF     0	207     0       111     96      207     54	46	96	15      10 GOOD    1
11	ADGC	NIA - LOAD  LOAD	80  113     363     367     109     476     77	23	109	258     11 GOOD    1
12	ADGC	Aging Proj  MAP     0	415     0       138     277     415     33	67	138	-139	12 GOOD    1
13	ADGC	Mayo Clini  MAYO	0	349     0       250     99      349     72	28	99	151     13 GOOD    1
14	ADGC	Miami Univ  MIA     61	181     19      186     14      200     93	7	14	172     14 MEH     1
15	ADGC	AD Genetic  MIR     0	284     47      316     15      331     95	5	15	301     15 MEH     0
16	ADGC	Mayo cl PD  MPD     0	20      0       0       20      20      0	100	0	-20     16 BAD     1
17	ADGC	NationC AD  NCRD	18	108     52      160     0       160     100	0	0	160     17 BAD     1
18	ADGC	Wash Unive  RAS     0	0       46      46      0       46      100	0	0	46      18 BAD     1
19	ADGC	Relig Ordr  ROS     0	351     0       154     197     351     44	56	154	-43     19 GOOD    1
20	CHARGE	RotterdamS  RS      0	1089	0       276     813     1089	25	75	276	-537	20 GOOD    0
21	ADGC	Texas AD S  TARC	0	144     0       132     12      144     92	8	12	120     21 MEH     0
22	ADGC	Un Toronto  TOR     0	0       9       9       0       9       100	0	0	9       22 BAD     0
23	ADGC	Vanderbilt  VAN     6	235     1       210     26      236     89	11	26	184     23 MEH     1
24	ADGC	WashNY Age  WCAP	0	147     3       34      116     150     23	77	34	-82     24 OKAY    0
25	ADGC 	Univ Penns  UPN     40 	0       3       0       0       0       0	0	0	0       25 NONE    0




GOOD = [1 2 6 9 10 11 12 13 19 20];
OKAY = [3 14 15 21 23 24];
BAD  = [4 5 7 816 17 18 22];

Along the way, the code will take steps to mitigate stratification
issues due to COHORT-specific effects. 

...I'm trying to think if there is anything else you should know...
I can't, sooo, I guess we should just go-ahead and start with step #1.
(aka the first chance for something to go wrong)...

%}

%--------------------------------------------------------------------------
% LOAD THE DATASET
%{
close all; clear; clc; rng('shuffle');
genosdir = fileparts(which('GENOS.m'));
cd(genosdir);

subfuncpath = [genosdir filesep 'genosfunctions'];
datasetpath = [genosdir filesep 'genosdata'];
gpath = [genosdir pathsep subfuncpath pathsep datasetpath];
addpath(gpath)



which('ADSP_MINI3.mat')
ADSP = load('ADSP_MINI3.mat');



disp('dataset loaded')
clearvars -except ADSP
%}
%--------------------------------------------------------------------------

% LOCI = ADSP.LOCI;
% CASE = ADSP.CASE;
% CTRL = ADSP.CTRL;
% USNP = ADSP.USNP;
% PHEN = ADSP.PHEN;
% 
% 
% clearvars -except ADSP LOCI CASE CTRL PHEN USNP
% clc; disp(LOCI(1:9,:))



% keyboard


%% GENERATE TRAINING AND TESTING GROUPS (CONDITIONAL RANDOM SAMPLE)
% 1. Identify cohorts in PHEN with a reasonable number of CASE and CTRL
% 2. Take a random selection of people from those cohorts so that
%    there are the same number of CASE and CTRL.
% 3. The rest of the people from the keeper cohorts can be used to
%    increase the size of the NN test group.
% 
% % ID CONSORT	STUDY_NAME  STUDY   FAM PEOPLE  ENR     CASE	CTRL	TOT    %CA %CO  EQL	CULL   COH GOODCOH HASBRAAK
% % 01  DGC     Adult Chng  ACT     0	1277	0       323     945     1268	25	75	323	-622	01 1 GOOD    0:6 (257)
% % 02  ADGC	AD Centers  ADC     0	3263	0       2438	817     3255	75	25	817	1621	02 1 GOOD    0:6 (1861)
% % 03  CHARGE	Athrosclro  ARIC	0	57      0       39      18      57      68	32	18	21      03 0 MEH     -
% % 04  CHARGE	Aus Stroke	ASPS	0	126     0       121     5       126     96	4	5	116     04 0 BAD     -
% % 05  ADGC	Chic Aging  CHAP	0	231     0       27      204     231     12	88	27	-177	05 0 OKAY    -
% % 06  CHARGE	CardioHlth  CHS     0	840     0       250     583     833     30	70	250	-333	06 1 GOOD    -
% % 07  ADGC	Hispanic S  CUHS	361	0       331     160     171     331     48	52	160	-11     07 0 BAD     -
% % 08  CHARGE	Erasmus Er  ERF     30	45      0       45      0       45      100	0	0	45      08 0 BAD     -
% % 09  CHARGE	Framingham  FHS     0	581     0       157     424     581     27	73	157	-267	09 1 GOOD    1:6 (47)
% % 10	ADGC	Gene Diffs  GDF     0	207     0       111     96      207     54	46	96	15      10 1 GOOD    1:3 (6)
% % 11	ADGC	NIA - LOAD  LOAD	80  113     363     367     109     476     77	23	109	258     11 1 GOOD    1:6 (65)
% % 12	ADGC	Aging Proj  MAP     0	415     0       138     277     415     33	67	138	-139	12 1 GOOD    0:6 (230)
% % 13	ADGC	Mayo Clini  MAYO	0	349     0       250     99      349     72	28	99	151     13 1 GOOD    4:6 (42)
% % 14	ADGC	Miami Univ  MIA     61	181     19      186     14      200     93	7	14	172     14 0 MEH     1:6 (42)
% % 15	ADGC	AD Genetic  MIR     0	284     47      316     15      331     95	5	15	301     15 0 MEH     -
% % 16	ADGC	Mayo cl PD  MPD     0	20      0       0       20      20      0	100	0	-20     16 0 BAD     0:4 (20)
% % 17	ADGC	NationC AD  NCRD	18	108     52      160     0       160     100	0	0	160     17 0 BAD     2:6 (50)
% % 18	ADGC	Wash Unive  RAS     0	0       46      46      0       46      100	0	0	46      18 0 BAD     5:6 (9)
% % 19	ADGC	Relig Ordr  ROS     0	351     0       154     197     351     44	56	154	-43     19 1 GOOD    0:6 (245)
% % 20	CHARGE	RotterdamS  RS      0	1089	0       276     813     1089	25	75	276	-537	20 1 GOOD    -
% % 21	ADGC	Texas AD S  TARC	0	144     0       132     12      144     92	8	12	120     21 0 MEH     -
% % 22	ADGC	Un Toronto  TOR     0	0       9       9       0       9       100	0	0	9       22 0 BAD     -
% % 23	ADGC	Vanderbilt  VAN     6	235     1       210     26      236     89	11	26	184     23 0 MEH     6 (1)
% % 24	ADGC	WashNY Age  WCAP	0	147     3       34      116     150     23	77	34	-82     24 0 OKAY    -
% % 25	ADGC 	Univ Penns  UPN     40 	0       3       0       0       0       0	0	0	0       25 0 NONE    -


clearvars -except ADSP LOCI CASE CTRL PHEN USNP
clc; close all;





COHSET = PHEN;     % 8824  TOTAL PEOPLE
unique(PHEN.COHORTNUM)



% RANDOMIZE PHENOTYPE TABLE
%-------------------------------------------------------------
NVARS  = size(COHSET,1);        % Total number of people
k      = randperm(NVARS)';      % Get N random ints in range 1:N
COHSET  = COHSET(k,:);          % Scramble Phenotype table




% REMOVE ANY COHORT WITH LESS THAN 10 CASES AND 10 CTRLS
%-------------------------------------------------------------
% cohNums(minCACO < 10) = [];
% nCACO(minCACO < 10,:) = [];
% minCACO(minCACO < 10) = [];





PHETRCASE = COHSET( (COHSET.AD == 1) & (COHSET.COHORTNUM <  5) , :);
PHETRCTRL = COHSET( (COHSET.AD == 0) & (COHSET.COHORTNUM <  5) , :);
PHETECASE = COHSET( (COHSET.AD == 1) & (COHSET.COHORTNUM >= 5) , :);
PHETECTRL = COHSET( (COHSET.AD == 0) & (COHSET.COHORTNUM >= 5) , :);



%% ASSIGN A TRTE LABEL



PHETRCASE.TRTECACO = zeros(size(PHETRCASE,1),1) + 1;
PHETRCTRL.TRTECACO = zeros(size(PHETRCTRL,1),1) + 2;
PHETECASE.TRTECACO = zeros(size(PHETECASE,1),1) + 3;
PHETECTRL.TRTECACO = zeros(size(PHETECTRL,1),1) + 4;



%%
% RANDOMIZE PHENOTYPE TABLES OF EACH SUBGROUP
%-------------------------------------------------------------
NVARS      = size(PHETRCASE,1);     % Total number of people
k          = randperm(NVARS)';      % Get N random ints in range 1:N
PHETRCASE  = PHETRCASE(k,:);        % Scramble Phenotype table

NVARS      = size(PHETECASE,1);     % Total number of people
k          = randperm(NVARS)';      % Get N random ints in range 1:N
PHETECASE  = PHETECASE(k,:);        % Scramble Phenotype table

NVARS      = size(PHETRCTRL,1);     % Total number of people
k          = randperm(NVARS)';      % Get N random ints in range 1:N
PHETRCTRL  = PHETRCTRL(k,:);        % Scramble Phenotype table

NVARS      = size(PHETECTRL,1);     % Total number of people
k          = randperm(NVARS)';      % Get N random ints in range 1:N
PHETECTRL  = PHETECTRL(k,:);        % Scramble Phenotype table



%-------------------------------------------------------------
clc; clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
PHETRCASE PHETRCTRL PHETECASE PHETECTRL

disp(PHETRCASE(1:9,:))
disp('--------------------'); disp('N TRAINING & TESTING EXAMPLES')
disp(' '); fprintf('PHETRCASE... %.0f \n',size(PHETRCASE,1));
disp(' '); fprintf('PHETRCTRL... %.0f \n',size(PHETRCTRL,1));
disp(' '); fprintf('PHETECASE... %.0f \n',size(PHETECASE,1));
disp(' '); fprintf('PHETECTRL... %.0f \n',size(PHETECTRL,1));
disp('--------------------');
cohcounts(PHETRCASE,PHETRCTRL,PHETECASE,PHETECTRL,1)









%% PERFORM SOME CHECKS ON THE FINAL TRAINING AND TESTING GROUPS
%{.
clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
PHETRCASE PHETRCTRL PHETECASE PHETECTRL




%  AGEZ & BRAGEZ EXPLAINED
%--------------------------------------------------------------------------
%{

When training a classifier to recognize contributors to Alzheimers,
we want the classifier to know that young CASES and old CONTROLS
should be weighted more heavily than old CASES and young CONTROLS.

Furthermore, we want the classifier to know that CASES with high BRAAK
scores and CONTROLS with low BRAAK scores should be weighted more heavily.

Since weights need to be a single value per person, we compute a weight
called 'BRAGE' that is based on both Age and BRAAK scores, such that...

    BRAGE = AGE - BRAK

...the worse their 'BRAK' value the more their age is reduced (the first
'BRAK' value below is 5.4):

    BRAGE(BRAAK==0)   : AGE + 5.4
    BRAGE(BRAAK==1)   : AGE + 3.4
    BRAGE(BRAAK==2)   : AGE + 1.4
    BRAGE(BRAAK==NaN) : AGE - 0.0
    BRAGE(BRAAK==3)   : AGE - 2.5
    BRAGE(BRAAK==4)   : AGE - 5.5
    BRAGE(BRAAK==5)   : AGE - 7.5
    BRAGE(BRAAK==6)   : AGE - 9.5

So for example if someone is 80 years old and has a BRAAK==0 their
BRAGE equals 85.4; if they had a BRAAK==6 their BRAGE equals 70.5

From there a Z-score is computed for both AGE & BRAGE, but for CASES
the Z +/- sign is flipped...

    CASE_AGEZ   = zscore(CASE_AGE)   * -1
    CASE_BRAGEZ = zscore(CASE_BRAGE) * -1

    CTRL_AGEZ   = zscore(CTRL_AGE)   
    CTRL_BRAGEZ = zscore(CTRL_BRAGE)

the sign for CASES is flipped because the youngest cases with the highest
BRAAK scores will have the lowest BRAGE, but we want those CASES to have
the greatest weight, so the sign is flipped.

%}
%--------------------------------------------------------------------------


% COMPUTE Z-SCORE FOR CASE-AGE & CTRL-AGE, SEPARATELY
%------------------------------------------------
PHECASE = [PHETRCASE; PHETECASE];
PHECASE.AGE(PHECASE.AGE<55) = 55;
PHECASE.AGE = round(PHECASE.AGE);
PHECASE.AGEZ = zeros(size(PHECASE,1),1);
PHECASE.AGEZ = zscore(PHECASE.AGE);       % KEEP SIGN OF CASE GROUP Z

PHECTRL = [PHETRCTRL; PHETECTRL];
PHECTRL.AGE(PHECTRL.AGE<60) = 60;
PHECTRL.AGE = round(PHECTRL.AGE);
PHECTRL.AGEZ = zeros(size(PHECTRL,1),1);
PHECTRL.AGEZ = zscore(PHECTRL.AGE);       % REMEMBER TO FLIP SIGN OF CTRL GROUP Z
% PHECTRL.AGEZ = zscore(PHECTRL.AGE) .* -1; % REMEMBER TO FLIP SIGN OF CTRL GROUP Z
%------------------------------------------------



% CREATE A BRAAK-BASED AGE: BRAGE = AGE - BRAAK
%------------------------------------------------
CASENAN = isnan(PHECASE.BRAAK);
CTRLNAN = isnan(PHECTRL.BRAAK);

PHECASE.BRAK = CASENAN .* -1;      % MAKE NAN EQUAL -1 FOR NOW, ZERO LATER
PHECTRL.BRAK = CTRLNAN .* -1;      % MAKE NAN EQUAL -1 FOR NOW, ZERO LATER

PHECASE.BRAK(~CASENAN) = PHECASE.BRAAK(~CASENAN);
PHECTRL.BRAK(~CTRLNAN) = PHECTRL.BRAAK(~CTRLNAN);

PHECASE.BRAK(PHECASE.BRAK==0) = -6.4;
PHECTRL.BRAK(PHECTRL.BRAK==0) = -6.4;

PHECASE.BRAK(PHECASE.BRAK==1) = -4.4;
PHECTRL.BRAK(PHECTRL.BRAK==1) = -4.4;

PHECASE.BRAK(PHECASE.BRAK==2) = -2.4;
PHECTRL.BRAK(PHECTRL.BRAK==2) = -2.4;

PHECASE.BRAK(PHECASE.BRAK==3) = 1.5;
PHECTRL.BRAK(PHECTRL.BRAK==3) = 1.5;

PHECASE.BRAK(PHECASE.BRAK==4) = 4.5;
PHECTRL.BRAK(PHECTRL.BRAK==4) = 4.5;

PHECASE.BRAK(PHECASE.BRAK==5) = 6.5;
PHECTRL.BRAK(PHECTRL.BRAK==5) = 6.5;

PHECASE.BRAK(PHECASE.BRAK>=6) = 8.5;
PHECTRL.BRAK(PHECTRL.BRAK>=6) = 8.5;


PHECASE.BRAK = PHECASE.BRAK+1;      % ADD +1 BECAUSE NAN WAS MADE -1 ABOVE
PHECTRL.BRAK = PHECTRL.BRAK+1;      % ADD +1 BECAUSE NAN WAS MADE -1 ABOVE

PHECASE.BRAGE = PHECASE.AGE - PHECASE.BRAK;
PHECTRL.BRAGE = PHECTRL.AGE - PHECTRL.BRAK;
%------------------------------------------------


% COMPUTE Z-SCORE FOR CASE-AGE & CTRL-AGE, SEPARATELY
%------------------------------------------------
PHECASE.BRAGEZ = zeros(size(PHECASE,1),1);
PHECASE.BRAGEZ = zscore(PHECASE.BRAGE);       % FLIP SIGN OF CASE Z BELOW

PHECTRL.BRAGEZ = zeros(size(PHECTRL,1),1);
PHECTRL.BRAGEZ = zscore(PHECTRL.BRAGE);       % KEEP SIGN OF CTRL Z BELOW
%------------------------------------------------



% PLOT HISTOGRAM OF AGE AND AGE-TRANSFORMATIONS
%------------------------------------------------
close all;
fh01 = figure('Units','normalized','OuterPosition',[.01 .05 .95 .90],...
              'Color','w','MenuBar','none');
ax01 = axes('Position',[.04 .56 .2 .4],'Color','none');
ax02 = axes('Position',[.04 .06 .2 .4],'Color','none');
ax03 = axes('Position',[.29 .56 .2 .4],'Color','none');
ax04 = axes('Position',[.29 .06 .2 .4],'Color','none');
ax05 = axes('Position',[.54 .56 .2 .4],'Color','none');
ax06 = axes('Position',[.54 .06 .2 .4],'Color','none');
ax07 = axes('Position',[.79 .56 .2 .4],'Color','none');
ax08 = axes('Position',[.79 .06 .2 .4],'Color','none');

rbins=50.5:1:90.5;
zbins=-4:.25:4;

axes(ax01); histogram(PHECASE.AGE,rbins); xlabel('CASE AGE')
axes(ax02); histogram(PHECTRL.AGE,rbins); xlabel('CTRL AGE')
axes(ax03); histogram(PHECASE.AGEZ,zbins);  xlabel('CASE AGE Z-SCORE')
axes(ax04); histogram(PHECTRL.AGEZ,zbins);  xlabel('CTRL AGE Z-SCORE')
axes(ax05); histogram(PHECASE.BRAGE,rbins); xlabel('CASE BRAGE')
axes(ax06); histogram(PHECTRL.BRAGE,rbins); xlabel('CTRL BRAGE')
axes(ax07); histogram(PHECASE.BRAGEZ,zbins); xlabel('CASE BRAGE Z-SCORE')
axes(ax08); histogram(PHECTRL.BRAGEZ,zbins); xlabel('CTRL BRAGE Z-SCORE')
ax01.XLim = [55 95]; ax01.YLim = [0 1400];
ax02.XLim = [55 95]; ax02.YLim = [0 1400];
ax03.XLim = [-4  4]; ax03.YLim = [0 1400];
ax04.XLim = [-4  4]; ax04.YLim = [0 1400];
ax05.XLim = [55 95]; ax05.YLim = [0 1400];
ax06.XLim = [55 95]; ax06.YLim = [0 1400];
ax07.XLim = [-4  4]; ax07.YLim = [0 1400];
ax08.XLim = [-4  4]; ax08.YLim = [0 1400];

ymaxr = max([ax01.YLim(2) ax02.YLim(2) ax05.YLim(2) ax06.YLim(2)]);
ymaxz = max([ax03.YLim(2) ax04.YLim(2) ax07.YLim(2) ax08.YLim(2)]);

xminr = min([ax01.XLim(1) ax02.XLim(1) ax05.XLim(1) ax06.XLim(1)]);
xminz = min([ax03.XLim(1) ax04.XLim(1) ax07.XLim(1) ax08.XLim(1)]);
xmaxr = max([ax01.XLim(2) ax02.XLim(2) ax05.XLim(2) ax06.XLim(2)]);
xmaxz = max([ax03.XLim(2) ax04.XLim(2) ax07.XLim(2) ax08.XLim(2)]);




clearvars -except ADSP LOCI CASE CTRL PHEN USNP PHECASE PHECTRL


%###############################################################
%%   FLIP SIGN OF CONTROL BRAGEZ & AGEZ VARIABLES
%###############################################################


PHECASE.AGEZ   = zscore(PHECASE.AGE) .* -1; 

PHECASE.BRAGEZ = zscore(PHECASE.BRAGE) .* -1;




%###############################################################
%%   SAVE BACK INTO TRAINING AND TESTING SETS
%###############################################################


PHETRCASE = PHECASE(PHECASE.TRTECACO==1 ,:);
PHETRCTRL = PHECTRL(PHECTRL.TRTECACO==2 ,:);
PHETECASE = PHECASE(PHECASE.TRTECACO==3 ,:);
PHETECTRL = PHECTRL(PHECTRL.TRTECACO==4 ,:);




clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
PHETRCASE PHETRCTRL PHETECASE PHETECTRL



%###############################################################
%%  USE PER-PERSON VARIANT COUNTS TO IDENTIFY OUTLIERS
%###############################################################

Qlo = quantile(PHETRCASE.TOTvars,.01);
Qhi = quantile(PHETRCASE.TOTvars,.99);
TRCASE  = PHETRCASE(((PHETRCASE.TOTvars > Qlo) & (PHETRCASE.TOTvars < Qhi)),:);


Qlo = quantile(PHETRCTRL.TOTvars,.01);
Qhi = quantile(PHETRCTRL.TOTvars,.99);
TRCTRL  = PHETRCTRL(((PHETRCTRL.TOTvars > Qlo) & (PHETRCTRL.TOTvars < Qhi)),:);


Qlo = quantile(PHETECASE.TOTvars,.01);
Qhi = quantile(PHETECASE.TOTvars,.99);
TECASE  = PHETECASE(((PHETECASE.TOTvars > Qlo) & (PHETECASE.TOTvars < Qhi)),:);


Qlo = quantile(PHETECTRL.TOTvars,.01);
Qhi = quantile(PHETECTRL.TOTvars,.99);
TECTRL  = PHETECTRL(((PHETECTRL.TOTvars > Qlo) & (PHETECTRL.TOTvars < Qhi)),:);


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


disp('%##################################################################')
disp('                TRANING & VALIDATIONS SETS CREATED!                ')
disp('%##################################################################')
%%
end