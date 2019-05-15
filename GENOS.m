%% GENOS: 
%          Omics & machine learning tools for genomics datasets
%
% OVERVIEW
%--------------------------------------
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
%--------------------------------------


%% STEP-1: LOAD THE DATASET
%{.
close all; clear; clc; rng('shuffle');
genosdir = fileparts(which('GENOS.m'));
matlabdir = fileparts(which('MATLABS.m'));
cd(genosdir);


subfuncpath = [genosdir filesep 'genosfunctions'];
datasetpath = [genosdir filesep 'genosdata'];
gpath = [genosdir pathsep subfuncpath pathsep datasetpath];
addpath(gpath)


which('GENOSDATA.mat')
ADSP = load('GENOSDATA.mat');



disp('dataset loaded')
clearvars -except ADSP


%}




%% CARBON COPY MAIN VARIABLES FROM ADSP.STRUCT

LOCI = ADSP.LOCI;
CASE = ADSP.CASE;
CTRL = ADSP.CTRL;
USNP = ADSP.USNP;
PHEN = ADSP.PHEN;


clearvars -except ADSP LOCI CASE CTRL PHEN USNP
clc; disp(LOCI(1:9,:))




%% IDENTIFY KEEPER PHEN AND MAKE CASE:CTRL EQUAL SIZE
% 1. Identify cohorts in PHEN with a reasonable number of CASE and CTRL
% 2. Take a random selection of people from those cohorts so that
%    there are the same number of CASE and CTRL.
% 3. The rest of the people from the keeper cohorts can be used to
%    increase the size of the NN test group.




% GET SET OF GOOD COHORTS
%-------------------------------------------------------------
COHSET = PHEN(PHEN.GOODCOH==1,:);     % 8824  TOTAL PEOPLE



% RANDOMIZE PHENOTYPE TABLE
%-------------------------------------------------------------
NVARS  = size(COHSET,1);        % Total number of people
k      = randperm(NVARS)';      % Get N random ints in range 1:N
COHSET  = COHSET(k,:);          % Scramble Phenotype table






% COUNT NUMBER OF CASE AND CTRL IN EACH COHORT
%-------------------------------------------------------------
cohNums = unique(COHSET.COHORTNUM);

nCA=[];nCO=[];
for nn = 1:numel(cohNums)
    nCA(nn) = sum(  (COHSET.COHORTNUM==cohNums(nn)) & (COHSET.AD==1)  );
    nCO(nn) = sum(  (COHSET.COHORTNUM==cohNums(nn)) & (COHSET.AD==0)  );
end

nCACO = [nCA' nCO'];
minCACO = min(nCACO,[],2);





% REMOVE ANY COHORT WITH LESS THAN 10 CASES AND 10 CTRLS
%-------------------------------------------------------------
cohNums(minCACO < 10) = [];
nCACO(minCACO < 10,:) = [];
minCACO(minCACO < 10) = [];




% CREATE PHENOTYPE TABLE FOR EACH COHORT
%-------------------------------------------------------------
COHCASE={};COHCTRL={};
for nn = 1:numel(cohNums)
    COHCASE{nn} = COHSET( (COHSET.COHORTNUM==cohNums(nn)) & (COHSET.AD==1) ,:);
    COHCTRL{nn} = COHSET( (COHSET.COHORTNUM==cohNums(nn)) & (COHSET.AD==0) ,:);
end





% GET RANDOM PARTICIPANT SET FROM EACH COHORT
%
% Get random permutation of M values between 1:N-5
%    where  M = min CASE\CTRL group size per cohort
%           N = total cohort size
%-------------------------------------------------------------
rca={};rco={};
for nn = 1:numel(cohNums)
    rca{nn} = randperm( nCACO(nn,1)  , round(minCACO(nn) * .7));
    rco{nn} = randperm( nCACO(nn,2)  , round(minCACO(nn) * .7));
end





% GET ROW INDEX NUMBER FOR PEOPLE NOT CHOSEN ABOVE
%
% Get index locations for people not part of the training subgroup
%-------------------------------------------------------------
ica={};ico={};
for nn = 1:numel(cohNums)
    [ia,ib] = ismember(  (1:nCACO(nn,1)) , rca{nn});
    [ic,id] = ismember(  (1:nCACO(nn,2)) , rco{nn});
    ica{nn} = find(~ia);
    ico{nn} = find(~ic);
end




% CREATE PHEN TABLES FOR TRAINING CASE/CTRL & TESTING CASE/CTRL
%-------------------------------------------------------------
COHTRANCA={};COHTRANCO={};
COHTESTCA={};COHTESTCO={};
for nn = 1:numel(cohNums)
    COHTRANCA{nn} = COHCASE{nn}(rca{nn},:);
    COHTRANCO{nn} = COHCTRL{nn}(rco{nn},:);

    COHTESTCA{nn} = COHCASE{nn}(ica{nn},:);
    COHTESTCO{nn} = COHCTRL{nn}(ico{nn},:);
end



% TRAINING & TESTING VARS ABOVE ARE CELL ARRAYS OF PHEN TABLES
% OF THE SELECTED INDIVIDUALS REPRESENTING EACH COHORT. 
% HERE THE CODE MERGES THESE INTO A TABLE FOR:
% (1) TRAINING-CASE 
% (2) TRAINING-CTRL 
% (3) TESTING-CASE 
% (4) TESTING-CTRL
%-------------------------------------------------------------
PHETRCASE = vertcat(COHTRANCA{:});
PHETRCTRL = vertcat(COHTRANCO{:});

PHETECASE = vertcat(COHTESTCA{:});
PHETECTRL = vertcat(COHTESTCO{:});




% RANDOMIZE PHENOTYPE TABLE
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
cohcounts(PHETRCASE,PHETRCTRL,PHETECASE,PHETECTRL)





%% PUT SOME OF THE TEST GROUP BACK INTO THE TRAINING GROUP


szCA = round(size(PHETECASE,1)/2);
szCO = round(size(PHETECTRL,1)/2);

PHETRCASE = [PHETRCASE; PHETECASE(1:szCA,:)];
PHETRCTRL = [PHETRCTRL; PHETECTRL(1:szCO,:)];

PHETECASE(1:szCA,:) = [];
PHETECTRL(1:szCO,:) = [];


clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
PHETRCASE PHETRCTRL PHETECASE PHETECTRL

disp('--------------------'); disp('N TRAINING & TESTING EXAMPLES')
disp(' '); fprintf('PHETRCASE... %.0f \n',size(PHETRCASE,1));
disp(' '); fprintf('PHETRCTRL... %.0f \n',size(PHETRCTRL,1));
disp(' '); fprintf('PHETECASE... %.0f \n',size(PHETECASE,1));
disp(' '); fprintf('PHETECTRL... %.0f \n',size(PHETECTRL,1));
disp('--------------------');
cohcounts(PHETRCASE,PHETRCTRL,PHETECASE,PHETECTRL)








%% MAKE SURE EVERYONE CTRL.AGE>72  &  CASE.AGE<91


PHETRCTRL(PHETRCTRL.AGE < 72 , :) = [];
PHETECTRL(PHETECTRL.AGE < 72 , :) = [];

PHETRCASE(PHETRCASE.AGE > 91 , :) = [];
PHETECASE(PHETECASE.AGE > 91 , :) = [];



clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
PHETRCASE PHETRCTRL PHETECASE PHETECTRL

disp('--------------------'); disp('N TRAINING & TESTING EXAMPLES')
disp(' '); fprintf('PHETRCASE... %.0f \n',size(PHETRCASE,1));
disp(' '); fprintf('PHETRCTRL... %.0f \n',size(PHETRCTRL,1));
disp(' '); fprintf('PHETECASE... %.0f \n',size(PHETECASE,1));
disp(' '); fprintf('PHETECTRL... %.0f \n',size(PHETECTRL,1));
disp('--------------------');
cohcounts(PHETRCASE,PHETRCTRL,PHETECASE,PHETECTRL)





%###############################################################
%%  USE PER-PERSON VARIANT COUNTS TO IDENTIFY OUTLIERS
%###############################################################
%
% The total number of variants has already been counted for each
% person in the dataset and listed in the 'PHEN' table under
% 'TOTvars'. If you ever need to recount, a custom function has
% been coded to do this for you...
%
%   [VNUM] = countvperper(SRR, AD, COHORTNUM, CASE, CTRL);
%
% Here we will use the per-person counts already performed and
% dismiss the individuals on either extreme.
disp('COUNTING VARIANTS PER PERSON')


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



clc; close all;
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






%% DISCARD VARIANT LOCI WITH KNOWN ULTRA-RARE ALT ALLELE COUNTS

PASS = (LOCI.CASEALT > 20) | (LOCI.CTRLALT > 20);

LOCI  = LOCI(PASS,:);
CASE  = CASE(PASS);
CTRL  = CTRL(PASS);
USNP  = USNP(PASS);


clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
TRCASE TRCTRL TECASE TECTRL

% CASE = cellfun(@uint32,CASE,'UniformOutput',false);
% CTRL = cellfun(@uint32,CTRL,'UniformOutput',false);
% USNP = cellfun(@uint32,USNP,'UniformOutput',false);
% 
% clear ADSP
% 
% ADSP.PHEN = PHEN;
% ADSP.LOCI = LOCI;
% ADSP.CASE = CASE;
% ADSP.CTRL = CTRL;
% ADSP.USNP = USNP;
% 
% save('ADSP_MINI.mat','ADSP');





%###############################################################
%%          COUNT NUMBER OF VARIANTS PER LOCI
%###############################################################
%
% The varsum() function will go through each known variant loci
% and check whether anyone's SRR ID from your subset of IDs match
% all known SRR IDs for that loci. It will then sum the total
% number of alleles (+1 for hetzy-alt, +2 for homzy-alt) for each
% loci and return the totals.
disp('COUNTING VARIANTS PER LOCI IN EACH TREATMENT GROUP')



[TRCASEN,TRCTRLN,TECASEN,TECTRLN,TRCASEUN,TRCTRLUN,TECASEUN,TECTRLUN] = ...
    snpsum(CASE,CTRL,USNP,TRCASE.SRR,TRCTRL.SRR,TECASE.SRR,TECTRL.SRR);


% SAVE COUNTS **PER CHROMOSOME**
LOCI.TRCASEREF = (numel(TRCASE.SRR)*2) - (TRCASEUN.*2) - TRCASEN;
LOCI.TRCTRLREF = (numel(TRCTRL.SRR)*2) - (TRCTRLUN.*2) - TRCTRLN;
LOCI.TRCASEALT = TRCASEN;
LOCI.TRCTRLALT = TRCTRLN;


LOCI.TECASEREF = (numel(TECASE.SRR)*2) - (TECASEUN.*2) - TECASEN;
LOCI.TECTRLREF = (numel(TECTRL.SRR)*2) - (TECTRLUN.*2) - TECTRLN;
LOCI.TECASEALT = TECASEN;
LOCI.TECTRLALT = TECTRLN;




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

pause(2);



% GET VARIANT COUNTS PER LOCI USING **PERSON INSTEAD OF CHROMOSOME**
[TRCASEN,TRCTRLN,TECASEN,TECTRLN,TRCASEUN,TRCTRLUN,TECASEUN,TECTRLUN] = ...
    ppsnpsum(CASE,CTRL,USNP,TRCASE.SRR,TRCTRL.SRR,TECASE.SRR,TECTRL.SRR);

% SAVE COUNTS **PER PERSON**
LOCI.PPTRCASEREF = (numel(TRCASE.SRR)) - (TRCASEUN) - TRCASEN;
LOCI.PPTRCTRLREF = (numel(TRCTRL.SRR)) - (TRCTRLUN) - TRCTRLN;
LOCI.PPTRCASEALT = TRCASEN;
LOCI.PPTRCTRLALT = TRCTRLN;


LOCI.PPTECASEREF = (numel(TECASE.SRR)) - (TECASEUN) - TECASEN;
LOCI.PPTECTRLREF = (numel(TECTRL.SRR)) - (TECTRLUN) - TECTRLN;
LOCI.PPTECASEALT = TECASEN;
LOCI.PPTECTRLALT = TECTRLN;


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








%########################################################################
%%               COMPUTE FISHER'S P-VALUE
%########################################################################
disp('COMPUTING FISHERS EXACT TEST STATISTICS PER DNA LOCUS')



%------------------------------------------------------------------------
% ***  PER CHROMOSOME  ***
%---------------------------

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





%------------------------------------------------------------------------
% ***  PER PERSON  ***
%-----------------------
% 
% % COMPUTE FISHERS STATISTICS FOR THE TRAINING GROUP
% [PPFISHP, PPFISHOR] = fishp_mex(LOCI.PPTRCASEREF,LOCI.PPTRCASEALT,...
%                             LOCI.PPTRCTRLREF,LOCI.PPTRCTRLALT);
% 
% LOCI.PPTRFISHP  = PPFISHP;
% LOCI.PPTRFISHOR = PPFISHOR;
% 
% 
% 
% % COMPUTE FISHERS STATISTICS FOR THE TRAINING GROUP
% [PPFISHP, PPFISHOR] = fishp_mex(LOCI.PPTECASEREF, LOCI.PPTECASEALT,...
%                             LOCI.PPTECTRLREF, LOCI.PPTECTRLALT);
% 
% LOCI.PPTEFISHP  = PPFISHP;
% LOCI.PPTEFISHOR = PPFISHOR;
%------------------------------------------------------------------------







%------------------------------------------------------------------------
% PLOT HISTOGRAM OF P-VALUES
%----------------------------

a = LOCI.TRFISHP;
b = LOCI.TEFISHP;
qa = quantile(a,[.001,.999]);
qb = quantile(b,[.001,.999]);

x = -log(a);
y = -log(b);
qx = quantile(x,[.001,.999]);
qy = quantile(y,[.001,.999]);

x = x(x>qx(1)&x<qx(2));
y = y(y>qy(1)&y<qy(2));



close all
fh01 = figure('Units','normalized','OuterPosition',[.01 .05 .95 .90],...
              'Color','w','MenuBar','none');
ax01 = axes('Position',[.06 .56 .4 .4],'Color','none');
ax02 = axes('Position',[.56 .56 .4 .4],'Color','none');
ax03 = axes('Position',[.06 .06 .4 .4],'Color','none');
ax04 = axes('Position',[.56 .06 .4 .4],'Color','none');

axes(ax01); histogram(a);
title('TRAINING SET FISHERS P-VALUE DISTRIBUTION')
axes(ax02); histogram(b);
title('TESTING SET FISHERS P-VALUE DISTRIBUTION')

axes(ax03); histogram(x);
title('TRAINING SET quantile(-log(FISHP),[.001 .999])')
axes(ax04); histogram(y);
title('TESTING SET quantile(-log(FISHP),[.001 .999])')
pause(2); close all;


clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
TRCASE TRCTRL TECASE TECTRL



%########################################################################
%% COMPUTE CHI SQUARE VALUE
%########################################################################
%{

%------------------------------------------------------------------------
% ***  PER CHROMOSOME  ***
%-----------------------------

% COMPUTE CHI SQUARED STATISTICS FOR THE TRAINING GROUP **per chromosome**
[CHIP, CHIOR] = chisq(LOCI.TRCASEREF,LOCI.TRCASEALT,...
                      LOCI.TRCTRLREF,LOCI.TRCTRLALT);

LOCI.TRCHIP  = CHIP;
LOCI.TRCHIOR = CHIOR;



% COMPUTE CHI SQUARED STATISTICS FOR THE TRAINING GROUP **per chromosome**
[CHIP, CHIOR] = chisq(LOCI.TECASEREF, LOCI.TECASEALT,...
                      LOCI.TECTRLREF, LOCI.TECTRLALT);

LOCI.TECHIP  = CHIP;
LOCI.TECHIOR = CHIOR;



%------------------------------------------------------------------------
% ***  PER PERSON  ***
%-----------------------------
% 
% % COMPUTE CHI SQUARED STATISTICS FOR THE TRAINING GROUP **per person**
% [CHIP, CHIOR] = chisq(LOCI.PPTRCASEREF,LOCI.PPTRCASEALT,...
%                       LOCI.PPTRCTRLREF,LOCI.PPTRCTRLALT);
% 
% LOCI.PPTRCHIP  = CHIP;
% LOCI.PPTRCHIOR = CHIOR;
% 
% 
% 
% % COMPUTE CHI SQUARED STATISTICS FOR THE TRAINING GROUP **per person**
% [CHIP, CHIOR] = chisq(LOCI.PPTECASEREF, LOCI.PPTECASEALT,...
%                       LOCI.PPTECTRLREF, LOCI.PPTECTRLALT);
% 
% LOCI.PPTECHIP  = CHIP;
% LOCI.PPTECHIOR = CHIOR;
% 
%------------------------------------------------------------------------





% PLOT HISTOGRAM OF CHI SQUARE P-VALUES
%------------------------------------

a = LOCI.TRCHIP;
b = LOCI.TECHIP;
qa = quantile(a,[.001,.999]);
qb = quantile(b,[.001,.999]);

x = -log(a);
y = -log(b);
qx = quantile(x,[.001,.999]);
qy = quantile(y,[.001,.999]);

x = x(x>qx(1)&x<qx(2));
y = y(y>qy(1)&y<qy(2));



close all
fh01 = figure('Units','normalized','OuterPosition',[.01 .05 .95 .90],...
              'Color','w','MenuBar','none');
ax01 = axes('Position',[.06 .56 .4 .4],'Color','none');
ax02 = axes('Position',[.56 .56 .4 .4],'Color','none');
ax03 = axes('Position',[.06 .06 .4 .4],'Color','none');
ax04 = axes('Position',[.56 .06 .4 .4],'Color','none');

axes(ax01); histogram(a);
title('TRAINING SET CHISQ P-VALUE DISTRIBUTION')
axes(ax02); histogram(b);
title('TESTING SET CHISQ P-VALUE DISTRIBUTION')

axes(ax03); histogram(x);
title('TRAINING SET quantile(-log(CHISQ),[.001 .999])')
axes(ax04); histogram(y);
title('TESTING SET quantile(-log(CHISQ),[.001 .999])')
pause(2)



fh1 = figure('Units','pixels','Position',[10 35 1100 750],...
    'Color','w','MenuBar','none');
ax1 = axes('Position',[.06 .06 .9 .9],'Color','none');

ph1 = scatter(LOCI.TRCHIP,LOCI.TRFISHP,'.');
xlim([0 1]); ylim([0 1])
xlabel('\fontsize{16} Chi Square P-value')
ylabel('\fontsize{16} Fishers Exact P-value')
pause(2)


% ph1 = scatter(LOCI.PPTRFISHP,LOCI.TRFISHP,'.');
% xlim([0 1]); ylim([0 1])
% xlabel('\fontsize{16} Fishers Exact P-value PER PERSON')
% ylabel('\fontsize{16} Fishers Exact P-value PER CHROMOSOME')
% pause(2)


clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
TRCASE TRCTRL TECASE TECTRL

%}







%########################################################################
%% PLOT SCATTER GRAPH COMPARISONS OF TRAINING & TESTING GROUPS
%########################################################################
%{
close all
fh01 = figure('Units','normalized','OuterPosition',[.01 .05 .95 .90],...
              'Color','w','MenuBar','none');
ax01 = axes('Position',[.06 .56 .4 .4],'Color','none');
ax02 = axes('Position',[.56 .56 .4 .4],'Color','none');
ax03 = axes('Position',[.06 .06 .4 .4],'Color','none');
ax04 = axes('Position',[.56 .06 .4 .4],'Color','none');

axes(ax01); scatter(LOCI.TRFISHP,LOCI.TEFISHP,'.');
xlim([0 1]); ylim([0 1])
title('TRAINING VS TESTING FISHERS P-VALUE')

axes(ax02); scatter(LOCI.TRCHIP,LOCI.TECHIP,'.');
xlim([0 1]); ylim([0 1])
title('TRAINING VS TESTING CHISQ P-VALUE')

axes(ax03); scatter(LOCI.PPTRFISHP,LOCI.PPTEFISHP,'.');
xlim([0 1]); ylim([0 1])
title('TRAINING VS TESTING FISHERS P-VALUE')

axes(ax04); scatter(LOCI.PPTRCHIP,LOCI.PPTECHIP,'.');
xlim([0 1]); ylim([0 1])
title('TRAINING VS TESTING CHISQ P-VALUE')

pause(2)
%}




return
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
clc; clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL

VLOCI     = LOCI;
VCASE     = CASE;
VCTRL     = CTRL;
VUSNP     = USNP;
VTRCASE   = TRCASE;
VTRCTRL   = TRCTRL;
VTECASE   = TECASE;
VTECTRL   = TECTRL;


% SET MAIN FISHP TO TRAINING GROUP FISHP
VLOCI.FISHP      = VLOCI.TRFISHP;
VLOCI.FISHOR     = VLOCI.TRFISHOR;
VLOCI.CASEREF    = VLOCI.TRCASEREF;
VLOCI.CASEALT    = VLOCI.TRCASEALT;
VLOCI.CTRLREF    = VLOCI.TRCTRLREF;
VLOCI.CTRLALT    = VLOCI.TRCTRLALT;



% SORT VARIANTS BY EITHER TRFISHP|CHRPOS
[~,i]  = sort(VLOCI.CHRPOS);
% [~,i]  = sort(VLOCI.TRFISHP);
VLOCI  = VLOCI(i,:);
VCASE  = VCASE(i);
VCTRL  = VCTRL(i);
VUSNP  = VUSNP(i);


clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL






%% REMOVE SELECT GENES BY NAME
% REMOVE KNOWN AD ALZ GENES
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
%     VCASE(x)   = [];
%     VCTRL(x)   = [];
%     VUSNP(x)   = [];
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
%     VCASE(x)   = [];
%     VCTRL(x)   = [];
%     VUSNP(x)   = [];
% 
% end
% VLOCI.VID  = (1:size(VLOCI,1))';
% 
% clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
% VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL
%}



%% CREATE MANHATTAN PLOT OF FISHER'S P-VALUES
% SEE GENOS_HAPLOTYPE.m
%{
clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL


% SORT VARIANTS BY FISHERS P-VALUE
[~,i]  = sort(VLOCI.TRFISHP);
% SORT VARIANTS BY CHRPOS
[~,i]  = sort(VLOCI.CHRPOS);
VLOCI  = VLOCI(i,:);
VCASE  = VCASE(i);
VCTRL  = VCTRL(i);
VUSNP  = VUSNP(i);


clc; close all;
fh1 = figure('Units','normalized','OuterPosition',[.01 .05 .9 .9],...
             'Color','w','MenuBar','none');
ax1 = axes('Position',[.05 .55 .9 .4],'Color','none',...
           'XColor','none','XLabel',[]); hold on;
ax2 = axes('Position',[.05 .05 .9 .4],'Color','none',...
           'XColor','none','XLabel',[]); hold on;
% ax1 = axes('Position',[.06 .05 .93 .9],'Color','none',...
% 'XColor','none','XLabel',[]); hold on;




axes(ax1)


AH = -log10(VLOCI.TRFISHP);
PAH = AH > 20;
AH(PAH) = 20-rand*2;
ph1 = scatter(1:size(VLOCI,1),  AH  ,300, VLOCI.CHR , '.');
% xlabel('Chr:Loc')
ylabel('\fontsize{22}-log(P)')
ax1.FontSize = 16;
% yt = ax1.YTickLabel;
% ax1.YTickLabel = ['\fontsize{24}' yt];
% title('Manhattan Plot of Fishers Exact Test Statistic')
bonfapha = .05/size(VLOCI,1); % -log10(5.9e-08)
line(ax1,[1 size(VLOCI,1)],-log10([bonfapha bonfapha]),...
'Color','red','LineStyle','--','LineWidth',2)
axis tight
ax1.YLim = [0 20];



axes(ax2)

MH = -log10(VLOCI.TEFISHP);
PMH = MH > 20;
MH(PMH) = 20-rand*2;
ph2 = scatter(1:size(VLOCI,1),  MH  ,300, VLOCI.CHR , '.');
% xlabel('Chr:Loc')
ylabel('\fontsize{22}-log(P)')
ax2.FontSize = 16;
% title('Manhattan Plot of Fishers Exact Test Statistic')
bonfapha = .05/size(VLOCI,1); % -log10(5.9e-08)
line(ax2,[1 size(VLOCI,1)],-log10([bonfapha bonfapha]),...
'Color','red','LineStyle','--','LineWidth',2)
axis tight
ax2.YLim = [0 20];


clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL fh1 ax1 ax2



%---- ADD GENE NAME TEXT LABELS TO GRAPH

nL = 90;
TRnlpCUTOFF = -log10(.00000001);
TEnlpCUTOFF = -log10(.0001);

VLOCI.ORD = (1:size(VLOCI,1))';


% SORT TABLE BY P-VALUE
[i,j] = sort(VLOCI.TRFISHP);
GENOGENE = VLOCI.GENE(j,:);
GENOORD  = VLOCI.ORD(j,:);
GENOP    = -log10(VLOCI.TRFISHP(j,:));
I = GENOP>20;
K = (exp(exp(log(log(GENOP(I)))))./20-1)./5;
GENOP(I) = 19 + K;


[i,j] = sort(VLOCI.TEFISHP);
GENODGENE = VLOCI.GENE(j,:);
GENODORD  = VLOCI.ORD(j,:);
GENODP    = -log10(VLOCI.TEFISHP(j,:));
I = GENODP>20;
K = exp(exp(log(log(GENODP(I)))))./20-1;
GENODP(I) = 19 + K;


TRi = GENOP  > TRnlpCUTOFF;
TEi = GENODP > TEnlpCUTOFF;


Gi = GENOGENE(TRi,:);
Xi = GENOORD(TRi,:);
Yi = GENOP(TRi,:);

Gj = GENODGENE(TEi,:);
Xj = GENODORD(TEi,:);
Yj = GENODP(TEi,:);

axes(ax1)
text(Xi , Yi , char(Gi) , 'HorizontalAlignment','right');

axes(ax2)
text(Xj , Yj , char(Gj) , 'HorizontalAlignment','right');



clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL...
fh1 ax1 ax2 Gi Xi Yi Gj Xj Yj





% SORT VARIANTS BY FISHERS/CHRPOS
% [~,i]  = sort(VLOCI.TRFISHP);
[~,i]  = sort(VLOCI.CHRPOS);
VLOCI  = VLOCI(i,:);
VCASE  = VCASE(i);
VCTRL  = VCTRL(i);
VUSNP  = VUSNP(i);

%---------- EXPORT IMAGE TO DESKTOP FOLDER -----------------------------
dt=char(datetime(datetime,'Format','yyyy-MM-dd-HH-mm-ss'));
print(['/Users/bradleymonk/Desktop/MANHATTAN ' dt '.png'],'-dpng');
%-----------------------------------------------------------------------



%------------------------------------------------------------------------
% 
%%                 [[[[  RESET BUTTON    ]]]]
% ------------------------------------------------------------------------
VLOCI     = LOCI;
VCASE     = CASE;
VCTRL     = CTRL;
VUSNP     = USNP;
VTRCASE   = TRCASE;
VTRCTRL   = TRCTRL;
VTECASE   = TECASE;
VTECTRL   = TECTRL;


% SET MAIN FISHP TO TRAINING GROUP FISHP
VLOCI.FISHP      = VLOCI.TRFISHP;
VLOCI.FISHOR     = VLOCI.TRFISHOR;
VLOCI.CASEREF    = VLOCI.TRCASEREF;
VLOCI.CASEALT    = VLOCI.TRCASEALT;
VLOCI.CTRLREF    = VLOCI.TRCTRLREF;
VLOCI.CTRLALT    = VLOCI.TRCTRLALT;


% % SORT VARIANTS BY FISHERS P-VALUE
% [~,i]  = sort(VLOCI.TRFISHP);
% VLOCI  = VLOCI(i,:);
% VCASE  = VCASE(i);
% VCTRL  = VCTRL(i);
% VUSNP  = VUSNP(i);


clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL








%% DETERMINE HAPLOTYPES WITH LINKAGE DISEQUILIBRIUM

clc; clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL

% SORT VARIANTS BY FISHERS P-VALUE
% [~,i]  = sort(VLOCI.TRFISHP);
i = VLOCI.TRFISHP < .005;
VLOCI  = VLOCI(i,:);
VCASE  = VCASE(i);
VCTRL  = VCTRL(i);
VUSNP  = VUSNP(i);



% COMPUTE WHETHER P-VALUES UPSTREAM AND DOWNSTREAM ARE DECREASING
PPOS = VLOCI.TRFISHP(2:end-1);
PPUP = VLOCI.TRFISHP(3:end);
PPDN = VLOCI.TRFISHP(1:end-2);


PREVP = PPOS < PPUP;
NEXTP = PPOS < PPDN;
VLOCI.PREVP = [0;PREVP;0];
VLOCI.NEXTP = [0;NEXTP;0];
VLOCI.MINIMA = VLOCI.PREVP & VLOCI.NEXTP;

VLOCI.MINIMA(1) = VLOCI.TRFISHP(1) < VLOCI.TRFISHP(2);
VLOCI.MINIMA(end) = VLOCI.TRFISHP(end) < VLOCI.TRFISHP(end-1);


LOLOC = VLOCI(VLOCI.MINIMA==1,:);





%% COMPUTE PROXIMITY TO UPSTREAM AND DOWNSTREAM VARIANT
% LPOS = LOC.CHRPOS(2:end-1);
% LPUP = LOC.CHRPOS(3:end);
% LPDN = LOC.CHRPOS(1:end-2);


APOE_TOMM40_DIST = 412079 - 395714;

% SORT BY CHRPOS
[~,i]  = sort(VLOCI.CHRPOS);
VLOCI  = VLOCI(i,:);
VCASE  = VCASE(i);
VCTRL  = VCTRL(i);
VUSNP  = VUSNP(i);



POS = VLOCI.POS(2:end-1);
NEXTPOS = VLOCI.POS(3:end);
PREVPOS = VLOCI.POS(1:end-2);

DISTPREV = POS - PREVPOS;
DISTNEXT = NEXTPOS - POS;

VLOCI.DISTPREV       = [0;DISTPREV;0];
VLOCI.DISTNEXT       = [0;DISTNEXT;0];
VLOCI.DISTNEXT(1)    = VLOCI.POS(2)-VLOCI.POS(1);

VLOCI.DISTPREV(end)  = VLOCI.POS(end)-VLOCI.POS(end-1);
VLOCI.DISTNEXT(1)    = VLOCI.POS(2)-VLOCI.POS(1);






OK1 = (VLOCI.DISTPREV < APOE_TOMM40_DIST);
OK2 = (VLOCI.DISTNEXT < APOE_TOMM40_DIST);
BAD = BAD1 | BAD2;

VLOCI(BAD,:) = [];
VCASE(i) = [];
VCTRL(i) = [];
VUSNP(i) = [];

% SORT BY CHRPOS
[~,i]  = sort(VLOCI.CHRPOS);
VLOCI  = VLOCI(i,:);
VCASE  = VCASE(i);
VCTRL  = VCTRL(i);
VUSNP  = VUSNP(i);



BAD = (VLOCI.DISTPREV < APOE_TOMM40_DIST) & (VLOCI.MINIMA == 0);
VLOCI(BAD,:) = [];
VCASE(i) = [];
VCTRL(i) = [];
VUSNP(i) = [];

% SORT BY CHRPOS
[~,i]  = sort(VLOCI.CHRPOS);
VLOCI  = VLOCI(i,:);
VCASE  = VCASE(i);
VCTRL  = VCTRL(i);
VUSNP  = VUSNP(i);






close all
% subplot(2,1,1); histogram(-log(1./double(LOC.DISTNEXT)))
% line([-log(1/600) -log(1/600)],[0 3000])


figure
subplot(2,1,1); histogram(log(double(VLOCI.DISTNEXT)))
% line([log(600) log(600)],[0 3000])
subplot(2,1,2); histogram(double(VLOCI.DISTNEXT(VLOCI.DISTNEXT<1e5)))






figure
subplot(2,1,1); histogram(log(double(LOX.DISTNEXT)))
% line([log(600) log(600)],[0 3000])
subplot(2,1,2); histogram(double(LOX.DISTNEXT(LOX.DISTNEXT<1e5)))




% SORT VARIANTS BY FISHERS P-VALUE
[~,i]  = sort(VLOCI.TRFISHP);
VLOCI  = VLOCI(i,:);
VCASE  = VCASE(i);
VCTRL  = VCTRL(i);
VUSNP  = VUSNP(i);


[~,i]  = sort(LOX.CHRPOS);
LOX    = LOX(i,:);
VCASE  = VCASE(i);
VCTRL  = VCTRL(i);
VUSNP  = VUSNP(i);



%}




%% ONLY KEEP ONE OR TWO VARIANTS PER GENE (THE ONES WITH LOWEST P-VALS)
%{

% SORT VARIANTS BY FISHERS P-VALUE
[~,i]  = sort(VLOCI.TRFISHP);
VLOCI  = VLOCI(i,:);
VCASE  = VCASE(i);
VCTRL  = VCTRL(i);
VUSNP  = VUSNP(i);



% GET FIRST TOP RATED VARIANT FOR EACH GENE
[~,i,~] = unique(string(VLOCI.GENE),'stable');

VLO1  = VLOCI(i,:);
VCA1  = VCASE(i);
VCO1  = VCTRL(i);
VUS1  = VUSNP(i);

VLOCI(i,:) = [];
VCASE(i)   = [];
VCTRL(i)   = [];
VUSNP(i)   = [];


% GET SECOND TOP RATED VARIANT FOR EACH GENE
[~,i,~] = unique(string(VLOCI.GENE),'stable');

VLO2  = VLOCI(i,:);
VCA2  = VCASE(i);
VCO2  = VCTRL(i);
VUS2  = VUSNP(i);

VLOCI = [VLO1 ; VLO2];
VCASE = [VCA1 ; VCA2];
VCTRL = [VCO1 ; VCO2];
VUSNP = [VUS1 ; VUS2];



% SORT VARIANTS BY FISHERS P-VALUE
[~,i]  = sort(VLOCI.TRFISHP);
VLOCI  = VLOCI(i,:);
VCASE  = VCASE(i);
VCTRL  = VCTRL(i);
VUSNP  = VUSNP(i);




clc; disp(VLOCI(1:10,:))
clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL
%}











%% TAKE THE TOP N or P<x GENES FOR NEURAL NET CLASSIFIER TRAINING



% SORT VARIANTS BY EITHER TRFISHP|CHRPOS
%[~,i]  = sort(VLOCI.CHRPOS);
[~,i]  = sort(VLOCI.TRFISHP);
VLOCI  = VLOCI(i,:);
VCASE  = VCASE(i);
VCTRL  = VCTRL(i);
VUSNP  = VUSNP(i);




Vpass = VLOCI.TRFISHP < .001;
SNPn = sum(Vpass);

% EXTRACT TOP-N NUMBER OF VARIANTS
VLOCI  = VLOCI(1:SNPn,:);
VCASE  = VCASE(1:SNPn);
VCTRL  = VCTRL(1:SNPn);
VUSNP  = VUSNP(1:SNPn);



% disp(VLOCI(1:9,:))
fprintf('\n LOCI COUNT FOR FEATURE MATRIX: %.0f \n\n',size(VLOCI,1))

clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL









%########################################################################
%% PLOT SCATTER GRAPH COMPARISONS OF TRAINING & TESTING GROUPS
%########################################################################
%{


nlogTRFISHP = -log(VLOCI.TRFISHP(1:end))+1;     
nlogTEFISHP = -log(VLOCI.TEFISHP(1:end))+1;     

nlogPPTRFISHP = -log(VLOCI.PPTRFISHP(1:end))+1;   
nlogPPTEFISHP = -log(VLOCI.PPTEFISHP(1:end))+1;   

nlogTRCHIP = -log(VLOCI.TRCHIP(1:end))+1;       
nlogTECHIP = -log(VLOCI.TECHIP(1:end))+1;       

nlogPPTRCHIP = -log(VLOCI.PPTRCHIP(1:end))+1;   
nlogPPTECHIP = -log(VLOCI.PPTECHIP(1:end))+1;   


TRCHIPinf = isinf(nlogTRCHIP);
TECHIPinf = isinf(nlogTECHIP);
TRCHIPnan = isnan(nlogTRCHIP);
TECHIPnan = isnan(nlogTECHIP);
nlogTRCHIP(TRCHIPinf) = max(nlogTRCHIP(~TRCHIPinf))+1;
nlogTECHIP(TECHIPinf) = max(nlogTECHIP(~TECHIPinf))+1;
nlogTRCHIP(TRCHIPnan) = 1;
nlogTECHIP(TECHIPnan) = 1;
nlogTRCHIP(nlogTRCHIP<1) = 1;
nlogTECHIP(nlogTECHIP<1) = 1;

PPTRCHIPinf = isinf(nlogPPTRCHIP);
PPTECHIPinf = isinf(nlogPPTECHIP);
PPTRCHIPnan = isnan(nlogPPTRCHIP);
PPTECHIPnan = isnan(nlogPPTECHIP);
nlogPPTRCHIP(PPTRCHIPinf) = max(nlogPPTRCHIP(~PPTRCHIPinf))+1;
nlogPPTECHIP(PPTECHIPinf) = max(nlogPPTECHIP(~PPTECHIPinf))+1;
nlogPPTRCHIP(PPTRCHIPnan) = 1;
nlogPPTECHIP(PPTECHIPnan) = 1;
nlogPPTRCHIP(nlogPPTRCHIP<1) = 1;
nlogPPTECHIP(nlogPPTECHIP<1) = 1;


disp('------------');disp(' ')
disp(min(nlogTRFISHP));   disp(max(nlogTRFISHP));
disp(min(nlogTEFISHP));   disp(max(nlogTEFISHP));
disp(' ')
disp(min(nlogPPTRFISHP)); disp(max(nlogPPTRFISHP));
disp(min(nlogPPTEFISHP)); disp(max(nlogPPTEFISHP));
disp(' ')
disp(min(nlogTRCHIP));    disp(max(nlogTRCHIP));
disp(min(nlogTECHIP));    disp(max(nlogTECHIP));
disp(' ')
disp(min(nlogPPTRCHIP));  disp(max(nlogPPTRCHIP));
disp(min(nlogPPTECHIP));  disp(max(nlogPPTECHIP));
disp('------------');disp(' ')




close all
fh01 = figure('Units','normalized','OuterPosition',[.01 .04 .90 .95],...
              'Color','w','MenuBar','none');
ax01 = axes('Position',[.06 .56 .4 .4],'Color','none');
ax02 = axes('Position',[.56 .56 .4 .4],'Color','none');
ax03 = axes('Position',[.06 .06 .4 .4],'Color','none');
ax04 = axes('Position',[.56 .06 .4 .4],'Color','none');


axes(ax01);
ph1 = loglog(nlogTRFISHP , nlogTEFISHP,'-mo',...
    'LineStyle','none',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','none',...
    'MarkerSize',6);
title('TRAINING VS TESTING FISHERS P-VALUE (PER CHROMOSOME)')


axes(ax02);
ph2 = loglog(nlogPPTRFISHP , nlogPPTEFISHP,'-mo',...
    'LineStyle','none',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','none',...
    'MarkerSize',6);
title('TRAINING VS TESTING FISHERS P-VALUE (PER PERSON)')


axes(ax03);
ph3 = loglog(nlogTRCHIP , nlogTECHIP,'-mo',...
    'LineStyle','none',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','none',...
    'MarkerSize',6);
title('TRAINING VS TESTING CHISQ P-VALUE (PER CHROMOSOME)')


axes(ax04);
ph4 = loglog(nlogPPTRCHIP , nlogPPTECHIP,'-mo',...
    'LineStyle','none',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','none',...
    'MarkerSize',6);
title('TRAINING VS TESTING CHISQ P-VALUE (PER PERSON)')





pause(.2)

% XYMAX = 550;
XYMAX = max([ax01.XLim ax02.XLim ax03.XLim ax04.XLim ...
             ax01.YLim ax02.YLim ax03.YLim ax04.YLim]);

axes(ax01); ax01.XLim = [1 XYMAX]; grid on;
axes(ax01); ax01.YLim = [1 XYMAX]; grid on;
axes(ax02); ax02.XLim = [1 XYMAX]; grid on;
axes(ax02); ax02.YLim = [1 XYMAX]; grid on;
axes(ax03); ax03.XLim = [1 XYMAX]; grid on;
axes(ax03); ax03.YLim = [1 XYMAX]; grid on;
axes(ax04); ax04.XLim = [1 XYMAX]; grid on;
axes(ax04); ax04.YLim = [1 XYMAX]; grid on;

pause(.2)



[FRHO,FPVAL] = corr(nlogTRFISHP,nlogTEFISHP);
[PPFRHO,PPFPVAL] = corr(nlogPPTRFISHP,nlogPPTEFISHP);
[CRHO,CPVAL] = corr(nlogTRCHIP,nlogTECHIP);
[PPCRHO,PPCPVAL] = corr(nlogPPTRCHIP,nlogPPTECHIP);


fprintf('\n FISHERS P CORRELATION (CHR): %.2f (P ~ %.1g) \n',FRHO,FPVAL)
fprintf('\n FISHERS P CORRELATION (PER): %.2f (P ~ %.1g) \n',PPFRHO,PPFPVAL)
fprintf('\n CHISQ   P CORRELATION (CHR): %.2f (P ~ %.1g) \n',CRHO,CPVAL)
fprintf('\n CHISQ   P CORRELATION (PER): %.2f (P ~ %.1g) \n',PPCRHO,PPCPVAL)



%------------------------------------------%
% set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf, ['/Users/bradleymonk/Desktop/CHIP_PER.png']);
%------------------------------------------%

clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL

%}



%% REMOVE TOP 100 VARIANTS (SPECIFIC TEST)
%{
% VLOCI(1:100,:)   = [];
% VCASE(1:100) = [];
% VCTRL(1:100) = [];
% VUSNP(1:100) = [];
% VLOCI.VID        = (1:size(VLOCI,1))';
%}







%###################################################################
%%      MAKE  RECTANGLE  NEURAL NET  VARIANT MATRIX
%###################################################################


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




% MAKE THE NEURAL NET TRAINING & TESTING MATRICES
[VTRX, TRX, TRL] = makeomicsnet(VLOCI,VCASE,VCTRL,VUSNP,TRPHE,[-1 -0 1 3]);
[VTEX, TEX, TEL] = makeomicsnet(VLOCI,VCASE,VCTRL,VUSNP,TEPHE,[-1 -0 1 3]);

[NTRX, NTRXT] = makenucleotidemx(VLOCI,VCASE,VCTRL,VUSNP,TRPHE,[-1 -0 1 3]);
[NTEX, NTEXT] = makenucleotidemx(VLOCI,VCASE,VCTRL,VUSNP,TEPHE,[-1 -0 1 3]);


% MAKE CASE=1 CTRL=-1
VTRX(:,2) = double(VTRX(:,2)) .*2 - 1;
VTEX(:,2) = double(VTEX(:,2)) .*2 - 1;
NTRX(:,2) = double(NTRX(:,2)) .*2 - 1;
NTEX(:,2) = double(NTEX(:,2)) .*2 - 1;


clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE...
VTRX VTEX TRX TRL TEX TEL NTRX NTRXT NTEX NTEXT




%% POST-PROCESS TRAINING & TESTING MATRICES


%######################################################################
%% POST-PROCESS TRAINING & TESTING MATRICES
%######################################################################
%
%-------------------------------------------------------------------------
% USING:   
%         REFREF = -1;    UNK = -1 or 0;    ALTALT = 1;    ALTALT = 3
% 
%--------------------------------------------------------------------------
% VTRX(: , 1)  =  TRPHE.SRR;        % COL1: ID
% VTRX(: , 2)  =  TRPHE.AD;         % COL2: AD
% VTRX(: , 3)  =  TRPHE.COHORTNUM;  % COL3: COHORT
% VTRX(: , 4)  =  TRPHE.AGE;        % COL4: AGE
% VTRX(: , 5)  =  TRPHE.APOE;       % COL5: APOE
% VTRX(: , 6)  =  TRPHE.BRAAK;      % COL6: BRAAK
%--------------------------------------------------------------------------


% DONT ADD AN INTERCEPT COLUMN FOR PCA ANALYSIS
VTRX = [VTRX(:,1:6) ones(size(VTRX,1),1) VTRX(:,7:end)];
VTEX = [VTEX(:,1:6) ones(size(VTEX,1),1) VTEX(:,7:end)];
NTRX = [NTRX(:,1:6) ones(size(NTRX,1),1) NTRX(:,7:end)];
NTEX = [NTEX(:,1:6) ones(size(NTEX,1),1) NTEX(:,7:end)];


VTRVXP = VTRX(:,7:end)';    % TRAINING VARIANT MATRIX - VARIANT X PERSON
VTRPXV = VTRVXP';           % TRAINING VARIANT MATRIX - PERSON X VARIANT

VTEVXP = VTEX(:,7:end)';    % TESTING VARIANT MATRIX - VARIANT X PERSON
VTEPXV = VTEVXP';           % TESTING VARIANT MATRIX - PERSON X VARIANT


NTRVXP = NTRX(:,7:end)';    % TRAINING NUCLEOTIDE MATRIX - VARIANT X PERSON
NTRPXV = NTRVXP';           % TRAINING NUCLEOTIDE MATRIX - PERSON X VARIANT

NTEVXP = NTEX(:,7:end)';    % TESTING NUCLEOTIDE MATRIX - VARIANT X PERSON
NTEPXV = NTEVXP';           % TESTING NUCLEOTIDE MATRIX - PERSON X VARIANT


fprintf('\n VTRPXV RANK %0.f; VTRPXV DIMS %0.f-by-%0.f \n',rank(VTRPXV),size(VTRPXV));
fprintf('\n NTRPXV RANK %0.f; NTRPXV DIMS %0.f-by-%0.f \n',rank(NTRPXV),size(NTRPXV));


% TRAINING VARIANT MATRIX CLASS LABELS
VTRL = zeros( numel(VTRX(: , 2)) ,2);  % Make 2-D class label matrix
VTRL(:,1) = VTRX(: , 2)==1;              % In col-1 set CASEs index 1's
VTRL(:,2) = VTRX(: , 2)~=1;              % In col-2 set CTRLs index 1's

% TESTING VARIANT MATRIX CLASS LABELS
VTEL = zeros( numel(VTEX(: , 2)) ,2);  % Make 2-D class label matrix
VTEL(:,1) = VTEX(: , 2)==1;              % In col-1 set CASEs index 1's
VTEL(:,2) = VTEX(: , 2)~=1;              % In col-2 set CTRLs index 1's


% TRAINING NUCLEOTIDE MATRIX CLASS LABELS
NTRL = zeros( numel(NTRX(: , 2)) ,2);  % Make 2-D class label matrix
NTRL(:,1) = NTRX(: , 2)==1;              % In col-1 set CASEs index 1's
NTRL(:,2) = NTRX(: , 2)~=1;              % In col-2 set CTRLs index 1's

% TESTING NUCLEOTIDE MATRIX CLASS LABELS
NTEL = zeros( numel(NTEX(: , 2)) ,2);  % Make 2-D class label matrix
NTEL(:,1) = NTEX(: , 2)==1;              % In col-1 set CASEs index 1's
NTEL(:,2) = NTEX(: , 2)~=1;              % In col-2 set CTRLs index 1's





clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE...
VTRX VTEX NTRX NTEX VTRVXP VTRPXV VTEVXP VTEPXV NTRVXP NTRPXV NTEVXP NTEPXV...
VTRL VTEL NTRL NTEL


%--------------------------------------------------------------------------
% VTRX(: , 1)  =  TRPHE.SRR;        % COL1: ID
% VTRX(: , 2)  =  TRPHE.AD;         % COL2: AD
% VTRX(: , 3)  =  TRPHE.COHORTNUM;  % COL3: COHORT
% VTRX(: , 4)  =  TRPHE.AGE;        % COL4: AGE
% VTRX(: , 5)  =  TRPHE.APOE;       % COL5: APOE
% VTRX(: , 6)  =  TRPHE.BRAAK;      % COL6: BRAAK
% VTRX(: , 7)  =  INTERCEPT
% VTRX(: , 8:end)  =  VARIANT CODE
%
% VTRL(:,1) = CASE LABELS
% VTRL(:,2) = CTRL LABELS
% 
% NTRX(: , 1)  =  TRPHE.SRR;        % COL1: ID
% NTRX(: , 2)  =  TRPHE.AD;         % COL2: AD
% NTRX(: , 3)  =  TRPHE.COHORTNUM;  % COL3: COHORT
% NTRX(: , 4)  =  TRPHE.AGE;        % COL4: AGE
% NTRX(: , 5)  =  TRPHE.APOE;       % COL5: APOE
% NTRX(: , 6)  =  TRPHE.BRAAK;      % COL6: BRAAK
% NTRX(: , 7)  =  INTERCEPT
% NTRX(: , 8:end)  =  NUCLEOTIDE CODE
%
% NTRL(:,1) = CASE LABELS
% NTRL(:,2) = CTRL LABELS
%--------------------------------------------------------------------------




% deTRX = decomposition(TRX,'cod');
% deTEX = decomposition(TEX,'cod');
% rcond(deTRX)
% rcond(deTEX)
% lscovTR = lscov(TRX',TRL(1,:)');
% guessTE = (lscovTR'*(TEX)) > 0;
% mean(guessTE == TEL(1,:))

% [V,D,W] = eig(TRX)


% TRXAGE = FULLTRX(2:end,3)';
% TEXAGE = FULLTEX(2:end,3)';
% TRX = [TRX; TRXAGE];
% TEX = [TEX; TEXAGE];




%%
disp('READY TO LEARN!'); 
return
%##########################################################################
%##########################################################################
%##########################################################################
%
%%                          MACHINE LEARNING
%
%##########################################################################
%##########################################################################
%##########################################################################


% save('nnvars.mat','VLOCI','VCASE','VCTRL','VUSNP','VTRCASE','VTRCTRL',...
% 'VTECASE','VTECTRL','TRPHE','TEPHE','FULLTRX','FULLTEX','TRX','TRL','TEX','TEL')
% 
% load('nnvars.mat');





%##########################################################################
%
%%          MACHINE LEARNING USING LOGISTIC REGRESSION
%
%##########################################################################
clc; clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE...
FULLTRX FULLTEX TRX TRL TEX TEL TrHAPLON TrHAPLOT TeHAPLON TeHAPLOT



% TRX & TEX ARE VARIANT-BY-PERSON; WE NEED PERSON-BY-VARIANT
TrPV = TRX';  
TrL  = TRL';
TePV = TEX';
TeL  = TEL';



% Create labels CASE=1 CTRL=0
TR_y = TrL(:,1);
TE_y = TeL(:,1);


%-------------------------------------------------
% USE THE TOP Vi VARIANTS FOR INTERACTION TERMS
doINTERACTIONS = 0

if doINTERACTIONS

    Vi = 25;
    TrPVi = TrPV(:,1:Vi);
    TePVi = TePV(:,1:Vi);



    % DETERMINE HOW MANY MORE ROWS THIS WILL ADD
    cs = cumsum(fliplr(1:(Vi-1)));
    added = cs(end); 
    fprintf('\n Including all two-way interactions for ');
    fprintf('\n   top *%0.f* variants will add *%0.f* rows. \n',Vi,added);
    ci = [0 cs]; 
    pause(1);




    % PREALLOCATE
    TrPVii = zeros(size(TrPV,1),added);
    TePVii = zeros(size(TePV,1),added);



    % GENERATE THE INTERACTION TERMS
    szx = 0;
    for n = 1:(size(TrPVi,2)-1)

        x = (TrPVi(:,n) .* TrPVi(:,(n+1):end));

        TrPVii(:,(ci(n)+1):(ci(n+1))) = (TrPVi(:,n) .* TrPVi(:,(n+1):end));
        TePVii(:,(ci(n)+1):(ci(n+1))) = (TePVi(:,n) .* TePVi(:,(n+1):end));


    szx = szx + size(x,2);
    end
    fprintf('\n Created %0.f (rows of) variant interaction terms. \n',szx)
    pause(1);



    % CONCATINATE ORIGINAL PV MATRIX WITH INTERACTION TERMS
    TrPV = [TrPV TrPVii];
    TePV = [TePV TePVii];


    % CONCATINATE PV MATRIX WITH QUADRADIC INTERACTION TERM
    TrPV = [TrPV TrPV.^2];
    TePV = [TePV TePV.^2];

end
%-------------------------------------------------



% ADD AN INTERCEPT COLUMN IF NEEDED
if mean(TrPV(:,1)) == 1
    % intercept was already added above
    TR_x = TrPV;
    TE_x = TePV;
else
    % intercept was not added above
    TR_x = [ones(size(TrPV,1),1) TrPV];
    TE_x = [ones(size(TePV,1),1) TePV];
end




%--------- PERFORM THE SO-CALLED MACHINE LEARNING STEP ------------------
% USE THE **NORMAL EQUATION** TO MINIMIZE ERROR (TRAINING SET ONLY)

BETA = pinv(TR_x' * TR_x) * (TR_x' * TR_y);

fprintf('\n Solved normal equation OLS for %0.f beta coefficients. \n\n',...
    size(BETA,1));  pause(.5)



% Get linear model predictions for each y:
TR_fx = sum( TR_x .* BETA' ,2);
TE_fx = sum( TE_x .* BETA' ,2);


% Descretize predictions
TR_yp = round(TR_fx);
TE_yp = round(TE_fx);


% GRADE THE PREDICTIONS
TRmu = mean(TR_yp == TR_y);
TEmu = mean(TE_yp == TE_y);




% Eval performance on high confidence outputs
TRhi = (TR_fx>.8) | (TR_fx<.2);
TRhiN  = TR_yp(TRhi) == TR_y(TRhi);
TRhiMu = mean(TRhiN);


TEhi = (TE_fx>.8) | (TE_fx<.2);
TEhiN  = TE_yp(TEhi) == TE_y(TEhi);
TEhiMu = mean(TEhiN);




disp('----------  TRAINING SET  ----------')
fprintf('TRAIN correct %0.0f%% overall(pop:100%%)\n' ,(TRmu .* 100))

fprintf('TRAIN correct %0.0f%% hicon  (pop: %0.0f%%)\n\n',TRhiMu.*100,mean(TRhi).*100)

disp('----------  TESTING SET  ----------')
fprintf('TEST correct %0.0f%% overall(pop:100%%)\n' ,(TEmu .* 100))

fprintf('TEST correct %0.0f%% hicon  (pop: %0.0f%%)\n\n',TEhiMu.*100,mean(TEhi).*100)



% SAVE MATRICIES THAT INCLUDE THE INTERACTION TERMS
TRXI = TR_x';
TEXI = TE_x';

clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE...
FULLTRX FULLTEX TRX TRL TEX TEL TrHAPLON TrHAPLOT TeHAPLON TeHAPLOT...
BETA TRXI TEXI 





%################################################################
%%  ASSIGN MACHINE LEARNING ACTIVATION VALUE TO PARTICIPANT PHE
%################################################################
%{
numnets = size(netz,2);
ACTR = zeros(size(TRX,2),numnets);
ACTE = zeros(size(TEX,2),numnets);

for nn = 1:numnets
    net = netz{nn};
    ACTI = net(TRX);
    ACTR(:,nn) = (ACTI(1,:))';
    ACTI = net(TEX);
    ACTE(:,nn) = (ACTI(1,:))';
end



BRAAK = TRPHE.BRAAK;

BRAAK(TRPHE.AUTOPSY==0) = NaN;


clc;

nanmean(BRAAK(TRPHE.MUACTR > .90))

nanmean(BRAAK( ((TRPHE.MUACTR<.90)&(TRPHE.MUACTR>.80))  ))

nanmean(BRAAK( ((TRPHE.MUACTR<.80)&(TRPHE.MUACTR>.70))  ))

nanmean(BRAAK( ((TRPHE.MUACTR<.70)&(TRPHE.MUACTR>.60))  ))

nanmean(BRAAK( ((TRPHE.MUACTR<.60)&(TRPHE.MUACTR>.50))  ))

nanmean(BRAAK(TRPHE.MUACTR < .5))
%}








%##########################################################################
%
%%                ARTIFICIAL NEURAL NETWORKS (MATLAB BUILT-IN)
%
%##########################################################################
rng('shuffle');

clc;clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE...
FULLTRX FULLTEX TRX TRL TEX TEL BETA TRXI TEXI



%-------------------------------------------------------
%   TRX & TEX CURRENTLY:      VARIANT-BY-PERSON 
%   MATLAB BUILTIN NN NEEDS:  VARIANT-BY-PERSON 
%    (*NO TRANSPOSE NEEDED*)
%-------------------------------------------------------


%--- DONT INCLUDE INTERACTION TERMS
% TRAINX  =  TRXI;  
% TRAINL  =  TRL;   
% TESTX   =  TEXI;  
% TESTL   =  TEL;


%--- DO INCLUDE INTERACTION TERMS
TRAINX  =  TRXI;
TRAINL  =  TRL;
TESTX   =  TEXI;
TESTL   =  TEL;



fprintf('\n size(TRAINX): %0.f x %0.f \n', size(TRAINX))
fprintf('\n size(TRAINL): %0.f x %0.f \n', size(TRAINL))
fprintf('\n size(TESTX):  %0.f x %0.f \n', size(TESTX))
fprintf('\n size(TESTL):  %0.f x %0.f \n', size(TESTL))




%--------------------------------------------------------
%---       TRAINING FUNCTION OPTION DEFINITIONS       ---
%--------------------------------------------------------
% NN.trainFcn = 'traingd';   %(GRADESC)
% NN.trainFcn = 'traingdm';  %(GRADESC + MOMENTUM)
% NN.trainFcn = 'traingda';  %(GRADESC + ADAPTIVE)
% NN.trainFcn = 'traingdx';  %(GRADESC + ADAPTIVE + MOMENTUM)


NN = patternnet([200 100 100]);

NN.trainFcn = 'traingdx';

NN.trainParam.max_fail = 20;

disp(NN)

% NN.layers{1}
% weights = getwb(NN)
% size(weights)
% 200 * 100 * 100
%  weights = NN.LW
% biases = NN.b
% view(NN)



%--------------------------------------------------------
%---   TRAIN THE ARTIFICIAL NEURAL NETWORK WEIGHTS    ---
%--------------------------------------------------------
net = train(NN,TRAINX,TRAINL);



%--- EVALUATE TRAINED NEURAL NETWORK PERFORMANCE ---
% [NNPerf] = netstats(net,TrVP,TrVPL,TeVP,TeVPL,.05,.85,.15,1);
[NNPerf] = nnetstats(net,TRAINX,TRAINL,TESTX,TESTL,.05,.85,.15,1);



%------------------------------------------%
set(gcf, 'PaperPositionMode', 'auto');
dt=char(datetime(datetime,'Format','yyyy-MM-dd-HH-mm-ss'));
saveas(gcf, ['/Users/bradleymonk/Desktop/FIG_' dt '.png']);
%------------------------------------------%










%##########################################################################
%
%%  MACHINE LEARNING USING HOMEBREW NEURAL NETS & MINIMIZER (FMINCG)
%
%##########################################################################
% LOGISTIC REGRESSION NEURAL NET WITH GRADIENT DESCENT


clc;clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE...
FULLTRX FULLTEX TRX TRL TEX TEL BETA TRXI TEXI





%-------------------------------------------------------
%   TRX & TEX CURRENTLY:     VARIANT-BY-PERSON 
%   THIS HOMEBREW NN NEEDS:  PERSON-BY-VARIANT 
%    (*TRANSPOSE NEEDED*)
%-------------------------------------------------------


%--- DONT INCLUDE INTERACTION TERMS
% TRAINX  =  TRX';
% TRAINL  =  (TRL(1,:)')+1;
% TESTX   =  TEX';
% TESTL   =  (TEL(1,:)')+1;



%--- DO INCLUDE INTERACTION TERMS
TRAINX  =  TRXI';
TRAINL  =  (TRL(1,:)')+1;
TESTX   =  TEXI';
TESTL   =  (TEL(1,:)')+1;



fprintf('\n size(TRAINX): %0.f x %0.f \n', size(TRAINX))
fprintf('\n size(TRAINL): %0.f x %0.f \n', size(TRAINL))
fprintf('\n size(TESTX):  %0.f x %0.f \n', size(TESTX))
fprintf('\n size(TESTL):  %0.f x %0.f \n', size(TESTL))



% RANDOMIZE ROWS
randp     = randperm(size(TRAINX,1));
TRAINX    = TRAINX(randp,:);
TRAINL    = TRAINL(randp,:);


% SET NUMBER OF L1 NEURONS EQUAL TO NUMBER OF VARIANTS
LAY1N  = size(TRAINX,2);
nLabels = 2;



% NEURAL NET PARAMETERS
%----------------------------------------------------------------------
lambda    = 0.005;       % .001 - .01 is a good operational window
epsInit   = 0.22;        % random initial theta weights
maxIters  = 100;         % 20-50 iterations should be sufficient
LAY2N     = 100;         % 10 - 50 neurons should be sufficient
LAY3N     = 50;          % 10 - 50 neurons should be sufficient
%----------------------------------------------------------------------

%function [J, grad] = omicscostfun3L(Tin, L1n, L2n, L3n, nLabs, X, y, lambda, AF)








% ESTABLISH NEURON ACTIVATION FUNCTION (choose 1)
%----------------------------------------------------------------------
% ACTIVATION FUNCTION LIST: wikipedia.org/wiki/Activation_function

lintrans = @(x,a,b,c,d) (c.*(1-(x-a)./(b-a)) + d.*((x-a)./(b-a)));
f  = @(z) 1 ./ (1 + exp(-z) );       % sigmoid activation function
f1 = @(z) 1 ./ (1 + exp(-z) );       % sigmoid activation function
f2 = @(g) f(g) .* (1-f(g));          % sigmoid gradient function
f3 = @(n) 2 ./ (1 + exp(-2*n)) - 1;  % sigmoid symetric transfer func
f4 = @(x) tanh(x);                   % TanH
f5 = @(x) atan(x);                   % ArcTan
f6 = @(x) x ./ sqrt(1+.01.*x.^2);    % ISRU

AF = f1;





%----------------------------------------------------------------------
% PREP NEURAL NET SYNAPTIC WEIGHTS (THETAS) & GRADIENT DESCENT COST FUNC
%----------------------------------------------------------------------




% INITIALIZE RANDOM THETA WEIGHTS
%----------------------------------------------------------------------
initTheta1 = rand(LAY2N, LAY1N+1) * 2 * epsInit - epsInit;
initTheta2 = rand(LAY3N, LAY2N+1) * 2 * epsInit - epsInit;
initTheta3 = rand(nLabels, LAY3N+1) * 2 * epsInit - epsInit;

Th1=initTheta1;
Th2=initTheta2;
Th3=initTheta3;





% UNROLL THETA MATRICES
%----------------------------------------------------------------------

T = [Th1(:) ; Th2(:); Th3(:)];




% PRACTICE REROLLING THETA MATRICES
%----------------------------------------------------------------------
nT1 = LAY2N .* (LAY1N+1);   %40100
nT2 = LAY3N .* (LAY2N+1);   %5050
nT3 = nLabels .* (LAY3N+1); %102
Theta1 = T(1:nT1);
Theta2 = T((nT1+1):(nT1+nT2));
Theta3 = T((nT1+nT2+1):(nT1+nT2+nT3));
T1 = reshape(Theta1, LAY2N, LAY1N+1);
T2 = reshape(Theta2, LAY3N, LAY2N+1);
T3 = reshape(Theta3, nLabels, LAY3N+1);



% ESTABLISH COST FUNCTION
%----------------------------------------------------------------------
% J = omicscostfun(nn_Thetas, LAY1N, LAY2N, nLabels, TRAINX, TRAINL, lambda, AF);
J = omicscostfun3L(T, LAY1N, LAY2N, LAY3N, nLabels, TRAINX, TRAINL, lambda, AF);





% TRAIN NEURAL NETWORK USING TRAINING DATASET
%----------------------------------------------------------------------
options = optimset('MaxIter', maxIters);

% GcostFun = @(T) omicscostfun(T, LAY1N, LAY2N, nLabels, TRAINX, TRAINL, lambda, AF);
GcostFun = @(T) omicscostfun3L(T, LAY1N, LAY2N, LAY3N, nLabels, TRAINX, TRAINL, lambda, AF);

% [nn_Thetas, cost] = omicsfmincg(GcostFun, initial_Thetas, options);
[nn_Thetas, cost] = omicsfmincg3L(GcostFun, T, options);

% [nn_Thetas, cost] = fmincon(GcostFun,T)
%----------------------------------------------------------------------



%   options = optimset('MaxIter', 500 , 'Display', 'iter', 'MaxFunEvals', 1000);
%   objFunc = @(t) lrCostFunction(t,X,y);
%   [result1] = fminsearch(objFunc, theta, options);
%   [result2]=  fmincg (objFunc, theta, options);





% REROLL THETA WEIGHTS
%----------------------------------------------------------------------
T = nn_Thetas;

nT1 = LAY2N .* (LAY1N+1);   %40100
nT2 = LAY3N .* (LAY2N+1);   %5050
nT3 = nLabels .* (LAY3N+1); %102
Theta1 = T(1:nT1);
Theta2 = T((nT1+1):(nT1+nT2));
Theta3 = T((nT1+nT2+1):(nT1+nT2+nT3));
T1 = reshape(Theta1, LAY2N, LAY1N+1);
T2 = reshape(Theta2, LAY3N, LAY2N+1);
T3 = reshape(Theta3, nLabels, LAY3N+1);


% close all;
% subplot(3,1,1); bar(Theta1)
% subplot(3,1,2); bar(Theta2)
% subplot(3,1,3); bar(Theta3)
% h1 = AF([ones(size(TESTX(1,:),1), 1) TESTX(1,:)] * T1');
% h2 = AF([ones(size(TESTX(1,:),1), 1) h1] * T2');
% h3 = AF([ones(size(TESTX(1,:),1), 1) h2] * T3');
% [a, p] = max(h3, [], 2);
% close all;
% subplot(3,1,1); bar(h1)
% subplot(3,1,2); bar(h2)
% subplot(3,1,3); bar(h3)




% EVALUATE PERFORMANCE ON **TRAINING** SET
%----------------------------------------------------------------------


h1 = AF([ones(size(TRAINX,1), 1) TRAINX] * T1');
h2 = AF([ones(size(TRAINX,1), 1) h1] * T2');
h3 = AF([ones(size(TRAINX,1), 1) h2] * T3');
[a, p] = max(h3, [], 2);

TRAINPCTCORRECT = mean(p == TRAINL);

disp('Percent accuracy on training data:')
disp( TRAINPCTCORRECT )




% EVALUATE NEURAL NETWORK PERFORMANCE ON **TEST** SET
%----------------------------------------------------------------------

h1 = AF([ones(size(TESTX,1), 1) TESTX] * T1');
h2 = AF([ones(size(TESTX,1), 1) h1] * T2');
h3 = AF([ones(size(TESTX,1), 1) h2] * T3');
[a, p] = max(h3, [], 2);

TESTPCTCORRECT = mean(p == TESTL);

disp('Percent accuracy on testing data:')
disp( TESTPCTCORRECT )
























% clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
% VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE...
% FULLTRX FULLTEX TRX TRL TEX TEL BETA TRXI TEXI



%##########################################################################
%
%%     MACHINE LEARNING WITH HOMEBREW HACKY VANILLA NEURAL NETS
%                  INCLUDES LIVE SIMULATION RESULTS
%
%##########################################################################

clc;clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE...
FULLTRX FULLTEX TRX TRL TEX TEL BETA TRXI TEXI





%-------------------------------------------------------
%   TRX & TEX CURRENTLY:     VARIANT-BY-PERSON 
%   THIS HOMEBREW NN NEEDS:  PERSON-BY-VARIANT 
%    (*TRANSPOSE NEEDED*)
%-------------------------------------------------------


%--- DONT INCLUDE INTERACTION TERMS
TRAINX  =  TRX';
TRAINL  =  TRL(1,:)';
TESTX   =  TEX';
TESTL   =  TEL(1,:)';



%--- DO INCLUDE INTERACTION TERMS
% TRAINX  =  TRXI';
% TRAINL  =  TRL(1,:)';
% TESTX   =  TEXI';
% TESTL   =  TEL(1,:)';



fprintf('\n size(TRAINX): %0.f x %0.f \n', size(TRAINX))
fprintf('\n size(TRAINL): %0.f x %0.f \n', size(TRAINL))
fprintf('\n size(TESTX):  %0.f x %0.f \n', size(TESTX))
fprintf('\n size(TESTL):  %0.f x %0.f \n', size(TESTL))




%--- PREP DATA FOR NN CLASSIFIER TRAINING
%---------------------------------------------

INPUTMXdims   = size(TRAINX);
INPUTYNdims   = size(TRAINL);

randp        = randperm(INPUTMXdims(1));
TRAINX       = TRAINX(randp,:);
TRAINL       = TRAINL(randp,:);


% PREVIEW LABEL & SAMPLE MATRICES
T=table(TRAINL(1:9), TRAINX(1:9,1:4));
T.Properties.VariableNames = {'LABELS','VARIANTS'};
head(T)



%-- CUSTOMIZABLE PARAMETERS
%---------------------------------------------
weightScaleFactor = 0.02;
L1NEURONS = 200;
L2NEURONS = 100;
pctThresh = 98;
nEpochs = 70;



%-- MAKE THREE-LAYER NEURAL NETWORK
%---------------------------------------------

L1WEIGHTS = rand( L1NEURONS , INPUTMXdims(2) ) * 2 - 1;
L2WEIGHTS = rand( L2NEURONS , L1NEURONS ) * 2 - 1;
L3WEIGHTS = rand( INPUTYNdims(2) , L2NEURONS  ) * 2 - 1;



%-- CREATE NEURON ACTIVATION FUNCTIONS
%---------------------------------------------
% ACTIVATION FUNCTION LIST:
%   https://en.wikipedia.org/wiki/Activation_function
%---------------------------------------------

lintrans = @(x,a,b,c,d) (c.*(1-(x-a)./(b-a)) + d.*((x-a)./(b-a)));
f0 = @(x) (1./(1+exp(-x)));          % sigmoid activation function
f1 = @(z) 1 ./ (1 + exp(-z) );       % sigmoid activation function
f2 = @(g) f(g) .* (1-f(g));          % sigmoid gradient function
f3 = @(n) 2 ./ (1 + exp(-2*n)) - 1;  % sigmoid symetric transfer func
f4 = @(x) tanh(x);                   % TanH
f5 = @(x) atan(x);                   % ArcTan
f6 = @(x) x ./ sqrt(1+.01.*x.^2);    % ISRU



%-- CHOOSE NEURON ACTIVATION FUNCTION
%---------------------------------------------
f=f0;

fLO=-1; 
fHI=2;




%-- PERFORM SINGLE INSTANCE OF NN TRAINING TO PREPOPULATE PLOTS
%---------------------------------------------

% Get first training example
In = TRAINX(1,:)';
YN = TRAINL(1);


% Propagate the signals through network
L1OUTS = f( L1WEIGHTS * In  );
L2OUTS = f( L2WEIGHTS * L1OUTS);
L3OUTS = f( L3WEIGHTS * L2OUTS);
% lintrans(L3OUTS,-10,10,0,1)





%-- SETUP FIGURES, AXES, & PLOTS
%---------------------------------------------
close all; clc
fh1  = figure('Units','normalized','OuterPosition',[.2 .1 .75 .88],'Color','w');


% GROUP SEPARATION GRAPH
hax1 = axes('Position',[.05 .35 .75 .25],'Color','none');
    hax1.XLim = [-.5 .5]; hax1.XLimMode = 'manual';
    colormap(hax1,cool)
    hold on;



% WEIGHTS FIRST LAYER
hax2 = axes('Position',[.05 .82 .9 .14],'Color','none');
    hax2.YLim = [fLO fHI]; hax2.YLimMode = 'manual';
    hold on;


% WEIGHTS 2ND LAYER
hax3 = axes('Position',[.05 .65 .9 .14],'Color','none');
    hax3.YLim = [fLO fHI]; hax3.YLimMode = 'manual';
    hold on;



% OVERALL PERCENT CORRECT
hax4 = axes('Position',[.85 .05 .1 .5],'Color','none');
    hax4.YLim = [0 100]; hax4.YLimMode = 'manual';
    hold on;
    

% GRADIENT DESCENT
hax5 = axes('Position',[.05 .05 .75 .25],'Color','none');
%     hax5.YLim = [0 1]; hax5.YLimMode = 'manual';
%     hold on;



axes(hax1)
ph1 = scatter( (rand(size(TRAINX,1),1)-.5).*.1 ,...
            (1:size(TRAINX,1)), 500 , TRAINL, '.');

axes(hax2)
ph2 = bar(hax2, L1OUTS );

axes(hax3)
ph3 = bar(hax3, L2OUTS );

axes(hax4)
ph4 = bar( 50 ); % Percent Correct


RMSmu = zeros(nEpochs,1);
axes(hax5)
ph5 = scatter( 1:nEpochs , RMSmu , 1500 , '.r');

% t4 = text(.6,1.05,'GUESS: CASE!','Color','r','FontSize',16);    
pause(1.5)





%----------------------------------------------------------------------    
%----------------------------------------------------------------------
% START MAIN TRAINING LOOP FOR 3-LAYER SUPERVISED LEARNING
%----------------------------------------------------------------------
%----------------------------------------------------------------------


% figure(fh1)
% y = zeros(nEpochs,1);


for nn=1:nEpochs
    fprintf('\nEpoch: % .0f\n',nn)

    RMS_Err = 0;

    LayerB_YN = ones(size(TRAINL));
    L3O       = zeros(size(TRAINL));
    


    %-- Iterate through all examples
    for i=1:size(TRAINL,1)

        
        %-- Input data from current example set
        In = TRAINX(i,:)';
        YN = TRAINL(i);
                

        
        %-- Propagate the signals through network
        L1OUTS = f(L1WEIGHTS * In    );
        L2OUTS = f(L2WEIGHTS * L1OUTS);
        L3OUTS = f(L3WEIGHTS * L2OUTS);
        
        
        
        %-- GET NETWORK DELTAS AND RMS ERROR
        delta_i = L3OUTS .* (1-L3OUTS) .* (YN-L3OUTS);            % NODE ERRORS IN LayerB
        
        delta_k = L2OUTS .* (1-L2OUTS) .* (L3WEIGHTS' * delta_i); % NODE ERRORS IN LayerAB

        delta_j = L1OUTS .* (1-L1OUTS) .* (L2WEIGHTS' * delta_k); % NODE ERRORS IN LayerA
                
        RMS_Err = RMS_Err + norm(YN - L3OUTS);                    % ROOT MSE
        
        
        LayerB_YN(i) = L3OUTS;
        


        %-- TWO LAYER NETWORK
        % U = U + weightScaleFactor .* delta_i * LayerA' ;
        % W = W + weightScaleFactor .* delta_j * In' ;
        
        %-- THREE LAYER NETWORK
        L2WEIGHTS = L2WEIGHTS + weightScaleFactor .* delta_k * L1OUTS' ;
        L3WEIGHTS = L3WEIGHTS + weightScaleFactor .* delta_i * L2OUTS' ;
        L1WEIGHTS = L1WEIGHTS + weightScaleFactor .* delta_j * In' ;
        
        L3O(i) = L3OUTS;
    end
    
    
    numCorrect = sum(TRAINL == round(LayerB_YN));
    pctCorrect = numCorrect / INPUTYNdims(1) * 100;
    
    fprintf('\nCorrect Guesses: %.0f of %.0f  (%.2f pct.) \n',...
            numCorrect, size(TRAINL,1) , pctCorrect)
    
    if pctCorrect > pctThresh
        fprintf('\nPercent correct is greater than: (%.2f pct.) \n', pctThresh)
        disp('Moving on to test phase...')
        pause(2)
        break
    end
    
    
    RMSmu(nn) = RMS_Err;
    
    %if i==1
    ph1.XData = L3O-.5;
    ph2.YData = L1OUTS;
    ph3.YData = L2OUTS;
    ph4.YData = pctCorrect;
    ph5.YData = RMSmu;
    pause(.05)
    %end

end



disp('------------------------------------------------')
disp('      EVALUATING NEURAL NET PERFORMANCE...      ')
disp('------------------------------------------------')


%--- PREP DATA FOR NN CLASSIFIER EVALUATION
%---------------------------------------------

TESTX   =  TEX';
TESTL   =  TEL(1,:)';

randp       = randperm(size(TESTL,1));
TESTX       = TESTX(randp,:);
TESTL       = TESTL(randp,:);


%-- PREVIEW LABEL & SAMPLE MATRICES
T=table(TRAINL(1:9), TRAINX(1:9,1:4));
T.Properties.VariableNames = {'LABELS','VARIANTS'};
head(T)



%----------------------------------------------------------------------
% RUN TEST DATASET ONCE THROUGH THE NN (WITHOUT UPDATING WEIGHTS)
% TO GET NEURAL NET ACTIVATION OUTPUT VALUES
%----------------------------------------------------------------------

for nn=1:1

    RMS_Err = 0;
    
    LayerB_YN = ones(size(TESTL));
    
    for i=1:size(TESTL,1)
        
        In = TESTX(i,:)';
        YN = TESTL(i);
        
        L1OUTS = f(L1WEIGHTS*In);
        L2OUTS = f(L2WEIGHTS*L1OUTS);
        L3OUTS = f(L3WEIGHTS*L2OUTS);
        delta_i = L3OUTS .* (1-L3OUTS) .* (YN-L3OUTS);
        delta_k = L2OUTS .* (1-L2OUTS) .* (L3WEIGHTS' * delta_i);
        delta_j = L1OUTS .* (1-L1OUTS) .* (L2WEIGHTS' * delta_k);
        RMS_Err = RMS_Err + norm(YN - L3OUTS);
        LayerB_YN(i) = L3OUTS;
        
    end
    
    numCorrect = sum(TESTL == round(LayerB_YN));
    pctCorrect = numCorrect / size(TESTL,1) * 100;
    
    fprintf('\n\nCorrect Guesses: %.0f of %.0f  (%.1f%%) \n\n', ...
            numCorrect, size(TESTL,1) , pctCorrect)
    
end

fprintf('Number of CASES: % .0f\n',sum(TESTL))
fprintf('Number of CTRLS: % .0f\n',size(TESTL,1) - sum(TESTL))







% clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
% VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE...
% FULLTRX FULLTEX TRX TRL TEX TEL BETA TRXI TEXI



%##########################################################################
%
%%  TEST IDENTITY MATRIX TO DETERMINE WHAT SNPS IMPACT NN OUTPUT MOST
%
%##########################################################################
%{
ADNN(2:end,3) =  1; % BIAS COL#3:
TRAINX       = ADNN(2:end,3:end);
TRAINL       = ADNN(2:end,2);

INPUTMXdims   = size(TRAINX);
INPUTYNdims   = size(TRAINL);


IdentityMx = eye(INPUTMXdims(2));
IdentityMx(:,end) = 1;


LayerB_Response = zeros(INPUTMXdims(2),1);

f = @(x) (1./(1+exp(-x)));
 
for i=1:INPUTMXdims(2)
    In = IdentityMx(i,:)';
    L1OUTS = f(L1WEIGHTS*In);
    L2OUTS = f(L2WEIGHTS*L1OUTS);
    L3OUTS = f(L3WEIGHTS*L2OUTS);
    LayerB_Response(i) = L3OUTS;
end

NNresponseMx = LayerB_Response;

[NNcoRes,NNcoInd] = sort(NNresponseMx);
NNcaRes = flipud(NNcoRes);
NNcaInd = flipud(NNcoInd);

format shortG
[NNcaRes(1:10) NNcaInd(1:10) NNcoRes(1:10) NNcoInd(1:10)]

Ng = 20;

CTRLgenes = NNcoInd(1:Ng);
CTRLgeneWeights = NNcoRes(1:Ng);

CASEgenes = NNcaInd(1:Ng);
CASEgeneWeights = NNcaRes(1:Ng);

INPUTX = ADNN(:,3:end);
%}
%################################################################



