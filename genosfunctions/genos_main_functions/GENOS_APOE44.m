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

ID CONSORT	STUDY_NAME  STUDY   FAM PEOPLE  ENR     CASE	CTRL	TOT    %CA %CO  EQL	CULL	NOTES
01  DGC     Adult Chng  ACT     0	1277	0       323     945     1268	25	75	323	-622	01 GOOD
02  ADGC	AD Centers  ADC     0	3263	0       2438	817     3255	75	25	817	1621	02 GOOD
03  CHARGE	Athrosclro  ARIC	0	57      0       39      18      57      68	32	18	21      03 MEH
04  CHARGE	Aus Stroke	ASPS	0	126     0       121     5       126     96	4	5	116     04 BAD
05  ADGC	Chic Aging  CHAP	0	231     0       27      204     231     12	88	27	-177	05 OKAY
06  CHARGE	CardioHlth  CHS     0	840     0       250     583     833     30	70	250	-333	06 GOOD
07  ADGC	Hispanic S  CUHS	361	0       331     160     171     331     48	52	160	-11     07 BAD
08  CHARGE	Erasmus Er  ERF     30	45      0       45      0       45      100	0	0	45      08 BAD
09  CHARGE	Framingham  FHS     0	581     0       157     424     581     27	73	157	-267	09 GOOD
10	ADGC	Gene Diffs  GDF     0	207     0       111     96      207     54	46	96	15      10 GOOD
11	ADGC	NIA - LOAD  LOAD	80  113     363     367     109     476     77	23	109	258     11 GOOD
12	ADGC	Aging Proj  MAP     0	415     0       138     277     415     33	67	138	-139	12 GOOD
13	ADGC	Mayo Clini  MAYO	0	349     0       250     99      349     72	28	99	151     13 GOOD
14	ADGC	Miami Univ  MIA     61	181     19      186     14      200     93	7	14	172     14 MEH
15	ADGC	AD Genetic  MIR     0	284     47      316     15      331     95	5	15	301     15 MEH
16	ADGC	Mayo cl PD  MPD     0	20      0       0       20      20      0	100	0	-20     16 BAD
17	ADGC	NationC AD  NCRD	18	108     52      160     0       160     100	0	0	160     17 BAD
18	ADGC	Wash Unive  RAS     0	0       46      46      0       46      100	0	0	46      18 BAD
19	ADGC	Relig Ordr  ROS     0	351     0       154     197     351     44	56	154	-43     19 GOOD
20	CHARGE	RotterdamS  RS      0	1089	0       276     813     1089	25	75	276	-537	20 GOOD
21	ADGC	Texas AD S  TARC	0	144     0       132     12      144     92	8	12	120     21 MEH
22	ADGC	Un Toronto  TOR     0	0       9       9       0       9       100	0	0	9       22 BAD
23	ADGC	Vanderbilt  VAN     6	235     1       210     26      236     89	11	26	184     23 MEH
24	ADGC	WashNY Age  WCAP	0	147     3       34      116     150     23	77	34	-82     24 OKAY
25	ADGC 	Univ Penns  UPN     40 	0       3       0       0       0       0	0	0	0       25 NONE





GOOD = [1 2 6 9 10 11 12 13 19 20];
OKAY = [5 24];
MEH  = [3 14 15 21 23];
BAD  = [4 7 816 17 18 22];

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
% 1. Identify cohorts in PHEN with a reasonable number of CASE and CTRL
% 2. Take a random selection of people from those cohorts so that
%    there are the same number of CASE and CTRL.
% 3. The rest of the people from the keeper cohorts can be used to
%    increase the size of the NN test group.




% GET SET OF GOOD COHORTS

% COHSET = PHEN(PHEN.GOODCOH==1,:);     % 8824  TOTAL PEOPLE
COHSET = PHEN;                          % 10910 TOTAL PEOPLE


GOOD = [1 2 6 9 10 11 12 13 19 20];
OKAY = [5 24];

i = ismember(COHSET.COHORTNUM, [GOOD OKAY]);
COHSET = COHSET(i,:);






% THE FOLLOWING IS PERFORMED FOR EACH GOOD COHORT SEPARATELY
% ENSURING EQUAL CASE AND CONTROL REPRESENTATION FROM EACH
% COHORT. THE FOLLOWING CODE WILL GENERATE A TRAINING AND
% TESTING GROUP WITH APPROXIMATELY EQUAL NUMBER OF PEOPLE
% FROM EACH GOOD COHORT.



% COUNT NUMBER OF CASE AND CTRL IN EACH COHORT

cohNums = unique(COHSET.COHORTNUM);

nCA=[];nCO=[];
for nn = 1:numel(cohNums)
    nCA(nn) = sum(  (COHSET.COHORTNUM==cohNums(nn)) & (COHSET.AD==1)  );
    nCO(nn) = sum(  (COHSET.COHORTNUM==cohNums(nn)) & (COHSET.AD==0)  );
end

nCACO = [nCA' nCO'];
minCACO = min(nCACO,[],2);





% REMOVE ANY COHORT WITH LESS THAN 10 CASES AND 10 CTRLS

cohNums(minCACO < 10) = [];
nCACO(minCACO < 10,:) = [];
minCACO(minCACO < 10) = [];




% CREATE PHENOTYPE TABLE FOR EACH COHORT

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

rca={};rco={};
for nn = 1:numel(cohNums)
    rca{nn} = randperm( nCACO(nn,1)  , round(minCACO(nn) * .9));
    rco{nn} = randperm( nCACO(nn,2)  , round(minCACO(nn) * .9));
end





% GET ROW INDEX NUMBER FOR PEOPLE NOT CHOSEN ABOVE
%
% Get index locations for people not part of the training subgroup

ica={};ico={};
for nn = 1:numel(cohNums)
    [ia,ib] = ismember(  (1:nCACO(nn,1)) , rca{nn});
    [ic,id] = ismember(  (1:nCACO(nn,2)) , rco{nn});
    ica{nn} = find(~ia);
    ico{nn} = find(~ic);
end




% CREATE PHEN TABLES FOR TRAINING CASE/CTRL & TESTING CASE/CTRL

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


PHETRCASE = vertcat(COHTRANCA{:});
PHETRCTRL = vertcat(COHTRANCO{:});

PHETECASE = vertcat(COHTESTCA{:});
PHETECTRL = vertcat(COHTESTCO{:});



clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
PHETRCASE PHETRCTRL PHETECASE PHETECTRL
disp(PHETRCASE(1:9,:))



%% IDENTIFY KEEPER PHEN AND MAKE CASE:CTRL EQUAL SIZE
%{
% 1. Identify cohorts in PHEN with a reasonable number of CASE and CTRL
% 2. Take a random selection of people from those cohorts so that
%    there are the same number of CASE and CTRL.
% 3. The rest of the people from the keeper cohorts can be used to
%    increase the size of the NN test group.




% GET SET OF GOOD COHORTS

COHSET = PHEN(PHEN.GOODCOH==1,:);     % 8824 TOTAL PEOPLE





% THE FOLLOWING IS PERFORMED FOR EACH GOOD COHORT SEPARATELY
% ENSURING EQUAL CASE AND CONTROL REPRESENTATION FROM EACH
% COHORT. THE FOLLOWING CODE WILL GENERATE A TRAINING AND
% TESTING GROUP WITH APPROXIMATELY EQUAL NUMBER OF PEOPLE
% FROM EACH GOOD COHORT.



% COUNT NUMBER OF CASE AND CTRL IN EACH COHORT

cohNums = unique(COHSET.COHORTNUM);

nCA=[];nCO=[];
for nn = 1:numel(cohNums)
    nCA(nn) = sum(  (COHSET.COHORTNUM==cohNums(nn)) & (COHSET.AD==1)  );
    nCO(nn) = sum(  (COHSET.COHORTNUM==cohNums(nn)) & (COHSET.AD==0)  );
end

nCACO = [nCA' nCO'];
minCACO = min(nCACO,[],2);





% CREATE PHENOTYPE TABLE FOR EACH COHORT

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

rca={};rco={};
for nn = 1:numel(cohNums)
    rca{nn} = randperm( nCACO(nn,1)  , minCACO(nn)-5);
    rco{nn} = randperm( nCACO(nn,2)  , minCACO(nn)-5);
end





% GET ROW INDEX NUMBER FOR PEOPLE NOT CHOSEN ABOVE
%
% Get index locations for people not part of the training subgroup

ica={};ico={};
for nn = 1:numel(cohNums)
    [ia,ib] = ismember(  (1:nCACO(nn,1)) , rca{nn});
    [ic,id] = ismember(  (1:nCACO(nn,2)) , rco{nn});
    ica{nn} = find(~ia);
    ico{nn} = find(~ic);
end




% CREATE PHEN TABLES FOR TRAINING CASE/CTRL & TESTING CASE/CTRL

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


PHETRCASE = vertcat(COHTRANCA{:});
PHETRCTRL = vertcat(COHTRANCO{:});

PHETECASE = vertcat(COHTESTCA{:});
PHETECTRL = vertcat(COHTESTCO{:});



clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
PHETRCASE PHETRCTRL PHETECASE PHETECTRL
disp(PHETRCASE(1:9,:))
%}




%% PUT HALF THE TEST GROUP BACK INTO THE TRAINING GROUP

szCA = round(size(PHETECASE,1)/2);
szCO = round(size(PHETECTRL,1)/2);

PHETRCASE = [PHETRCASE; PHETECASE(1:szCA,:)];
PHETRCTRL = [PHETRCTRL; PHETECTRL(1:szCO,:)];

PHETECASE(1:szCA,:) = [];
PHETECTRL(1:szCO,:) = [];


clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
PHETRCASE PHETRCTRL PHETECASE PHETECTRL




%% IDENTIFY APOE44 INDIVIDUALS

% APOE44_CTRL = [
%      1125223
%      1125232
%      1127783
%      1285294
%      1293115
%      1438861
%      1275619
%      1275795
%      1267683
%      1279016
%      1284834
%      1506358
%      1198253
%      1501559
% ];


APOE44TF = PHEN.APOE == 44;

APOE44TAB = PHEN(APOE44TF,:);

APOE44TAB = sortrows(APOE44TAB,{'AD','GOODCOH'},{'ascend','descend'});

APOE44CTRL = APOE44TAB((APOE44TAB.AD == 0) & (APOE44TAB.GOODCOH == 1),:);



clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
PHETRCASE PHETRCTRL PHETECASE PHETECTRL APOE44CTRL
clc; disp(APOE44CTRL)



%---- REMOVE APOE44 CTRLs FROM TRAIN/TEST TABLES 
%---- THEN ADD ALL APOE44 CTRLs TO TEST TABLE


[Ai,Bi] = ismember(PHETRCTRL.SRR,APOE44CTRL.SRR);
A = PHETRCTRL(Ai,:);
PHETRCTRL(Ai,:) = [];

[Ai,Bi] = ismember(PHETECTRL.SRR,APOE44CTRL.SRR);
B = PHETECTRL(Ai,:);
PHETECTRL(Ai,:) = [];


PHETECTRL = [APOE44CTRL; PHETECTRL];


clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
PHETRCASE PHETRCTRL PHETECASE PHETECTRL APOE44CTRL
clc; disp(PHETECTRL(1:20,1:9))




%% SIMPLY USE 15% OF TRAINING GROUP TO MAKE A TEST GROUP
%{
% szCA = round(size(PHETRCASE,1) * .15);
% szCO = round(size(PHETRCTRL,1) * .15);
% 
% PHETECASE = PHETRCASE(1:szCA,:);
% PHETECTRL = PHETRCTRL(1:szCO,:);
% 
% PHETRCASE(1:szCA,:) = [];
% PHETRCTRL(1:szCO,:) = [];
% 
% 
% clearvars -except ADSP GENB LOCI CASE CTRL PHEN USNP...
% PHETRCASE PHETRCTRL PHETECASE PHETECTRL
%}





%% DICTATE EXACTLY HOW MANY PEOPLE IN TRAINING GROUP
%{
% S = 1500;
% 
% szCA = round(size(PHETRCASE,1));
% szCO = round(size(PHETRCTRL,1));
% 
% dCA = szCA - S;
% dCO = szCO - S;
% 
% 
% PHETECASE = [PHETECASE; PHETRCASE(1:dCA,:)];
% PHETECTRL = [PHETECTRL; PHETRCTRL(1:dCO,:)];
% 
% PHETRCASE(1:dCA,:) = [];
% PHETRCTRL(1:dCO,:) = [];
% 
% 
% clearvars -except ADSP GENB LOCI CASE CTRL PHEN USNP...
% PHETRCASE PHETRCTRL PHETECASE PHETECTRL
%}






%% MAKE SURE EVERYONE CTRL.AGE>75  &  CASE.AGE<90

% THE YOUNGEST APOE44 CTRLS ARE 74, 79, 80, 80, 81...

PHETRCTRL(PHETRCTRL.AGE < 80 , :) = [];
PHETECTRL(PHETECTRL.AGE < 80 , :) = [];

PHETRCASE(PHETRCASE.AGE > 85 , :) = [];
PHETECASE(PHETECASE.AGE > 85 , :) = [];






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

PASS = (LOCI.CASEALT > 10) | (LOCI.CTRLALT > 10);

LOCI  = LOCI(PASS,:);
CASE  = CASE(PASS);
CTRL  = CTRL(PASS);
USNP  = USNP(PASS);


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







%% REMOVE VARIANTS WITH VERY LOW COUNTS ALT ALLELE COUNTS
%{
% PASS = (TRCASEN > 5) & (TRCTRLN > 5);
% TRCASEN = TRCASEN(PASS);
% TRCTRLN = TRCTRLN(PASS);
% LOCI  = LOCI(PASS,:);
% CASE  = CASE(PASS);
% CTRL  = CTRL(PASS);
% USNP  = USNP(PASS);
% LOCI.VID  = (1:size(LOCI,1))';
% 
% clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
% TRCASE TRCTRL TECASE TECTRL
%}



%% REMOVE VARIANTS WHERE ALT > REF
%{
% PASS = (LOCI.TRCASEREF > LOCI.TRCASEALT./1.5) |...
%        (LOCI.TRCTRLREF > LOCI.TRCTRLALT./1.5);
% sum(~PASS)
% 
% LOCI  = LOCI(PASS,:);
% CASE  = CASE(PASS);
% CTRL  = CTRL(PASS);
% USNP  = USNP(PASS);
% LOCI.VID  = (1:size(LOCI,1))';
% 
% clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
% TRCASE TRCTRL TECASE TECTRL
%}








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


% SET MAIN FISHP TO TRAINING GROUP FISHP
VLOCI.FISHP      = VLOCI.TRFISHP;
VLOCI.FISHOR     = VLOCI.TRFISHOR;
VLOCI.CASEREF    = VLOCI.TRCASEREF;
VLOCI.CASEALT    = VLOCI.TRCASEALT;
VLOCI.CTRLREF    = VLOCI.TRCTRLREF;
VLOCI.CTRLALT    = VLOCI.TRCTRLALT;


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



%% ONLY KEEP TWO VARIANTS PER GENE (THE ONES WITH LOWEST P-VAL)


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






%% TAKE THE TOP N GENES FOR NEURAL NET CLASSIFIER TRAINING



% FIRST SORT BY TRAINING GROUP FISHP
[~,i]  = sort(VLOCI.TRFISHP);
VLOCI  = VLOCI(i,:);
VCASE  = VCASE(i);
VCTRL  = VCTRL(i);
VUSNP  = VUSNP(i);





% EXTRACT TOP-N NUMBER OF VARIANTS
N = 120;

VLOCI  = VLOCI(1:N,:);
VCASE  = VCASE(1:N);
VCTRL  = VCTRL(1:N);
VUSNP  = VUSNP(1:N);



clc; disp(VLOCI(1:9,:))
fprintf('\n LOCI COUNT FOR FEATURE MATRIX: %.0f \n\n',size(VLOCI,1))

clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL







%% REMOVE TOP 100 VARIANTS (SPECIFIC TEST)
%{
% VLOCI(1:100,:)   = [];
% VCASE(1:100) = [];
% VCTRL(1:100) = [];
% VUSNP(1:100) = [];
% VLOCI.VID        = (1:size(VLOCI,1))';
%}








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




%##########################################################################
%
%%                  IMPORT 23 AND ME DATA
%
%##########################################################################
%{
clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE...
FULLTRX FULLTEX TRX TRL TEX TEL


ttme_path = '/Users/bradleymonk/Documents/MATLAB/GIT/genomics/genos/genosdata/genome23andme_noxy.txt';

opts = detectImportOptions(ttme_path); % opts = setvaropts;
TTME = readtable(ttme_path,opts);

TTME.POS = uint64(TTME.POS);
TTME.CHR = uint64(TTME.CHR);

POS = string(TTME.POS);
POS = convertStringsToChars(POS);
N = cellfun(@(x) numel(x),POS);
M = 10-N;

CHR = double(TTME.CHR);
CHR = CHR.*(10.^M);
CHR = uint64(CHR);
CHR = string(CHR);
CHR = convertStringsToChars(CHR);

% POS = convertCharsToStrings(POS);
% CHR = convertCharsToStrings(CHR);

POS = char(POS);
CHR = char(CHR);

CHRPOS = horzcat(CHR,POS);
CHRPOS = string(CHRPOS);
CHRPOS = strrep(CHRPOS, ' ', '');

TTME.CHRPOS = str2num(char(CHRPOS));
TTME.CHRPOS = uint64(TTME.CHRPOS);


clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE...
FULLTRX FULLTEX TRX TRL TEX TEL TTME



%% CHECK IF ANY CHRPOS FROM 23ANDME MATCHES VLOCI CHRPOS


TCHRPOS = TTME.CHRPOS;
VCHRPOS = VLOCI.CHRPOS;


% [Ai,Bi] = ismember(TCHRPOS,VCHRPOS);
[Ai,Bi] = ismember(VCHRPOS,TCHRPOS);

[i,j] = find(Ai);

sum(Ai)
sum(Bi)

TTSEQ = zeros(size(VCHRPOS));
TTSEQ = Ai;

BP = TTME.BP;



A = VCHRPOS;
B = TCHRPOS;
C = cellstr(num2str(uint64(zeros(size(A)))));


[Ai,Bi] = ismember(A,B);
Bj = Bi(Bi>0);
% B(Bj)
C(Bi>0)=BP(Bj);

C = strrep(C, '0', 'NN');
C = strrep(C, '-', 'N');
C = strrep(C, 'I', 'I');
C = strrep(C, 'D', 'D');

VLOCI.TTME = C;


% VN = VLOCI.Properties.VariableNames;

clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE...
FULLTRX FULLTEX TRX TRL TEX TEL TTME



%% DETERMINE THE REF / ALT STATUS OF 23ANDME PATIENT


REF = VLOCI.REF;
ALT = VLOCI.ALT;
TME = VLOCI.TTME;

TM = split(TME,"");
TMa = TM(:,2);
TMb = TM(:,3);


REFa = strcmp(TMa,REF).*-1;
REFb = strcmp(TMb,REF).*-1;
ALTa = strcmp(TMa,ALT).*3;
ALTb = strcmp(TMb,ALT).*3;

TTMESNP = REFa + REFb + ALTa + ALTb;

% SNP(SNP==0) = -2;     % HOMREF: -2
% SNP(SNP==1) =  0;     % UNCALL:  0
% SNP(SNP==2) =  2;     % HETALT:  2
% SNP(SNP==3) =  6;     % HOMALT:  6



clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE...
TTME TTMESNP
%}






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

[FULLTRX, TRX, TRL] = makeapoenet(VLOCI,VCASE,VCTRL,VUSNP,TRPHE,1);

[FULLTEX, TEX, TEL] = makeapoenet(VLOCI,VCASE,VCTRL,VUSNP,TEPHE,1);



% TRX & TEX are  [Variant (rows) x Person (cols)]
% 
% FULLTRX & FULLTEX are  [Person (rows) x Variant (cols)]
%     Top row of FULL(1,:) is CHRPOS
%     First 4 cols of FULL(:,1:4) are SRR|AD|AGE|APOE





clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE...
FULLTRX FULLTEX TRX TRL TEX TEL TTME TTMESNP









%###################################################################
%%           BUILD NONLINEAR INTERACTIONS MATRIX
%###################################################################
% TRX(TRX==0) = -1;
% TEX(TEX==0) = -1;
% 
% TRX(1,:) = (TRPHE.AD==1)' .* 4 - 2;
% TEX(1,:) = (TEPHE.AD==1)' .* 4 - 2;
% 
% disp((FULLTRX(2:11,2))')
% disp(TRX(1:10,1:10))



% TRX    &  TEX    are  [Variant (rows) x Person  (cols)]
% TRPV   &  TEPV   are  [Person  (rows) x Variant (cols)]
% TRPVNL &  TEPVNL are  [Person  (rows) x Variant (cols)]

TRPV = TRX';
TEPV = TEX';

TRPVNL = TRPV;
TEPVNL = TEPV;


% TRAINING MATRIX WILL BE EXPANDED TO...
cs = cumsum(1:size(TRPV,2));
disp(cs(end))

% For each variant...
for n = 1:(size(TRPV,2)-1)

    TRPVNL = [TRPVNL ( TRPV(:,n) .* TRPV(:,(n+1):end) ) ];
    TEPVNL = [TEPVNL ( TEPV(:,n) .* TEPV(:,(n+1):end) ) ];

end

% Include X^2 Effect
TRPVNL = [TRPVNL TRPVNL.^2];
TEPVNL = [TEPVNL TEPVNL.^2];



% Normalize matrix so min value is 1
% TRPVNL = TRPVNL - min(TRPVNL(:)) + 1;
% TEPVNL = TEPVNL - min(TEPVNL(:)) + 1;
% TRPVNL = log(TRPVNL);
% TEPVNL = log(TEPVNL);
% imagesc(log(TRPVNL)')

TRX = TRPVNL';
TEX = TEPVNL';


clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE...
FULLTRX FULLTEX TRX TRL TEX TEL







%%
return
%##########################################################################
%
%%                  PERFORM LINEAR REGRESSION
%
%##########################################################################

% TRX  &  TEX  are  [Variant (rows) x Person  (cols)]
% RX   &  RL   are  [Person  (rows) x Variant (cols)]

RX = TRX';
EX = TEX';

RL = TRL';
EL = TEL';
% RL = TRL' - .5;
% EL = TEL' - .5;


% Create labels CASE=1 CTRL=0
TR_y = RL(:,1);


% Add intercept column
TR_x = [ones(size(RX,1),1) RX];


% Perform the normal function calculation
t = pinv(TR_x' * TR_x) * (TR_x' * TR_y);


% Get linear model predictions for each y:
TR_fx = sum( TR_x .* t' ,2);


% Descretize predictions
TR_yp = round(TR_fx);
% TR_yp = TR_fx > 0;


TRmu = mean(TR_yp == TR_y);
% TRmu = mean(TR_yp == (TR_y>0));


% Eval performance on high confidence outputs
TRHi = (TR_fx>.8) | (TR_fx<.2);
% TRHi = abs(TR_fx)>.4;

TRmuHi = mean(TR_yp(TRHi) == TR_y(TRHi) );
% TRmuHi = mean(TR_yp(TRHi) == (TR_y(TRHi)>0) );






%% Evaluate performance of hold-out test set
EX = TEX';
EL = TEL';


TE_x = [ones(size(EX,1),1) EX];
TE_y = EL(:,1);

TE_fx = sum( TE_x .* t' ,2);
TE_yp = TE_fx>.5;

TEmu = mean(TE_yp == TE_y)


% Eval performance on high confidence outputs
TEHi = (TE_fx>.8) | (TE_fx<.2);

TEmuHi = mean(TE_yp(TEHi) == TE_y(TEHi))



%% PRINT OUTPUTS TO CON
clc;
fprintf('\n TRAIN Pct correct all: %.0f ',ceil(TRmu*100))
fprintf('\n TRAIN Pct pop all: %.0f \n',100)

fprintf('\n TRAIN Pct correct hicon: %.0f ',ceil(TRmuHi*100))
fprintf('\n TRAIN Pct pop hicon: %.0f \n\n',ceil(mean(TRHi)*100))

fprintf('\n TEST Pct correct all: %.0f ',ceil(TEmu*100))
fprintf('\n TEST Pct pop all: %.0f \n',100)

fprintf('\n TEST Pct correct hicon: %.0f ',ceil(TEmuHi*100))
fprintf('\n TEST Pct pop hicon: %.0f \n\n',ceil(mean(TEHi)*100))


clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE...
FULLTRX FULLTEX TRX TRL TEX TEL TTME TTMESNP TR_fx TE_fx t ...
TRmu TRmuHi TRHi TEmu TEmuHi TEHi ME_x ME_Act ME_fx ME_yp




















%------------------------------------------------
%%     EVALUATE THE APOE44 CTRLS
%------------------------------------------------

% Check if people in PHEN file reported as e44 are actually e44
FRX = FULLTRX(2:end,:)';
A = FRX(4,:)==44;
FRX = FRX(:,A);
FRX_NOT44 = (FRX(4,:) ~= 44) | (FRX(5,:) ~= 6);
disp(sum(  FRX_NOT44  ))

FRX(:,FRX_NOT44) = [];

TRPHE44 = TRPHE(A,:);
TRPHE44(FRX_NOT44,:) = [];




FEX = FULLTEX(2:end,:)';
A = FEX(4,:)==44;
FEX = FEX(:,A);
FEX_NOT44 = (FEX(4,:) ~= 44) | (FEX(5,:) ~= 6);
disp(sum(  FEX_NOT44  ))

FEX(:,FEX_NOT44) = [];

TEPHE44 = TEPHE(A,:);
TEPHE44(FEX_NOT44,:) = [];




% SEPARATE OUT APOE44 CASES AND CTRLS.
FRX44CA = FRX(:,FRX(2,:)==1);
FRX44CO = FRX(:,FRX(2,:)==0);

TRPHE44CA = TRPHE44(TRPHE44.AD==1,:);
TRPHE44CO = TRPHE44(TRPHE44.AD==0,:);



FEX44CA = FEX(:,FEX(2,:)==1);
FEX44CO = FEX(:,FEX(2,:)==0);

TEPHE44CA = TEPHE44(TEPHE44.AD==1,:);
TEPHE44CO = TEPHE44(TEPHE44.AD==0,:);


disp(' ')
disp(FEX44CA(1:5,:))
disp(' ');disp(' ');
disp(FEX44CO(1:5,:))


TE44CA = FEX44CA(5:end,:)';
TE44CO = FEX44CO(5:end,:)';


% Evaluate regression classifier on TE44CO

TE44CA_x = [ones(size(TE44CA,1),1) TE44CA];
TE44CO_x = [ones(size(TE44CO,1),1) TE44CO];


TE44CA_y =  ones(size(TE44CA,1),1);
TE44CO_y = zeros(size(TE44CO,1),1);


TE44CA_Act = TE44CA_x .* t';
TE44CO_Act = TE44CO_x .* t';


TE44CA_fx = sum( TE44CA_Act ,2);
TE44CO_fx = sum( TE44CO_Act ,2);


TE44CA_yp = round(TE44CA_fx);
TE44CO_yp = round(TE44CO_fx);


TE44CA_mu = mean(TE44CA_yp == TE44CA_y);
TE44CO_mu = mean(TE44CO_yp == TE44CO_y);



cai = find(~TE44CA_yp)
sum(TE44CA_Act(cai,:))
bar(TE44CA_Act(cai,:))



% VLOCI.Beta  = t(2:end);


disp(' ')
disp([FEX44CO(1:5,:); TE44CO_fx' ])




%%
clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE...
FULLTRX FULLTEX TRX TRL TEX TEL TTME TTMESNP TR_fx TE_fx ...
TRmu TRmuHi TRHi TEmu TEmuHi TEHi ME_x ME_Act ME_fx ME_yp t...
TE44CO TE44CA TE44CA_Act TE44CO_Act TE44CA_fx TE44CO_fx...
TE44CA_yp TE44CO_yp TE44CA_mu TE44CO_mu





%% DISPLAY AND PLOT PERFORMANCE OF E44 CONTROLS
close all;
fh1 = figure('Units','normalized','OuterPosition',[.01 .06 .9 .88],'Color','w','MenuBar','none');
ax1 = axes('Position',[.05 .53 .9 .42],'Color','none'); box off; hold on;
ax2 = axes('Position',[.05 .05 .9 .42],'Color','none'); box off; hold on;

axes(ax1)
bar(t(3:120))
xlim([-2 122])

axes(ax2)
bar((TE44CO_Act(:,3:120))')
xlim([-2 122])


clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE...
FULLTRX FULLTEX TRX TRL TEX TEL TTME TTMESNP TR_fx TE_fx ...
TRmu TRmuHi TRHi TEmu TEmuHi TEHi ME_x ME_Act ME_fx ME_yp t...
TE44CO TE44CA TE44CA_Act TE44CO_Act TE44CA_fx TE44CO_fx...
TE44CA_yp TE44CO_yp TE44CA_mu TE44CO_mu



%% DISPLAY AND PLOT PERFORMANCE OF E44 CASES VS E44 CONTROLS
close all;
fh1 = figure('Units','normalized','OuterPosition',[.01 .06 .9 .88],'Color','w','MenuBar','none');
ax1 = axes('Position',[.05 .53 .9 .42],'Color','none'); box off; hold on;
ax2 = axes('Position',[.05 .05 .9 .42],'Color','none'); box off; hold on;


ns = 3;
ne = 150;

axes(ax1)
bar((TE44CO_Act(:,ns:ne))')
xlim([-2 ne+2])

axes(ax2)
bar((TE44CA_Act(:,ns:ne))')
xlim([-2 ne+2])


clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE...
FULLTRX FULLTEX TRX TRL TEX TEL TTME TTMESNP TR_fx TE_fx ...
TRmu TRmuHi TRHi TEmu TEmuHi TEHi ME_x ME_Act ME_fx ME_yp t...
TE44CO TE44CA TE44CA_Act TE44CO_Act TE44CA_fx TE44CO_fx...
TE44CA_yp TE44CO_yp TE44CA_mu TE44CO_mu



%% DETERMINE WHAT IS UNIQUE ABOUT LOW-CONF APOE44 INDIVIDUALS


APOEstd = std(TE44CO)';

[x,i] = sort(APOEstd,'descend');

LC  = TE44CO_fx < .7
LCi = find(LC)


LOCO = TE44CO(find(LC),:);
HICO = TE44CO(find(~LC),:);

LOMU = mean(LOCO)
HIMU = mean(HICO)

MUDIFF = round(HIMU - LOMU);


% 
% 
% imagesc(TE44CO)
% t
% ME_x
% ME_Act
% ME_fx




% DISPLAY AND PLOT PERFORMANCE
clc;close all;
fh1 = figure('Units','normalized','OuterPosition',[.01 .06 .9 .6],'Color','w','MenuBar','none');
ax1 = axes('Position',[.05 .53 .9 .42],'Color','none');
ax2 = axes('Position',[.05 .05 .9 .42],'Color','none');

TECO = TE44CO(:,1:100);

axes(ax1)
imagesc(TECO)
ylim([.5 11.5])

ax1.YTick = .5:size(TECO,1)+1;
YL1 = cellstr(['\fontsize{1}  0']);
YL2 = cellstr([repmat('\fontsize{14}    ',size(TECO,1),1), num2str((1:11)')])
YLabs = [YL1; YL2]
ax1.YTickLabel = cellstr(YLabs);
% ax1.FontName = 'Courier';
% ax1.TickDir = 'out'
ax1.YAxis.TickDirection = 'out';
ax1.YTickLabelRotation = 90

ax1.YGrid = 'on';
ax1.GridColor = [.99 .33 .33]
ax1.GridAlpha = .8
ax1.GridLineStyle = ':';
ax1.LineWidth = 1



TECO = TE44CO(:,101:200);

axes(ax2)
imagesc(TECO)
ylim([.5 11.5])

ax2.YTick = .5:size(TECO,1)+1;
YL1 = cellstr(['\fontsize{1}  0']);
YL2 = cellstr([repmat('\fontsize{14}    ',size(TECO,1),1), num2str((1:11)')])
YLabs = [YL1; YL2]
ax2.YTickLabel = cellstr(YLabs);
% ax1.FontName = 'Courier';
% ax2.TickDir = 'out'
ax2.YAxis.TickDirection = 'out';
ax2.YTickLabelRotation = 90

ax2.YGrid = 'on';
ax2.GridColor = [.99 .33 .33]
ax2.GridAlpha = .8
ax2.GridLineStyle = ':';
ax2.LineWidth = 1

% ax2.XTick = 1:100;
% XL1 = cellstr([repmat('\fontsize{9}    ',100,1), num2str((1:100)')])
% ax2.XTickLabel = cellstr(XL1);



clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE...
FULLTRX FULLTEX TRX TRL TEX TEL TTME TTMESNP TR_fx TE_fx ...
TRmu TRmuHi TRHi TEmu TEmuHi TEHi ME_x ME_Act ME_fx ME_yp t...
TE44CO TE44CA TE44CA_Act TE44CO_Act TE44CA_fx TE44CO_fx...
TE44CA_yp TE44CO_yp TE44CA_mu TE44CO_mu



%% FIT GLM

RX = TRX';
RL = TRL';
EX = TEX';
EL = TEL';

y = nnlab(RL);

[b,dev,stats] = glmfit(RX,y,'binomial');







% Get linear model predictions for each y:
TR_x = [ones(size(RX,1),1) RX];

TR_fx = sum( TR_x .* b' ,2);


% Descretize predictions
TR_yp = round(TR_fx);

TRmu = mean((TR_fx>.12) == y)


% Eval performance on high confidence outputs
TRHi = (TR_fx>.8) | (TR_fx<.2);

TRmuHi = mean(TR_yp(TRHi) == TR_y(TRHi));








%################################################################
%%  ASSIGN NEURAL NET ACTIVATION VALUE TO PARTICIPANT PHE
%################################################################

% numnets = size(netz,2);
% ACTR = zeros(size(TRX,2),numnets);
% ACTE = zeros(size(TEX,2),numnets);
% 
% for nn = 1:numnets
%     net = netz{nn};
%     ACTI = net(TRX);
%     ACTR(:,nn) = (ACTI(1,:))';
%     ACTI = net(TEX);
%     ACTE(:,nn) = (ACTI(1,:))';
% end



BRAAK = TRPHE.BRAAK;

BRAAK(TRPHE.AUTOPSY==0) = NaN;


clc;

nanmean(BRAAK(TRPHE.MUACTR > .90))

nanmean(BRAAK( ((TRPHE.MUACTR<.90)&(TRPHE.MUACTR>.80))  ))

nanmean(BRAAK( ((TRPHE.MUACTR<.80)&(TRPHE.MUACTR>.70))  ))

nanmean(BRAAK( ((TRPHE.MUACTR<.70)&(TRPHE.MUACTR>.60))  ))

nanmean(BRAAK( ((TRPHE.MUACTR<.60)&(TRPHE.MUACTR>.50))  ))

nanmean(BRAAK(TRPHE.MUACTR < .5))









%##########################################################################
%
%%          MACHINE LEARNING USING ARTIFICIAL NEURAL NETWORKS
%
%##########################################################################



NN = patternnet([200 100 100]);

net = train(NN,TRX,TRL);

[NNPerf] = netstats(net,TRX,TRL,TEX,TEL,.05,.85,.15,1);




% clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
% VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE...
% FULLTRX FULLTEX TRX TRL TEX TEL net netz


clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE...
FULLTRX FULLTEX TRX TRL TEX TEL net netz NETi ...
TTME TTMESNP TR_fx TE_fx t ...
TRmu TRmuHi TRHi TEmu TEmuHi TEHi





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
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE...
FULLTRX FULLTEX TRX TRL TEX TEL net netz NETi ...
TTME TTMESNP TR_fx TE_fx t ...
TRmu TRmuHi TRHi TEmu TEmuHi TEHi











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
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE...
FULLTRX FULLTEX TRX TRL TEX TEL net netz NETi ...
TTME TTMESNP TR_fx TE_fx t ...
TRmu TRmuHi TRHi TEmu TEmuHi TEHi








%% PLOT ROC CURVE


simTRAIN = sim(net,TRX);
simTEST  = sim(net,TEX);

[TPR,FPR,THR] = roc(TRL,simTRAIN);
% [TPR,FPR,THR] = roc(TEL,simTEST);



TPR = cell2mat(TPR');
FPR = cell2mat(FPR');
THR = cell2mat(THR');



%--- Plot ROC
clc; close all;
fh1 = figure('Units','normalized','Position',[.05 .04 .9 .9],'Color','w','MenuBar','none');
ax1 = axes('Position',[.06 .06 .9 .9],'Color','none');

ph1 = plot(FPR',TPR','LineWidth',3);

%--- Add Lines and Grid to plot
ph2 = line([0 1],[0 1],'Color','k','LineStyle','--','LineWidth',1);
box on; grid on


%--- Add Legend
legend(ax1,{'Case','Control'})




clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE...
FULLTRX FULLTEX TRX TRL TEX TEL net netz NETi TRC TEC TRL_COH TEL_COH...
TTME TTMESNP TR_fx TE_fx t ...
TRmu TRmuHi TRHi TEmu TEmuHi TEHi






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



%% MAKE TEST SET X% CASES Y% CTRL
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

