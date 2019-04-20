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

%ID CONSORT	STUDY_NAME  STUDY   FAM PEOPLE  ENR     CASE	CTRL	TOT    %CA %CO  EQL	CULL	NOTES
01  DGC     Adult Chng  ACT     0	1277	0       323     945     1268	25	75	323	-622	KEEP
02  ADGC	AD Centers  ADC     0	3263	0       2438	817     3255	75	25	817	1621	KEEP
03  CHARGE	Athrosclro  ARIC	0	57      0       39      18      57      68	32	18	21      REMOVE
04  CHARGE	Aus Stroke	ASPS	0	126     0       121     5       126     96	4	5	116     REMOVE
05  ADGC	Chic Aging  CHAP	0	231     0       27      204     231     12	88	27	-177	REMOVE *5*
06  CHARGE	CardioHlth  CHS     0	840     0       250     583     833     30	70	250	-333	KEEP
07  ADGC	Hispanic S  CUHS	361	0       331     160     171     331     48	52	160	-11     LATINO *7*
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
23	ADGC	Vanderbilt  VAN     6	235     1       210     26      236     89	11	26	184     REMOVE *23*
24	ADGC	WashNY Age  WCAP	0	147     3       34      116     150     23	77	34	-82     REMOVE *24*
25	ADGC 	Univ Penns  UPN     40 	0       3       0       0       0       0	0	0	0       REMOVE



Along the way, the code will take steps to mitigate stratification
issues due to COHORT-specific effects. 

...I'm trying to think if there is anything else you should know...
I can't, sooo, I guess we should just go-ahead and start with step #1.
(aka the first chance for something to go wrong)...

%}

%% STEP-1: LOAD THE DATASET

close all; clear; clc; rng('shuffle');
genosdir = fileparts(which('GENOS_NEURALNET_POSTHOC.m'));
cd(genosdir);


subfuncpath = [genosdir '/genosfunctions'];
datasetpath = [genosdir '/genosdata'];
gpath = [genosdir ':' subfuncpath ':' datasetpath];
addpath(gpath)


% MATDATA = 'ADSP_ALT_UNC.mat';
MATDATA = 'ADSP.mat';
which(MATDATA)
load(MATDATA)


clc; rng('shuffle');
clearvars -except ADSP



%% CARBON COPY MAIN VARIABLES FROM ADSP.STRUCT


LOCI = ADSP.LOCI;
CASE = ADSP.CASE;
CTRL = ADSP.CTRL;
USNP = ADSP.UNSNP;
PHEN = ADSP.PHEN;

clearvars -except ADSP LOCI CASE CTRL USNP PHEN



%% CLEAN-UP VARIANT TABLE


clc; clearvars -except ADSP LOCI CASE CTRL USNP PHEN
disp(LOCI(1:9,:))




%% IDENTIFY KEEPER PHEN AND MAKE CASE:CTRL EQUAL SIZE
% 1. Identify cohorts in PHEN with a reasonable number of CASE and CTRL
% 2. Take a random selection of people from those cohorts so that
%    there are the same number of CASE and CTRL.
% 3. The rest of the people from the keeper cohorts can be used to
%    increase the size of the NN test group.



% TAG GOOD PHEN
GOODCOH = [1 2 6 9 10 11 12 13 19 20]';
MEHHCOH = [5 7 23 24]'; % 7 IS HISPANIC STUDY


[ai,bi] = ismember(PHEN.COHORTNUM,GOODCOH);
% [ai,bi] = ismember(PHEN.COHORTNUM,[GOODCOH; MEHHCOH]);



PHEN.GOODCOH = ai;


% PUT GOOD PHEN INTO SEPARATE SET
COHSET = PHEN(PHEN.GOODCOH==1,:);     % 8824 TOTAL PEOPLE




% COUNT NUMBER OF CASE AND CTRL IN GOOD PHEN
cohNums = unique(COHSET.COHORTNUM);

nCA=[];nCO=[];
for nn = 1:numel(cohNums)
    nCA(nn) = sum(  (COHSET.COHORTNUM==cohNums(nn)) & (COHSET.AD==1)  );
    nCO(nn) = sum(  (COHSET.COHORTNUM==cohNums(nn)) & (COHSET.AD==0)  );
end

nCACO = [nCA' nCO'];
minCACO = min(nCACO,[],2);



% CREATE COHORT TABLE FOR EACH COHORT
COHCASE={};COHCTRL={};
for nn = 1:numel(cohNums)
    COHCASE{nn} = COHSET( (COHSET.COHORTNUM==cohNums(nn)) & (COHSET.AD==1) ,:);
    COHCTRL{nn} = COHSET( (COHSET.COHORTNUM==cohNums(nn)) & (COHSET.AD==0) ,:);
end



% GET RANDOM ROW NUMBER FOR EACH COHORT, EQUAL CASE AND CTRL COUNT PER COHORT
% Get random permutation of M values between 1:N-5
%    where  M = min CASE\CTRL group size per cohort
%           N = total cohort size
rca={};rco={};
for nn = 1:numel(cohNums)
rca{nn} = randperm( nCACO(nn,1)  , minCACO(nn)-3);
rco{nn} = randperm( nCACO(nn,2)  , minCACO(nn)-3);
end



% GET ROW INDEX NUMBER FOR PEOPLE NOT CHOSEN ABOVE
% Get index locations for people not part of the training subgroup
ica={};ico={};
for nn = 1:numel(cohNums)
[ia,ib] = ismember(  (1:nCACO(nn,1)) , rca{nn});
[ic,id] = ismember(  (1:nCACO(nn,2)) , rco{nn});
ica{nn} = find(~ia);
ico{nn} = find(~ic);
end



% Make the final set of people in the training and testing groups
COHTRANCA={};COHTRANCO={};
COHTESTCA={};COHTESTCO={};
for nn = 1:numel(cohNums)
    COHTRANCA{nn} = COHCASE{nn}(rca{nn},:);
    COHTRANCO{nn} = COHCTRL{nn}(rco{nn},:);

    COHTESTCA{nn} = COHCASE{nn}(ica{nn},:);
    COHTESTCO{nn} = COHCTRL{nn}(ico{nn},:);
end



% Turn the cell array back into a single table
COHTRANCASE = vertcat(COHTRANCA{:});
COHTRANCTRL = vertcat(COHTRANCO{:});

COHTESTCASE = vertcat(COHTESTCA{:});
COHTESTCTRL = vertcat(COHTESTCO{:});



clearvars -except ADSP LOCI CASE CTRL USNP PHEN...
COHTRANCASE COHTRANCTRL COHTESTCASE COHTESTCTRL
disp(COHTRANCASE(1:9,:))






%% PUT HALF THE TEST GROUP BACK INTO THE TRAINING GROUP

szCA = round(size(COHTESTCASE,1)/2);
szCO = round(size(COHTESTCTRL,1)/2);

COHTRANCASE = [COHTRANCASE; COHTESTCASE(1:szCA,:)];
COHTRANCTRL = [COHTRANCTRL; COHTESTCTRL(1:szCO,:)];

COHTESTCASE(1:szCA,:) = [];
COHTESTCTRL(1:szCO,:) = [];


clearvars -except ADSP LOCI CASE CTRL USNP PHEN...
COHTRANCASE COHTRANCTRL COHTESTCASE COHTESTCTRL




%% SIMPLY USE 15% OF TRAINING GROUP TO MAKE A TEST GROUP

% szCA = round(size(COHTRANCASE,1) * .15);
% szCO = round(size(COHTRANCTRL,1) * .15);
% 
% COHTESTCASE = COHTRANCASE(1:szCA,:);
% COHTESTCTRL = COHTRANCTRL(1:szCO,:);
% 
% COHTRANCASE(1:szCA,:) = [];
% COHTRANCTRL(1:szCO,:) = [];
% 
% 
% clearvars -except ADSP GENB LOCI CASE CTRL USNP PHEN...
% COHTRANCASE COHTRANCTRL COHTESTCASE COHTESTCTRL



%% DICTATE EXACTLY HOW MANY PEOPLE IN TRAINING GROUP

% S = 1500;
% 
% szCA = round(size(COHTRANCASE,1));
% szCO = round(size(COHTRANCTRL,1));
% 
% dCA = szCA - S;
% dCO = szCO - S;
% 
% 
% COHTESTCASE = [COHTESTCASE; COHTRANCASE(1:dCA,:)];
% COHTESTCTRL = [COHTESTCTRL; COHTRANCTRL(1:dCO,:)];
% 
% COHTRANCASE(1:dCA,:) = [];
% COHTRANCTRL(1:dCO,:) = [];
% 
% 
% clearvars -except ADSP GENB LOCI CASE CTRL USNP PHEN...
% COHTRANCASE COHTRANCTRL COHTESTCASE COHTESTCTRL





%% MAKE SURE EVERYONE IN CTRL GROUP IS OVER 75 YEARS OLD

COHTRANCTRL(COHTRANCTRL.AGE < 75 , :) = [];
COHTESTCTRL(COHTESTCTRL.AGE < 75 , :) = [];







%###############################################################
%%  USE PER-PERSON VARIANT COUNTS TO IDENTIFY OUTLIERS
%###############################################################

% The total number of variants has already been counted for each
% person in the dataset and listed in the 'PHEN' table under
% 'TOTvars'. If you ever need to recount, a custom function has
% been coded to do this for you...
%
%   [VNUM] = countvperper(SRR, AD, COHORTNUM, CASE, CTRL);
%
% Here we will use the per-person counts already performed and
% dismiss the individuals on either extreme.


Qlo = quantile(COHTRANCASE.TOTvars,.01);
Qhi = quantile(COHTRANCASE.TOTvars,.99);
TRCASE  = COHTRANCASE(((COHTRANCASE.TOTvars > Qlo) & (COHTRANCASE.TOTvars < Qhi)),:);


Qlo = quantile(COHTRANCTRL.TOTvars,.01);
Qhi = quantile(COHTRANCTRL.TOTvars,.99);
TRCTRL  = COHTRANCTRL(((COHTRANCTRL.TOTvars > Qlo) & (COHTRANCTRL.TOTvars < Qhi)),:);


Qlo = quantile(COHTESTCASE.TOTvars,.01);
Qhi = quantile(COHTESTCASE.TOTvars,.99);
TECASE  = COHTESTCASE(((COHTESTCASE.TOTvars > Qlo) & (COHTESTCASE.TOTvars < Qhi)),:);


Qlo = quantile(COHTESTCTRL.TOTvars,.01);
Qhi = quantile(COHTESTCTRL.TOTvars,.99);
TECTRL  = COHTESTCTRL(((COHTESTCTRL.TOTvars > Qlo) & (COHTESTCTRL.TOTvars < Qhi)),:);



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

clearvars -except ADSP  LOCI CASE CTRL PHEN USNP...
COHTRANCASE COHTRANCTRL COHTESTCASE COHTESTCTRL...
TRCASE TRCTRL TECASE TECTRL








%% FILTER VARIANTS WHERE CASEALT<5 & CTRLALT<5

PASS = (LOCI.CASEALT > 5) &...
       (LOCI.CTRLALT > 5);
sum(~PASS)

LOCI  = LOCI(PASS,:);
CASE  = CASE(PASS);
CTRL  = CTRL(PASS);
USNP  = USNP(PASS);
LOCI.VID  = (1:size(LOCI,1))';


clearvars -except ADSP GENB LOCI CASE CTRL PHEN USNP...
COHTRANCASE COHTRANCTRL COHTESTCASE COHTESTCTRL...
TRCASE TRCTRL TECASE TECTRL TRCASEN TRCTRLN TECASEN TECTRLN





%###############################################################
%%          COUNT NUMBER OF VARIANTS PER LOCI
%###############################################################
% The varsum() function will go through each known variant loci
% and check whether anyone's SRR ID from your subset of IDs match
% all known SRR IDs for that loci. It will then sum the total
% number of alleles (+1 for hetzy-alt, +2 for homzy-alt) for each
% loci and return the totals.


[TRCASEN, TRCTRLN] = varsum(CASE, TRCASE.SRR, CTRL, TRCTRL.SRR);

[TECASEN, TECTRLN] = varsum(CASE, TECASE.SRR, CTRL, TECTRL.SRR);


[TRCASEUN, TRCTRLUN, TECASEUN, TECTRLUN] = uncsum(USNP,...
    TRCASE.SRR, TRCTRL.SRR,TECASE.SRR, TECTRL.SRR);



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
COHTRANCASE COHTRANCTRL COHTESTCASE COHTESTCTRL...
TRCASE TRCTRL TECASE TECTRL TRCASEN TRCTRLN TECASEN TECTRLN










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



[~,i] = sort(LOCI.TRFISHP);
LOCI  = LOCI(i,:);
CASE  = CASE(i);
CTRL  = CTRL(i);
USNP  = USNP(i);
LOCI.VID = (1:size(LOCI,1))';


clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
TRCASE TRCTRL TECASE TECTRL




%% COMPUTE CHI SQUARE VALUE



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




%###############################################################
%###############################################################
%###############################################################
%###############################################################
%% PREPARE DATA FOR NEURAL NET CLASSIFIER SUPERVISED LEARNING
%###############################################################
%###############################################################
%###############################################################
%###############################################################
clc;




%% STORE VARIABLES FOR NEURAL NETWORK TRAINING AS 'AMX'

AMX         = LOCI;
AMXCASE     = CASE;
AMXCTRL     = CTRL;
AMXUSNP     = USNP;
AMXTRCASE   = TRCASE;
AMXTRCTRL   = TRCTRL;
AMXTECASE   = TECASE;
AMXTECTRL   = TECTRL;



clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
AMX AMXCASE AMXCTRL AMXUSNP AMXTRCASE AMXTRCTRL AMXTECASE AMXTECTRL




%--------------------------------------------------------------------------
% REMOVE KNOWN AD ALZ GENES BY NAME


% ALZGENES = string({"APOE";"BIN1";"CLU";"ABCA7";"CR1";...
%                    "PICALM";"MS4A6A";"CD33";"MS4A4E";"CD2AP"});
% 
% for nn = 1:numel(ALZGENES)
% 
%     x = strcmp(AMX.GENE,ALZGENES{nn});
%     sum(x)
% 
%     AMX(x,:) = [];
%     AMXCASE(x) = [];
%     AMXCTRL(x) = [];
%     AMXUSNP(x) = [];
% 
% end
% AMX.VID  = (1:size(AMX,1))';



% ALZGENES = string({"APOE"});
% 
% for nn = 1:numel(ALZGENES)
% 
%     x = strcmp(AMX.GENE,ALZGENES{nn});
%     sum(x)
% 
%     AMX(x,:) = [];
%     AMXCASE(x) = [];
%     AMXCTRL(x) = [];
%     AMXUSNP(x) = [];
% 
% end
% AMX.VID  = (1:size(AMX,1))';






%--------------------------------------------------------------------------
% REMOVE SELECT GENES OR JUNK GENES BY NAME

REMGENES = string({"TYRO3";"MUC21";"MUC8";...
                   "MUC6";"MUC16";"MUC40";"MUC4";"MUC20"});

% REMGENES = string({"TOMM40";"TYRO3";"MUC21";"MUC8";...
%                    "MUC6";"MUC16";"MUC40";"MUC4";"MUC20"});

for nn = 1:numel(REMGENES)

    x = strcmp(AMX.GENE,REMGENES{nn});
    sum(x)

    AMX(x,:) = [];
    AMXCASE(x) = [];
    AMXCTRL(x) = [];
    AMXUSNP(x) = [];

end
AMX.VID  = (1:size(AMX,1))';


clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
AMX AMXCASE AMXCTRL AMXUSNP AMXTRCASE AMXTRCTRL AMXTECASE AMXTECTRL









%--------------------------------------------------------------------------
% TAKE THE TOP N GENES FOR NEURAL NET CLASSIFIER TRAINING


% FIRST SORT BY TRAINING GROUP FISHP
[~,i] = sort(AMX.TRFISHP);
AMX      = AMX(i,:);
AMXCASE  = AMXCASE(i);
AMXCTRL  = AMXCTRL(i);
AMXUSNP  = AMXUSNP(i);
AMX.VID  = (1:size(AMX,1))';




N = 400;

AMX      = AMX(1:N,:);
AMXCASE  = AMXCASE(1:N);
AMXCTRL  = AMXCTRL(1:N);
AMXUSNP  = AMXUSNP(1:N);
AMX.VID  = (1:size(AMX,1))';

disp(AMX(1:9,:))

fprintf('\n %.0f final loci count \n\n',size(AMX,1))


% SET MAIN FISHP TO TRAINING GROUP FISHP
AMX.FISHP      = AMX.TRFISHP;
AMX.FISHOR     = AMX.TRFISHOR;
AMX.CASEREF    = AMX.TRCASEREF;
AMX.CASEALT    = AMX.TRCASEALT;
AMX.CTRLREF    = AMX.TRCTRLREF;
AMX.CTRLALT    = AMX.TRCTRLALT;


clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
AMX AMXCASE AMXCTRL AMXUSNP AMXTRCASE AMXTRCTRL AMXTECASE AMXTECTRL





%--------------------------------------------------------------------------
% REMOVE TOP 100 VARIANTS (SPECIFIC TEST)
% AMX(1:100,:)   = [];
% AMXCASE(1:100) = [];
% AMXCTRL(1:100) = [];
% AMXUSNP(1:100) = [];
% AMX.VID        = (1:size(AMX,1))';








%--------------------------------------------------------------------------
% VISUALIZE CORRELATION BETWEEN TRAINING AND TESTING FISHER'S P-VALUES

clc; close all;
fh1 = figure('Units','normalized','OuterPosition',[.1 .05 .7 .92],'Color','w','MenuBar','none');
ax1 = axes('Position',[.06 .06 .9 .9],'Color','none');


ph1 = loglog(-log(AMX.TRFISHP(1:end))+1 , -log(AMX.TEFISHP(1:end))+1,'-mo',...
    'LineStyle','none',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.49 1 .63],...
    'MarkerSize',10);
grid on; axis tight


disp(AMX(1:50,:))
[RHO,PVAL] = corr(-log(AMX.TRFISHP),-log(AMX.TEFISHP));
disp(RHO);
disp(PVAL)

pause(2); close all;


clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
AMX AMXCASE AMXCTRL AMXUSNP AMXTRCASE AMXTRCTRL AMXTECASE AMXTECTRL






%--------------------------------------------------------------------------
% MAKE  RECTANGLE  NEURAL NET  VARIANT MATRIX
clc;


AMXTRCACO = [AMXTRCASE; AMXTRCTRL];
AMXTECACO = [AMXTECASE; AMXTECTRL];


% [TRNN, TRCA, TRCO] = varmx(AMX,AMXCASE,AMXCTRL,AMXTRCACO);
% [TENN, TECA, TECO] = varmx(AMX,AMXCASE,AMXCTRL,AMXTECACO);


[TRNN, TRCA, TRCO] = unnet(AMX,AMXCASE,AMXCTRL,AMXUSNP,AMXTRCACO);
[TENN, TECA, TECO] = unnet(AMX,AMXCASE,AMXCTRL,AMXUSNP,AMXTECACO);



clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
AMX AMXCASE AMXCTRL AMXUSNP AMXTRCASE AMXTRCTRL AMXTECASE AMXTECTRL TRNN TENN





%--------------------------------------------------------------------------
% MAKE  TRAINING AND TESTING SETS

% TRAIN DATA MATRIX needs to have:
%   variables (variants) in rows
%   examples  (people)   in columns
%
% TRAIN LABEL MATRIX will have as many rows as there are classes. 
% Here there are just two classes, CASE and CTRL. If a person is a CASE,
% put a '1' in the first row, and '0' in all other rows; if a person is
% a CTRL put a '1' in the second row, and '0' everywhere else. This label
% matrix will have as many columns as there are people; and will have
% the same number of columns as the training data matrix.
%
% For example if there are 1000 variants, and 300 cases and 200 ctrls,
% the training data matrix should be 1000x500 and the training label
% matrix should be 2x500.
% 
% TEST DATA MATRIX should have the same number of rows (variants)
% as the TRAIN DATA MATRIX, and those rows should be in the same
% exact order in terms of variable (variant) they represent. The
% columns do not need to be the same length, since those represent
% people, and are treated as independent examples.


[TRAINX,TRAINL,TESTX,TESTL,SRR] = traintestmx(TRNN, TENN);

TRAINMX  = TRAINX';
TRAINLAB = TRAINL';

TESTMX  = TESTX';
TESTLAB = TESTL';



clc
disp('Each of these pairs should match...')
disp([sum(TRNN(2:end,4:20)); sum(TRAINX(:,1:17))]')
disp('Each of these pairs should match...')
disp([sum(TENN(2:end,4:20)); sum(TESTX(:,1:17))]')



clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
AMX AMXCASE AMXCTRL AMXUSNP AMXTRCASE AMXTRCTRL AMXTECASE AMXTECTRL...
TRNN TENN TRAINMX TRAINLAB TESTMX TESTLAB TRAINX TESTX TRAINL TESTL SRR







%--------------------------------------------------------------------------
% EXPORT CSV / TXT FILES FOR WORK IN TENSORFLOW
%{
%% MAKE BETTER INPUTS TO ALLOW FOR NONLINEAR DYNAMICS

% % INSTEAD OF VARIANTS BEING {-1,0,1,2} MAKE THEM {-1,0,1,4}
% 
% TRAINMX(TRAINMX==2) = 4;
% TESTMX(TESTMX==2)   = 4;


%######################################################################
LVQNET()
an LVQ (learning vector quantization) neural networks consist of two
layers. The first layer maps input vectors into clusters that are found
by the network during training. The second layer merges groups of first
layer clusters into the classes defined by the target data.

The total number of first layer clusters is determined by the number of
hidden neurons. The larger the hidden layer the more clusters the first
layer can learn, and the more complex mapping of input to target classes
can be made. The relative number of first layer clusters assigned to
each target class are determined according to the distribution of target
classes at the time of network initialization. This occurs when the
network is automatically configured the first time train is called, or
manually configured with the function configure, or manually initialized
with the function init is called.

lvqnet(hiddenSize,lvqLR,lvqLF) takes these arguments,

hiddenSize 		Size of hidden layer (default = 10) 

lvqLR 			LVQ learning rate (default = 0.01) 

lvqLF 			LVQ learning function (default = 'learnlv1') and returns an LVQ neural network.

The other option for the lvq learning function is learnlv2.




%######################################################################
NETWORK()
network creates new custom networks. It is used to create networks that are then 
customized by functions such as feedforwardnet and narxnet.

net = network without arguments returns a new neural network with no inputs, 
layers or outputs.

net = network(numInputs,numLayers,biasConnect,inputConnect,layerConnect,outputConnect) 
takes these optional arguments (shown with default values):

numInputs		Number of inputs, 0
numLayers		Number of layers, 0
biasConnect		numLayers-by-1 Boolean vector, zeros
inputConnect	numLayers-by-numInputs Boolean matrix, zeros
layerConnect	numLayers-by-numLayers Boolean matrix, zeros
outputConnect	1-by-numLayers Boolean vector, zeros

%######################################################################
PATTERNNET()

[x,t] = iris_dataset;
net = patternnet(10);
net = train(net,x,t);
view(net)
y = net(x);
perf = perform(net,t,y);
classes = vec2ind(y);



%######################################################################
COMPETLAYER()

inputs = iris_dataset;
net = competlayer(2);
net = train(net,inputs);
view(net)
outputs = net(inputs);
classes = vec2ind(outputs);

%######################################################################
SELFORGMAP()

x = simplecluster_dataset;
net = selforgmap([8 8]);
net = train(net,x);
view(net)
y = net(x);
classes = vec2ind(y);



%######################################################################
WORKBENCH

% net = feedforwardnet(5);
% net = train(net,TRAINMX,TRAINLAB);
% nnstart; view(net); net.numLayers = 4;

netP = patternnet([55 25]);
netC = competlayer(2);
netS = selforgmap([2 2]);


% NEURALNET = [netA; netB];
NEURALNET = netP;
view(NEURALNET)


net = train(NEURALNET,TRAINMX,TRAINLAB);


[NNPerf] = netstats(net,TRAINMX,TRAINLAB,TESTMX,TESTLAB,.05,.85,.15);



net = train(netS,TRAINMX);



% save('NNready.mat')
% load('NNready.mat')
%}






% v = ['g4v' num2str(size(TRAINMX,1)) '_'];
% writetable(AMX,  [v 'LOCI.csv'])
% PHENTRAIN = [TRCASE; TRCTRL];
% PHENTEST  = [TECASE; TECTRL];
% writetable(PHENTRAIN  ,[v 'PHENTRAIN.csv'])
% writetable(PHENTEST   ,[v 'PHENTEST.csv'])
% writetable(array2table(TRAINMX)   ,[v 'TRAINMX.csv'])
% writetable(array2table(TRAINLAB)  ,[v 'TRAINLAB.csv'])
% writetable(array2table(TESTMX)    ,[v 'TESTMX.csv'])
% writetable(array2table(TESTLAB)   ,[v 'TESTLAB.csv'])
% 
% % dlmwrite('TRAINMX.csv',TRAINMX,'delimiter',',','precision',1)
% % dlmwrite('TRAINLAB.csv',TRAINLAB,'delimiter',',','precision',1)
% % dlmwrite('TESTMX.csv',TESTMX,'delimiter',',','precision',1)
% % dlmwrite('TESTLAB.csv',TESTLAB,'delimiter',',','precision',1)





return
% #######################################################################
%%         MACHINE LEARNING USING ARTIFICIAL NEURAL NETWORKS
% #######################################################################
% TRAIN THE NEURAL NETWORK
clc;
% save('ReadyForNN.mat')
% load('ReadyForNN.mat')





NN = patternnet([300 200 100]);

net = train(NN,TRAINMX,TRAINLAB);

[NNPerf] = netstats(net,TRAINMX,TRAINLAB,TESTMX,TESTLAB,.05,.85,.15,1);



clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
AMX AMXCASE AMXCTRL AMXUSNP AMXTRCASE AMXTRCTRL AMXTECASE AMXTECTRL...
TRNN TENN TRAINMX TRAINLAB TESTMX TESTLAB TRAINX TESTX...
NN net NNPerf







%% RUN THE TRAINING A BUNCH OF TIMES (DUAL HIDDEN LAYER)


nLoops = 5;

netz={};
NETi = zeros(nLoops,6);
for nn = 1:nLoops
    netz{nn} = patternnet([100 100 100]);
    netz{nn} = train(netz{nn},TRAINMX,TRAINLAB);

    [NNPerf] = netstats(netz{nn},TRAINMX,TRAINLAB,TESTMX,TESTLAB,.05,.85,.15,0);
    NETi(nn,:) = NNPerf;


disp(' ')
disp('    All TRAIN     ALL TEST    MID CORRECT    MID POP   HIGH CORRECT   HIGH POP')
disp(NETi)
disp(' ')
end


FAIL = NETi(:,6)<2;

NETi( FAIL , :) = [];
netz(FAIL)      = [];




clc
disp('Dual-hidden layer L1:500n L2:300n L3:200')
disp(' ')
disp('    All TRAIN     ALL TEST    MID CORRECT    MID POP   HIGH CORRECT   HIGH POP')
disp(NETi)
disp(' ')
disp('Mean:'); disp(mean(NETi))
disp('ALL VALUES ABOVE ARE PERCENTAGES')




clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
AMX AMXCASE AMXCTRL AMXUSNP AMXTRCASE AMXTRCTRL AMXTECASE AMXTECTRL...
TRNN TENN TRAINMX TRAINLAB TESTMX TESTLAB TRAINX TESTX...
NN net NNPerf netz NETi







%% EXAMINE & REPORT ACCURACY FOR HIGH CONFIDENCE DATA
%----------------------------------------------------------------------

NS = NETi ./ 10;


BMOscore = NS(:,5) .* (NS(:,6) ./3) + NS(:,3) + NS(:,2) + NS(:,5);
[~,topBMO] = max(BMOscore);
[~,i] = sort(BMOscore,'descend');

clc; disp('  ConfThresh    Pct_Correct    Pct_Register')
for mm = 1:10

net = netz{1,i(mm)};


TESTLOGIT = net(TESTMX);

TESTCONF  = abs(TESTLOGIT(1,:) - .5);


TEclass = vec2ind(TESTLOGIT);
TEGUESS = round(TESTLOGIT);

TECASEn = sum(TESTLAB(1,:));
TECTRLn = sum(TESTLAB(2,:));
TOTn    = size(TESTLAB,2);

HITS = sum(TEGUESS == TESTLAB)>0;


ConfThresh = 0:.05:.45;

NNC.ConfThresh   = zeros(size(ConfThresh));
NNC.N_hits       = zeros(size(ConfThresh));
NNC.N_tot        = zeros(size(ConfThresh));
NNC.Pct_Correct  = zeros(size(ConfThresh));
NNC.Pct_Register = zeros(size(ConfThresh));



for nn = 1:numel(ConfThresh)

    NNC.ConfThresh(nn) = ConfThresh(nn);

    NNC.N_hits(nn) = sum(HITS(TESTCONF>ConfThresh(nn)));

    NNC.N_tot(nn) = sum(TESTCONF>ConfThresh(nn));

    NNC.Pct_Correct(nn) = (NNC.N_hits(nn) / NNC.N_tot(nn)) * 100;

    NNC.Pct_Register(nn) = (NNC.N_tot(nn) / TOTn) * 100;


end

NETj = [NNC.ConfThresh' NNC.Pct_Correct' NNC.Pct_Register'];


disp(' ');disp(' ');
% disp('  ConfThresh    Pct_Correct    Pct_Register')
disp(NETj)


end



clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
AMX AMXCASE AMXCTRL AMXUSNP AMXTRCASE AMXTRCTRL AMXTECASE AMXTECTRL...
TRNN TENN TRAINMX TRAINLAB TESTMX TESTLAB TRAINX TESTX...
NN net NNPerf netz NNC NETi NETj


% TESTMX
% TESTLAB
% [NNPerf] = netperformance(net,TESTMX,TESTLAB,.05,.85,.15);


















return
%################################################################
%%  ASSESS VARIANT PROFILE OF HIGH-CONFIDENCE POPULATION
%################################################################


NN = patternnet([500 40 30]);

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
fh1 = figure('Units','normalized','OuterPosition',[.01 .05 .95 .92],'Color','w','MenuBar','none');
ax1 = axes('Position',[.06 .05 .43 .90],'Color','none');
ax2 = axes('Position',[.55 .05 .43 .90],'Color','none');


axes(ax1)
imagesc(HICASE)
title('High Confidence Cases')
ax1.YTick = 1:size(HICASE,1);
YLabs = cellstr([repmat('\fontsize{9} ',length(AMX.GENE),1), char(AMX.GENE)]);
ax1.YTickLabel = cellstr(YLabs);

axes(ax2)
imagesc(HICTRL)
title('High Confidence Controls')
ax2.YTick = 1:size(HICTRL,1);
ax2.YTickLabel = cellstr(YLabs);


map =  [.1  .1  .6 ;
        .2  .2  .8 ;
        .95  .4  .2 ;
        1   1   .6];
colormap(map)

% C = jet();
% colormap(C(1:end-6,:))
% colorbar





%%

HICASE = TRAINMX(:,HICONFCASE);
HICTRL = TRAINMX(:,HICONFCTRL);

HICASE = mean(HICASE,2);
HICTRL = mean(HICTRL,2);



close all
fh1 = figure('Units','normalized','OuterPosition',[.01 .05 .95 .92],'Color','w','MenuBar','none');
ax1 = axes('Position',[.06 .05 .25 .90],'Color','none');
ax2 = axes('Position',[.33 .05 .27 .90],'Color','none');
ax3 = axes('Position',[.66 .05 .25 .90],'Color','none');


axes(ax1)
imagesc(HICASE)
title('Average High Confidence Case')
ax1.YTick = 1:size(HICASE,1);
YLabs = cellstr([repmat('\fontsize{9} ',length(AMX.GENE),1), char(AMX.GENE)])
ax1.YTickLabel = cellstr(YLabs);
ax1.XTickLabel = [];

axes(ax2)
imagesc(HICTRL)
title('Average High Confidence Control')
ax2.YTick = 1:size(HICTRL,1);
ax2.YTickLabel = [];
ax2.XTickLabel = [];

C = jet();
colormap(C(1:end-6,:))
colorbar


axes(ax3)
imagesc( abs(HICASE-HICTRL))
title('abs |CASE mean - CTRL mean|')
ax3.YTick = 1:size(HICTRL,1);
ax3.YTickLabel = cellstr(YLabs);
ax3.XTickLabel = [];


C = jet();
colormap(ax3, C(12:end-4,:))
colorbar



%%
clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
AMX AMXCASE AMXCTRL AMXUSNP AMXTRCASE AMXTRCTRL AMXTECASE AMXTECTRL...
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
fh1 = figure('Units','normalized','OuterPosition',[.01 .06 .9 .7],'Color','w','MenuBar','none');
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
AMX AMXCASE AMXCTRL AMXUSNP AMXTRCASE AMXTRCTRL AMXTECASE AMXTECTRL...
TRNN TENN TRAINMX TRAINLAB TESTMX TESTLAB TRAINX TESTX...
NN net NNPerf HC




%################################################################
%%  PERFORM MACHINE LEARNING USING HOMEBREW METHODS
%################################################################
%{

% LOGISTIC REGRESSION NEURAL NET WITH GRADIENT DECENT

Niters = 20;
RES = zeros(Niters,5);

% szTx = size(TRAINX,1);
% nP = fliplr(round(linspace(100,szTx,Niters)))';



for nn = 1:Niters


clearvars -except ADSP GENB LOCI CASE CTRL PHEN RES nn...
TRCASE TRCTRL TECASE TECTRL TRCASEN TRCTRLN TECASEN TECTRLN...
AMX AMXCASE AMXCTRL AMXTRCASE AMXTRCTRL AMXTECASE AMXTECTRL...
TRAINMX TRAINLAB TESTMX TESTLAB TRNN TENN nP


% The *TRNN* matrix is formatted as follows
%
% ROW-1: AMX.VID of variant for TRNN(1,4:end)
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
AMX AMXCASE AMXCTRL AMXTRCASE AMXTRCTRL AMXTECASE AMXTECTRL...
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
AMX AMXCASE AMXCTRL AMXTRCASE AMXTRCTRL AMXTECASE AMXTECTRL...
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
AMX AMXCASE AMXCTRL AMXTRCASE AMXTRCTRL AMXTECASE AMXTECTRL...
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

AMX.NNactiv  = h2(:,2);
AMX.NNguess  = p-1;
NNEYE = AMX(:,[7 26 27]);
disp(NNEYE(1:10,:))

%}
