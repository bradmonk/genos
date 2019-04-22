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
% genosdir = fileparts(which('GENOS_PCARM.m'));
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
PHEN = ADSP.PHEN;
USNP = ADSP.UNSNP;


clearvars -except ADSP LOCI CASE CTRL PHEN USNP
clc; disp(LOCI(1:9,:))




%% IDENTIFY KEEPER PHEN AND MAKE CASE:CTRL EQUAL SIZE
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





%% PUT HALF THE TEST GROUP BACK INTO THE TRAINING GROUP

szCA = round(size(PHETECASE,1)/2);
szCO = round(size(PHETECTRL,1)/2);

PHETRCASE = [PHETRCASE; PHETECASE(1:szCA,:)];
PHETRCTRL = [PHETRCTRL; PHETECTRL(1:szCO,:)];

PHETECASE(1:szCA,:) = [];
PHETECTRL(1:szCO,:) = [];


clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
PHETRCASE PHETRCTRL PHETECASE PHETECTRL






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






%% MAKE SURE EVERYONE IN CTRL GROUP IS OVER 75 YEARS OLD

PHETRCTRL(PHETRCTRL.AGE < 75 , :) = [];
PHETECTRL(PHETECTRL.AGE < 75 , :) = [];







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

PASS = (LOCI.CASEALT > 3) | (LOCI.CTRLALT > 3);

LOCI  = LOCI(PASS,:);
CASE  = CASE(PASS);
CTRL  = CTRL(PASS);
USNP  = USNP(PASS);




[TRCASEN,TRCTRLN,TECASEN,TECTRLN,TRCASEUN,TRCTRLUN,TECASEUN,TECTRLUN] =...
    snpsum(CASE,CTRL,USNP,TRCASE.SRR,TRCTRL.SRR,TECASE.SRR,TECTRL.SRR);


%###############################################################
%%          COUNT NUMBER OF VARIANTS PER LOCI
%###############################################################

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
TRCASE TRCTRL TECASE TECTRL

return





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
N = 400;

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























% #######################################################################
%%         MACHINE LEARNING USING ARTIFICIAL NEURAL NETWORKS
% #######################################################################



% NN = patternnet([300 200 100]);
% 
% net = train(NN,TRX,TRL);
% 
% [NNPerf] = netstats(net,TRX,TRL,TEX,TEL,.05,.85,.15,1);




clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE...
FULLTRX FULLTEX TRX TRL TEX TEL net netz




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
disp(['         NETn    All TRAIN     ALL TEST     '...
      'MID CORRECT    MID POP   HIGH CORRECT   HIGH POP'])
disp(NETi)
disp(' ')
end


clc
disp('Dual-hidden layer L1:300n L2:100n L3:100')
disp(' ')
disp(['         NETn    All TRAIN     ALL TEST     '...
      'MID CORRECT    MID POP   HIGH CORRECT   HIGH POP'])
disp(NETi)
disp(' ')
disp('Mean:'); disp(mean(NETi))
disp('ALL VALUES ABOVE ARE PERCENTAGES')




clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE TR TE...
TRX TRL TEX TEL net netz NETi






%################################################################
%%  ASSESS VARIANT PROFILE OF HIGH-CONFIDENCE POPULATION
%################################################################


% SET THE CURRENT NEURAL NET TO THE BEST TRAINED NET FROM ABOVE
NS = NETi ./ size(netz,2);
BMOscore = NS(:,6) .* (NS(:,7) ./3) + NS(:,4) + NS(:,3) + NS(:,6);
[~,topBMO] = max(BMOscore);
[~,i] = sort(BMOscore,'descend');
net = netz{1,i(1)};



% GET INDEX OF HIGH CONFIDENCE CORRECT DECISIONS
ACTIVATION = net(TRX);
CLASSES    = vec2ind(ACTIVATION);
ACTIVATION = ACTIVATION(1,:);
GUESSES    = round(ACTIVATION);


iCORRECT = TRL(1,:) == GUESSES;

iHICASE  = ACTIVATION>.9;
iHICTRL  = ACTIVATION<.1;

iCORHICASE = iCORRECT & iHICASE;
iCORHICTRL = iCORRECT & iHICTRL;






clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE TR TE...
TRX TRL TEX TEL net netz NETi iCORHICASE iCORHICTRL









return
%################################################################
%%  PERFORM PCA ON HIGH CONFIDENCE CASE AND CONTROL MATRICES
%################################################################


% CREATE HIGH CONFIDENCE VARIANT-BY-PERSON MATRICES

HICMXCASE = TRX(:,iCORHICASE);

HICMXCTRL = TRX(:,iCORHICTRL);



% SET VALUES BELOW ZERO TO ZERO

HICMXCASE(HICMXCASE<0) = 0;

HICMXCTRL(HICMXCTRL<0) = 0;



% PCA HIGH CONFIDENCE (VARIANT-BY-PERSON)' MATRICES

[CASEcoeff,CASEscore,~,~,~,~] = pca(HICMXCASE');

[CTRLcoeff,CTRLscore,~,~,~,~] = pca(HICMXCTRL');




% GET STANDARD DEVIATION AND SCORE TRANSPOSE

CTRLstdR = std(CTRLscore,1,1);

CTRLscorep = CTRLscore';



% PLOT DATA

close all; figure; 
plot(CTRLstdR); 
title('std score ctrl10pct')





%% GET SORTED SCORES

[scoreSortedB,scoreSortedI]= sort(abs(CTRLscore),1,'descend');


for n=1:size(scoreSortedB,2)

    GtS(:,n)=(scoreSortedB(:,n)>2.5*CTRLstdR(n)); 

end



[r,c]=find(GtS==2);

SumGtS=sum(GtS,1);





%% GET THE GENE LIST FOR THE SORTED SCORES

GENELIST = VLOCI.GENE;

for n=1:size(scoreSortedB,2)

    GeneGropsGS{n} = GENELIST(scoreSortedI(GtS(:,n),n));

end

 




 

%% RUN LOOP TO GET TARGET VARIANTS

observedN = []; expectedN = []; fmean = [];

for i = 1:size(GeneGropsGS,2)

    g = GeneGropsGS{i};

    ai = ismember(GENELIST,g);

    bi = (HICMXCTRL>0);

    d = []; f = [];

    for j = 1:size(HICMXCTRL,1)

        c = bi(j,:) + (ai');

        d(j) = sum(c==2) >= sum(ai);

        f(j) = sum(c==2);

    end


    observedN(i) = sum(d);

    fmean(i) = mean(f(:));

    varmu = mean(bi,1);

    varprod = varmu(ai);

    cp = cumprod(varprod);

    cp = cp(end);

    expectedN(i) = cp .* size(HICMXCTRL,1);

end


outcome1 = observedN ./ expectedN;

outcome2 = fmean ./ expectedN;

close all
histogram(f)
