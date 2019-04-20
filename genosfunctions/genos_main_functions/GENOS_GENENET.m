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
PHEN = ADSP.PHEN;


clearvars -except ADSP LOCI CASE CTRL PHEN
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
% HERE THE CODE MERGES THESE INTO A SINGLE TABLE FOR:
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

% PHETRCTRL(PHETRCTRL.AGE < 75 , :) = [];
% PHETECTRL(PHETECTRL.AGE < 75 , :) = [];







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

clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
TRCASE TRCTRL TECASE TECTRL








%###############################################################
%%          COUNT NUMBER OF VARIANTS PER LOCI
%###############################################################
%
% The varsum() function will go through each known variant loci
% and check whether anyone's SRR ID from your subset of IDs match
% all known SRR IDs for that loci. It will then sum the total
% number of alleles (+1 for hetzy-alt, +2 for homzy-alt) for each
% loci and return the totals.

USNP = ADSP.UNSNP;


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







%% REMOVE VARIANTS WITH VERY LOW COUNTS ALT ALLELE COUNTS

PASS = (LOCI.TRCASEALT > 5) & (LOCI.TRCTRLALT > 5);

LOCI  = LOCI(PASS,:);
CASE  = CASE(PASS);
CTRL  = CTRL(PASS);
USNP  = USNP(PASS);


clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
TRCASE TRCTRL TECASE TECTRL




%% REMOVE VARIANTS WHERE ALT > REF

PASS = (LOCI.TRCASEREF > (LOCI.TRCASEALT./1.5)) &...
       (LOCI.TRCTRLREF > (LOCI.TRCTRLALT./1.5));


LOCI  = LOCI(PASS,:);
CASE  = CASE(PASS);
CTRL  = CTRL(PASS);
USNP  = USNP(PASS);


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


% [~,i] = sort(LOCI.TRFISHP);
% LOCI  = LOCI(i,:);
% CASE  = CASE(i);
% CTRL  = CTRL(i);
% USNP  = USNP(i);
% LOCI.VID = (1:size(LOCI,1))';


clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
TRCASE TRCTRL TECASE TECTRL








%########################################################################
%% STORE VARIABLES FOR NEURAL NETWORK TRAINING AS 'VLOCI'
%########################################################################


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




%% REMOVE NON-ASYMETRIC VARIANTS


VLOCI.CASEREF = VLOCI.TRCASEREF;
VLOCI.CASEALT = VLOCI.TRCASEALT;
VLOCI.CTRLREF = VLOCI.TRCTRLREF;
VLOCI.CTRLALT = VLOCI.TRCTRLALT;
VLOCI.FISHP   = VLOCI.TRFISHP;
VLOCI.FISHOR  = VLOCI.TRFISHOR;



PASS = (VLOCI.FISHP < .01) |...
       (VLOCI.FISHOR > 2)  |...
       (VLOCI.FISHOR < .5);

VLOCI  = VLOCI(PASS,:);
VCASE  = VCASE(PASS);
VCTRL  = VCTRL(PASS);
VUSNP  = VUSNP(PASS);





% SORT VARIANT TABLE BY FISHP
[~,i] = sortrows(VLOCI.FISHP);
VLOCI = VLOCI(i,:);
VCASE = VCASE(i);
VCTRL = VCTRL(i);
VUSNP = VUSNP(i);




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
pause(3); close all; disp(' ')

clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL






%###################################################################
%%      MAKE  RECTANGLE  NEURAL NET  VARIANT MATRIX
%###################################################################

clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL


VTRCACO = [VTRCASE; VTRCTRL];
VTECACO = [VTECASE; VTECTRL];


doQ = 0;
[TRNN, ~, ~] = genennmake(VLOCI,VCASE,VCTRL,VUSNP,VTRCACO,doQ);
disp(TRNN(1:20,1:10))

[TENN, ~, ~] = genennmake(VLOCI,VCASE,VCTRL,VUSNP,VTECACO,doQ);
disp(TENN(1:20,1:10))


clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRNN TENN





Gi = TRNN(1,4:end)';
VX = TRNN(2:end,4:end);
[Ui , Ug] = unique(Gi);
TRNX = zeros(size(VX,1), size(Ui,1));
for i=1:numel(Ui)
    TRNX(:,i) = sum(  VX(: , Gi == Ui(i) ) > 0   , 2);
end

disp(TRNX(1:10,1:10))
imagesc(TRNX)

ADNN = padarray(TRNX,[1 3],0,'pre');
ADNN(: , 1)  =  TRNN(:,1);    % COL1: SRR
ADNN(: , 2)  =  TRNN(:,2);    % COL2: CACO
ADNN(: , 3)  =  TRNN(:,3);    % COL3: BIAS
ADNN(1 , 4:end)  =  Ui';    % ROW1: GENEi

TRNN = ADNN;


clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRNN TENN


Gi = TENN(1,4:end)';
VX = TENN(2:end,4:end);
[Ui , Ug] = unique(Gi);
TENX = zeros(size(VX,1), size(Ui,1));
for i=1:numel(Ui)
    TENX(:,i) = sum(  VX(: , Gi == Ui(i) ) > 0   , 2);
end

ADNN = padarray(TENX,[1 3],0,'pre');
ADNN(: , 1)  =  TENN(:,1);    % COL1: SRR
ADNN(: , 2)  =  TENN(:,2);    % COL2: CACO
ADNN(: , 3)  =  TENN(:,3);    % COL3: BIAS
ADNN(1 , 4:end)  =  Ui';    % ROW1: GENEi

TENN = ADNN;


clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRNN TENN


disp(' ');disp(' ');
disp(TRNN(1:10,1:10))
disp(' ');disp(' ');
disp(TENN(1:10,1:10))


imagesc(TRNN(2:end,4:end))
pause(3); close all;



%% PREP MATRICES FOR NEURAL NET TRAINING

[TRAINX,TRAINL,TESTX,TESTL] = nnmatrix(TRNN, TENN);

TRAINMX  = TRAINX';
TRAINLAB = TRAINL';

TESTMX  = TESTX';
TESTLAB = TESTL';





disp(' ')
disp('Each of these pairs should match...')
disp([sum(TRNN(2:end,4:6)); sum(TRAINX(:,1:3))]')
disp('Each of these pairs should match...')
disp([sum(TENN(2:end,4:6)); sum(TESTX(:,1:3))]')



clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRNN TENN...
TRAINMX TRAINLAB TESTMX TESTLAB TRAINX TESTX


return
% #######################################################################
%%         MACHINE LEARNING USING ARTIFICIAL NEURAL NETWORKS
% #######################################################################



NN = patternnet([300 200 100]);

net = train(NN,TRAINMX,TRAINLAB);

[NNPerf] = netstats(net,TRAINMX,TRAINLAB,TESTMX,TESTLAB,.05,.85,.15,1);



clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL...
TRNN TENN TRAINMX TRAINLAB TESTMX TESTLAB TRAINX TESTX...
NN net NNPerf







%% RUN THE TRAINING A BUNCH OF TIMES (DUAL HIDDEN LAYER)


nLoops = 10;

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
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL...
TRNN TENN TRAINMX TRAINLAB TESTMX TESTLAB TRAINX TESTX...
NN net NNPerf netz NETi







%% EXAMINE & REPORT ACCURACY FOR HIGH CONFIDENCE DATA
%----------------------------------------------------------------------

NS = NETi ./ size(netz,2);


BMOscore = NS(:,5) .* (NS(:,6) ./3) + NS(:,3) + NS(:,2) + NS(:,5);
[~,topBMO] = max(BMOscore);
[~,i] = sort(BMOscore,'descend');

clc; disp('  ConfThresh    Pct_Correct    Pct_Register')
for mm = 1:size(netz,2)

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
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL...
TRNN TENN TRAINMX TRAINLAB TESTMX TESTLAB TRAINX TESTX...
NN net NNPerf netz NNC NETi NETj


% TESTMX
% TESTLAB
% [NNPerf] = netperformance(net,TESTMX,TESTLAB,.05,.85,.15);












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

