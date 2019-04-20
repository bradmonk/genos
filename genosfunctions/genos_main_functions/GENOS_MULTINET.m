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







%###############################################################
%%  USE PER-PERSON VARIANT COUNTS TO IDENTIFY OUTLIERS
%###############################################################


Qlo = quantile(PHETRCASE.TOTvars,.01);
Qhi = quantile(PHETRCASE.TOTvars,.99);
PTRCASE  = PHETRCASE(((PHETRCASE.TOTvars > Qlo) & (PHETRCASE.TOTvars < Qhi)),:);


Qlo = quantile(PHETRCTRL.TOTvars,.01);
Qhi = quantile(PHETRCTRL.TOTvars,.99);
PTRCTRL  = PHETRCTRL(((PHETRCTRL.TOTvars > Qlo) & (PHETRCTRL.TOTvars < Qhi)),:);


Qlo = quantile(PHETECASE.TOTvars,.01);
Qhi = quantile(PHETECASE.TOTvars,.99);
PTECASE  = PHETECASE(((PHETECASE.TOTvars > Qlo) & (PHETECASE.TOTvars < Qhi)),:);


Qlo = quantile(PHETECTRL.TOTvars,.01);
Qhi = quantile(PHETECTRL.TOTvars,.99);
PTECTRL  = PHETECTRL(((PHETECTRL.TOTvars > Qlo) & (PHETECTRL.TOTvars < Qhi)),:);



clc; close all;
subplot(2,2,1), histogram(PTRCASE.TOTvars); title('TRAIN CASE')
subplot(2,2,2), histogram(PTRCTRL.TOTvars); title('TRAIN CTRL')
subplot(2,2,3), histogram(PTECASE.TOTvars); title('TEST  CASE')
subplot(2,2,4), histogram(PTECTRL.TOTvars); title('TEST  CTRL')
disp('------------------------------------')
disp('TOTAL VARIANTS PER-PERSON')
disp('                MIN    MAX')
fprintf('TRAINING CASE: %.0f  %.0f \n',min(PTRCASE.TOTvars),  max(PTRCASE.TOTvars))
fprintf('TRAINING CTRL: %.0f  %.0f \n',min(PTRCTRL.TOTvars),  max(PTRCTRL.TOTvars))
fprintf('TESTING  CASE: %.0f  %.0f \n',min(PTECASE.TOTvars),  max(PTECASE.TOTvars))
fprintf('TESTING  CTRL: %.0f  %.0f \n',min(PTECTRL.TOTvars),  max(PTECTRL.TOTvars))
disp('------------------------------------')

clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
PTRCASE PTRCTRL PTECASE PTECTRL








%###############################################################
%%          COUNT NUMBER OF VARIANTS PER LOCI
%###############################################################

USNP = ADSP.UNSNP;


[TRCASEN, TRCTRLN] = varsum(CASE, PTRCASE.SRR, CTRL, PTRCTRL.SRR);

[TECASEN, TECTRLN] = varsum(CASE, PTECASE.SRR, CTRL, PTECTRL.SRR);


[TRCASEUN, TRCTRLUN, TECASEUN, TECTRLUN] = uncsum(USNP,...
    PTRCASE.SRR, PTRCTRL.SRR,PTECASE.SRR, PTECTRL.SRR);



% SAVE COUNTS AS NEW TABLE COLUMNS
LOCI.TRCASEREF = (numel(PTRCASE.SRR)*2) - (TRCASEUN.*2) - TRCASEN;
LOCI.TRCTRLREF = (numel(PTRCTRL.SRR)*2) - (TRCTRLUN.*2) - TRCTRLN;
LOCI.TRCASEALT = TRCASEN;
LOCI.TRCTRLALT = TRCTRLN;


LOCI.TECASEREF = (numel(PTECASE.SRR)*2) - (TECASEUN.*2) - TECASEN;
LOCI.TECTRLREF = (numel(PTECTRL.SRR)*2) - (TECTRLUN.*2) - TECTRLN;
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
PTRCASE PTRCTRL PTECASE PTECTRL







%% ONLY KEEP LOCI WHERE MORE THAN 5 PEOPLE HAVE THE ALT ALLELE


P = (LOCI.TRCASEALT > 5) &...
    (LOCI.TRCTRLALT > 5);

LOCI  = LOCI(P,:);
CASE  = CASE(P);
CTRL  = CTRL(P);
USNP  = USNP(P);



clearvars -except ADSP LOCI CASE CTRL PHEN USNP PTRCASE PTRCTRL PTECASE PTECTRL


%% ONLY KEEP VARIANTS WHERE ALT < REF*1.5

P = (LOCI.TRCASEALT < (LOCI.TRCASEALT.*1.5)) &...
    (LOCI.TRCTRLALT < (LOCI.TRCTRLALT.*1.5));

LOCI  = LOCI(P,:);
CASE  = CASE(P);
CTRL  = CTRL(P);
USNP  = USNP(P);

clearvars -except ADSP LOCI CASE CTRL PHEN USNP PTRCASE PTRCTRL PTECASE PTECTRL










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
pause(3); close all;


clearvars -except ADSP LOCI CASE CTRL PHEN USNP PTRCASE PTRCTRL PTECASE PTECTRL











%%

PTRCASE.TRTECACO = repmat("TRAINCASE",size(PTRCASE,1),1);
PTRCTRL.TRTECACO = repmat("TRAINCTRL",size(PTRCTRL,1),1);
PTECASE.TRTECACO = repmat("TESTCASE", size(PTECASE,1),1);
PTECTRL.TRTECACO = repmat("TESTCTRL", size(PTECTRL,1),1);

ADSP.PHE = [PTRCASE; PTRCTRL; PTECASE; PTECTRL];



% ########################################################################
%% STORE VARIABLES FOR NEURAL NETWORK TRAINING AS 'VLOCI'
% ########################################################################

clearvars -except ADSP LOCI CASE CTRL PHEN USNP


VLOCI     = LOCI;
VCASE     = CASE;
VCTRL     = CTRL;
VUSNP     = USNP;
VPHE      = ADSP.PHE;


clearvars -except ADSP LOCI CASE CTRL PHEN USNP VLOCI VCASE VCTRL VUSNP VPHE





%% TAKE THE TOP N GENES FOR NEURAL NET CLASSIFIER TRAINING



% SET MAIN FISHP TO TRAINING GROUP FISHP
VLOCI.FISHP      = VLOCI.TRFISHP;
VLOCI.FISHOR     = VLOCI.TRFISHOR;
VLOCI.CASEREF    = VLOCI.TRCASEREF;
VLOCI.CASEALT    = VLOCI.TRCASEALT;
VLOCI.CTRLREF    = VLOCI.TRCTRLREF;
VLOCI.CTRLALT    = VLOCI.TRCTRLALT;



% SORT BY TRAINING GROUP FISHP
[~,i]  = sort(VLOCI.TRFISHP);
VLOCI  = VLOCI(i,:);
VCASE  = VCASE(i);
VCTRL  = VCTRL(i);
VUSNP  = VUSNP(i);




% EXTRACT TOP-N NUMBER OF VARIANTS
N = 600;

VLOCI  = VLOCI(1:N,:);
VCASE  = VCASE(1:N);
VCTRL  = VCTRL(1:N);
VUSNP  = VUSNP(1:N);




clc; disp(VLOCI(1:9,:))
fprintf('\n LOCI COUNT FOR FEATURE MATRIX: %.0f \n\n',size(VLOCI,1))

clearvars -except ADSP LOCI CASE CTRL PHEN USNP VLOCI VCASE VCTRL VUSNP VPHE






%###################################################################
%%      MAKE  RECTANGLE  NEURAL NET  VARIANT MATRIX
%###################################################################

AMXTRCACO = VPHE(strcmp(VPHE.TRTECACO,'TRAINCASE')|strcmp(VPHE.TRTECACO,'TRAINCTRL'),:);
AMXTECACO = VPHE(strcmp(VPHE.TRTECACO,'TESTCASE')|strcmp(VPHE.TRTECACO,'TESTCTRL'),:);


doQ = 1;
[TRNN, ~, ~] = nnmake(VLOCI,VCASE,VCTRL,VUSNP,AMXTRCACO,doQ);
clc; disp(TRNN(1:7,1:10))

[TENN, ~, ~] = nnmake(VLOCI,VCASE,VCTRL,VUSNP,AMXTECACO,doQ);
disp(TENN(1:7,1:10))


[TRAINX,TRAINL,TESTX,TESTL] = nnmatrix(TRNN, TENN);

TRMX = TRAINX';
TRLX = TRAINL';

TEMX = TESTX';
TELX = TESTL';


disp(' ')
disp('Each of these pairs should match...')
disp([sum(TRNN(2:end,4:6)); sum(TRAINX(:,1:3))]')
disp('Each of these pairs should match...')
disp([sum(TENN(2:end,4:6)); sum(TESTX(:,1:3))]')



clearvars -except ADSP LOCI CASE CTRL PHEN USNP VLOCI VCASE VCTRL VUSNP VPHE...
TRMX TRLX TEMX TELX











return
% ######################################################################
%%            GET THINGS READY FOR AGE-BASED TRAINING
% ######################################################################


PHE = VPHE;


% ADD YEARS TO ALL CTRL AGES
PHE.AGE(PHE.AD==0) = PHE.AGE(PHE.AD==0) + 16;


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
histogram(PHE.AGE(PHE.AD==1))
hold on
histogram(PHE.AGE(PHE.AD==0))
disp(PHE(1:9,:)); pause(1)




clc; clearvars -except ADSP LOCI CASE CTRL PHEN USNP VLOCI VCASE VCTRL VUSNP...
VPHE PHE TRMX TRLX TEMX TELX






%% MAKE NEURAL NET MATRICES

TRPHE = PHE(strcmp(PHE.TRTECACO,'TRAINCASE')|strcmp(PHE.TRTECACO,'TRAINCTRL'),:);
TEPHE = PHE(strcmp(PHE.TRTECACO,'TESTCASE') |strcmp(PHE.TRTECACO,'TESTCTRL'),:);



% SCRAMBLE TRAINING PHENOTYPE ORDER
N      = size(TRPHE,1);  % Total number of people
k      = randperm(N)';   % Get N random ints in range 1:N
TRPHE  = TRPHE(k,:);     % Scramble Phenotype table

% SCRAMBLE TESTING PHENOTYPE ORDER
N      = size(TEPHE,1);  % Total number of people
k      = randperm(N)';   % Get N random ints in range 1:N
TEPHE  = TEPHE(k,:);     % Scramble Phenotype table





% MAKE THE NEURAL NET TRAINING MATRIX
[TR, vMX] = makeAgeNN(VLOCI,VCASE,VCTRL,VUSNP,TRPHE,0);
disp(TR(1:10,1:10))

TRAIN = TR(2:end,:);
% TRCHR = TRMX(1,5:end);



% MAKE THE NEURAL NET TESTING MATRIX
[TE, vMX] = makeAgeNN(VLOCI,VCASE,VCTRL,VUSNP,TEPHE,0);
disp(TE(1:10,1:10))

TEST  = TE(2:end,:);
% TECHR = TEMX(1,5:end);



% CREATE FEATURE MATRIX AND LABEL MATRIX
TRX   = TRAIN(:,5:end)';
TRL   = TRAIN(:,2)';

TEX   = TEST(:,5:end)';
TEL   = TEST(:,2)';





disp(TRX(1:7,1:7))
disp(' ')
disp(TRL(1,1:7))

clearvars -except ADSP LOCI CASE CTRL PHEN USNP VLOCI VCASE VCTRL VUSNP...
VPHE PHE TRMX TRLX TEMX TELX TRX TRL TEX TEL






%% MAKE LABEL MATRICES HAVE N COLUMNS EQUAL TO N CLASSES


u = unique(TRL);
TRLx = zeros(numel(u),size(TRL,2));
for i=1:size(TRL,2)
    TRLx(:,i) = (u==TRL(i))';
end
TRL = TRLx>0;


TELx = zeros(numel(u),size(TEL,2));
for i=1:size(TEL,2)
    TELx(:,i) = (u==TEL(i))';
end
TEL = TELx>0;



clc; clearvars -except ADSP LOCI CASE CTRL PHEN USNP VLOCI VCASE VCTRL VUSNP...
VPHE PHE TRMX TRLX TEMX TELX TRX TRL TEX TEL









% #######################################################################
%%         MACHINE LEARNING USING ARTIFICIAL NEURAL NETWORKS
% #######################################################################
clc; close all;



% VNN = patternnet([300 100 50]);
% 
% varnet = train(VNN,TRMX,TRLX);



% ANN = patternnet([300 100 50]);
% 
% agenet = train(ANN,TRX,TRL);


%% CHECK PERFORMANCE
clc

[NNPerf] = multinetstats(varnet,TRMX,TRLX,TEMX,TELX,.05,.85,.15,1,...
                         agenet,TRX,TRL,TEX,TEL);


clearvars -except ADSP LOCI CASE CTRL PHEN USNP VLOCI VCASE VCTRL VUSNP...
VPHE PHE TRMX TRLX TEMX TELX TRX TRL TEX TEL agenet varnet



TLAB  =  vec2ind(TELX);
TGUS  =  vec2ind(varnet(TEMX));
disp('PCT CORRECT CASE/CONTROL TEST-SET WITH VARNET (CHANCE=50%):'); 
mean(TLAB==TGUS)


TLAB  =  vec2ind(TEL);
TGUS  =  vec2ind(agenet(TEX));
disp('PCT CORRECT AGE TEST-SET WITH AGENET (CHANCE=25%):'); 
mean(TLAB==TGUS)


% GUESS = round(outputs);
% errors = gsubtract(TESTL,outputs);
% performance = perform(net,TESTL,outputs)





%% CHECK PERFORMANCE OF THE AGENET ON TRAINING DATA

TLAB  = vec2ind(TRL);

TGUS = vec2ind(agenet(TRX));


u = unique(PHE.AGE);

for i = 1:numel(u);  TLAB(TLAB==i) = u(i);  end

for i = 1:numel(u);  TGUS(TGUS==i) = u(i);  end


MU = zeros(size(u));
for i = 1:numel(u);  MU(i) = mean(TGUS(TLAB==u(i)));  end


AGEjit = TLAB + ((rand(numel(TLAB),1)-.5).*8)';
GESjit = TGUS + ((rand(numel(TGUS),1)-.5).*8)';



%---------------------------
close all; figure('Units','normalized','Position',[.05 .04 .6 .8],'Color','w');
ax1 = axes('Position',[.13 .13 .80 .80],'Color','none');


scatter(GESjit, AGEjit, 50, TLAB, '.');

%---------------------------
map =  [.0 .7 .1;  .9 .5 .1;  .9 .1 .6;  .1 .5 .6];
colormap(map)
xlabel({'\fontsize{22} AGE ', '\fontsize{16} (GUESS BIN)'});
ylabel({'\fontsize{22} AGE ', '\fontsize{16} (ACTUAL BIN)'});
ax1.XTick = unique(vec2ind(TRL));
ax1.YTick = unique(vec2ind(TRL));
ax1.FontSize = 16;
cb1 = colorbar;
cb1.Ticks = 1:numel(unique(TLAB));
%---------------------------


% PLOT AVERAGE GUESS
hold on;
ph2=plot(MU,u,'o-.');
ph2.MarkerSize=35;
ph2.LineWidth=4;
ph2.Color=[.01 .5 .5];
ph2.MarkerEdgeColor='none';
ph2.MarkerFaceColor=[.01 .5 .5];


clc; clearvars -except ADSP LOCI CASE CTRL PHEN USNP VLOCI VCASE VCTRL VUSNP...
VPHE PHE TRMX TRLX TEMX TELX TRX TRL TEX TEL agenet varnet




%% CHECK PERFORMANCE OF THE AGENET ON TEST DATA

TLAB = vec2ind(TEL);

TGUS = vec2ind(agenet(TEX));


u = unique(PHE.AGE);

for i = 1:numel(u);  TLAB(TLAB==i) = u(i);  end

for i = 1:numel(u);  TGUS(TGUS==i) = u(i);  end


MU = zeros(size(u));
for i = 1:numel(u);  MU(i) = mean(TGUS(TLAB==u(i)));  end



AGEjit = TLAB + ((rand(numel(TLAB),1)-.5).*8)';
GESjit = TGUS + ((rand(numel(TGUS),1)-.5).*8)';



%---------------------------
close all; figure('Units','normalized','Position',[.05 .04 .6 .8],'Color','w');
ax1 = axes('Position',[.13 .13 .80 .80],'Color','none');


scatter(GESjit, AGEjit, 50, TLAB, '.');

%---------------------------
map =  [.0 .7 .1;  .9 .5 .1;  .9 .1 .6;  .1 .5 .6];
colormap(map)
xlabel({'\fontsize{22} AGE ', '\fontsize{16} (GUESS BIN)'});
ylabel({'\fontsize{22} AGE ', '\fontsize{16} (ACTUAL BIN)'});
ax1.XTick = unique(vec2ind(TRL));
ax1.YTick = unique(vec2ind(TRL));
ax1.FontSize = 16;
cb1 = colorbar;
cb1.Ticks = 1:numel(unique(TLAB));
%---------------------------



% PLOT AVERAGE GUESS
hold on;
ph2=plot(MU,u,'o-.');
ph2.MarkerSize=35;
ph2.LineWidth=4;
ph2.Color=[.01 .5 .5];
ph2.MarkerEdgeColor='none';
ph2.MarkerFaceColor=[.01 .5 .5];


clc; clearvars -except ADSP LOCI CASE CTRL PHEN USNP VLOCI VCASE VCTRL VUSNP...
VPHE PHE TRMX TRLX TEMX TELX TRX TRL TEX TEL agenet varnet




%%
clc;

VLAB  =  vec2ind(TRLX);
VGUS  =  vec2ind(varnet(TRMX));
disp('TRAINING PCT CORRECT CASE/CONTROL WITH VARNET:'); 
disp( mean(VLAB==VGUS) )


ALAB  =  vec2ind(TRL);
AGUS  =  vec2ind(agenet(TRX));
disp('TRAINING PCT CORRECT AGE WITH AGENET:'); 
disp( mean(ALAB==AGUS) )



AVLAB  =  ALAB;
AVLAB(AVLAB==4) = -2;
AVLAB(AVLAB>0)  = 1;
AVLAB = abs(AVLAB);

AVGUS  =  AGUS;
AVGUS(AVGUS==4) = -2;
AVGUS(AVGUS>0)  = 1;
AVGUS = abs(AVGUS);

disp('TRAINING PCT CORRECT CASE/CTRL WITH AGENET:'); 
disp( mean(AVLAB==AVGUS) )



u = unique(PHE.AGE);

for i = 1:numel(u);  ALAB(ALAB==i) = u(i);  end

for i = 1:numel(u);  AGUS(AGUS==i) = u(i);  end


MU = zeros(size(u));
for i = 1:numel(u);  MU(i) = mean(AGUS(ALAB==u(i)));  end

ADJ = SIG(  (MU - mean(MU)) ./ std(MU)  );
ADJ = (ADJ - ADJ(3)) .* -1;

disp(' '); 
disp('SUGGESTED AGENET-OUTPUT ADJUSTMENTS ON VARNET-OUTPUT'); 
disp(ADJ)




clearvars -except ADSP LOCI CASE CTRL PHEN USNP VLOCI VCASE VCTRL VUSNP...
VPHE PHE TRMX TRLX TEMX TELX TRX TRL TEX TEL agenet varnet ADJ



%%
clc;

VLAB  =  vec2ind(TELX);
VGUS  =  vec2ind(varnet(TEMX));
disp('TEST PCT CORRECT CASE/CONTROL WITH VARNET:'); 
disp( mean(VLAB==VGUS) )


ALAB  =  vec2ind(TEL);
AGUS  =  vec2ind(agenet(TEX));
disp('TEST PCT CORRECT AGE WITH AGENET:'); 
disp( mean(ALAB==AGUS) )



AVLAB  =  ALAB;
AVLAB(AVLAB==4) = -2;
AVLAB(AVLAB>0)  = 1;
AVLAB = abs(AVLAB);

AVGUS  =  AGUS;
AVGUS(AVGUS==4) = -2;
AVGUS(AVGUS>0)  = 1;
AVGUS = abs(AVGUS);

disp('TEST PCT CORRECT CASE/CTRL WITH AGENET:'); 
disp( mean(AVLAB==AVGUS) )



clearvars -except ADSP LOCI CASE CTRL PHEN USNP VLOCI VCASE VCTRL VUSNP...
VPHE PHE TRMX TRLX TEMX TELX TRX TRL TEX TEL agenet varnet ADJ
















%% MAKE ADJUSTMENTS TO GUESS WEIGHTS AND EVALUATE PERFORMANCE
clc;


VGUS  =  varnet(TEMX);
VGUS  =  VGUS(1,:);


AGUS  =  mean(agenet(TEMX) .* [2;1;.5;-2]);


VLAB  =  TELX(1,:);

VAGUS  = VGUS + AGUS;

VR = VAGUS;
VR(VR<.5)=0;
VR(VR>=.5)=1;

z=NaN(size(VLAB));

P = [VGUS; AGUS; z;z; VAGUS; z;z; VR; VLAB];
disp(P)

mean(P(8,:)==P(9,:))





%% MAKE ADJUSTMENTS TO GUESS WEIGHTS AND EVALUATE PERFORMANCE
clc;
clearvars -except ADSP LOCI CASE CTRL PHEN USNP VLOCI VCASE VCTRL VUSNP...
VPHE PHE TRMX TRLX TEMX TELX TRX TRL TEX TEL agenet varnet ADJ





VGUS  =  varnet(TEMX);
VGUS  =  VGUS(1,:);

VLAB  =  TELX(1,:);

AGU = agenet(TEMX);


nLoops = 100000;
BetaMax = [2;1;.5;-2];
AlphaMax = 0;
x = rand(nLoops,1).*5;
y = rand(nLoops,1).*2;
k = x-y;
z = rand(nLoops,1).*-3;

tic
for nn = 1:nLoops

    %x = rand*5;
    %y = rand*2;
    %z = rand*-3;
    %Beta = [x; x-y; y; z];

    Beta = [x(nn); k(nn); y(nn); z(nn)];

    AGUS  =  mean(AGU .* Beta);
    %AGUS  =  SIG(mean(AGU .* Beta))-.5;

    VR  = VGUS + AGUS;
    

    VR(VR<.5)=0;
    VR(VR>=.5)=1;

    Alpha = mean(VR==VLAB);

    if Alpha > AlphaMax

        BetaMax = Beta;
        AlphaMax = Alpha;

    end

end
toc

DeltaMax = BetaMax;

disp(BetaMax)
disp(' ')
disp(AlphaMax)


clearvars -except ADSP LOCI CASE CTRL PHEN USNP VLOCI VCASE VCTRL VUSNP...
VPHE PHE TRMX TRLX TEMX TELX TRX TRL TEX TEL agenet varnet...
DeltaMax LambdaMax BetaWeights VGHIC





%% MAKE ADJUSTMENTS TO HIGH CONF GUESS WEIGHTS AND EVALUATE PERFORMANCE
clc;
clearvars -except ADSP LOCI CASE CTRL PHEN USNP VLOCI VCASE VCTRL VUSNP...
VPHE PHE TRMX TRLX TEMX TELX TRX TRL TEX TEL agenet varnet...
BetaMax LambdaMax DeltaMax BetaWeights VGHIC





VGUS  =  varnet(TEMX);
VGUS  =  VGUS(1,:);

VGHIC = VGUS>.85 | VGUS<.15;

VLAB  =  TELX(1,:);

AGU = agenet(TEMX);

VGUS = VGUS(VGHIC);
VLAB = VLAB(VGHIC);
AGU  = AGU(:,VGHIC);





nLoops = 200000;
BetaMax = [2;1;.5;-2];
AlphaMax = 0;
x = rand(nLoops,1).*5;
y = rand(nLoops,1).*2;
k = x-y;
z = rand(nLoops,1).*-3;

tic
for nn = 1:nLoops

    %x = rand*5;
    %y = rand*2;
    %z = rand*-3;
    %Beta = [x; x-y; y; z];

    Beta = [x(nn); k(nn); y(nn); z(nn)];

    AGUS  =  mean(AGU .* Beta);
    %AGUS  =  SIG(mean(AGU .* Beta))-.5;

    VR  = VGUS + AGUS;

    VR(VR<.5)=0;
    VR(VR>=.5)=1;

    Alpha = mean(VR==VLAB);

    if Alpha > AlphaMax

        BetaMax = Beta;
        AlphaMax = Alpha;

    end

end
toc

LambdaMax = BetaMax;

disp(LambdaMax)
disp(' ')
disp(AlphaMax)


clearvars -except ADSP LOCI CASE CTRL PHEN USNP VLOCI VCASE VCTRL VUSNP...
VPHE PHE TRMX TRLX TEMX TELX TRX TRL TEX TEL agenet varnet...
DeltaMax LambdaMax BetaWeights VGHIC


%%
clc;
BetaWeights = repmat(VGHIC,4,1) .* 0.0;

xx = repmat(LambdaMax,1,size(BetaWeights(:,VGHIC),2));
yy = repmat(DeltaMax,1,size(BetaWeights(:,~VGHIC),2));

BetaWeights(:,VGHIC) = xx;

BetaWeights(:,~VGHIC) = yy;


clearvars -except ADSP LOCI CASE CTRL PHEN USNP VLOCI VCASE VCTRL VUSNP...
VPHE PHE TRMX TRLX TEMX TELX TRX TRL TEX TEL agenet varnet...
DeltaMax LambdaMax BetaWeights VGHIC




VGUS  =  varnet(TEMX);
VGUS  =  VGUS(1,:);

VLAB  =  TELX(1,:);

AGU   = agenet(TEMX);

AGUS  =  mean(AGU .* BetaWeights);

VR  = VGUS + AGUS;

VR(VR<.5)=0;
VR(VR>=.5)=1;

PCTCORRECT = mean(VR==VLAB)















%%

% GUESS = round(outputs);
% errors = gsubtract(TESTL,outputs);
% performance = perform(net,TESTL,outputs)



% varout = varnet(TEMX);
% 
% ageout = agenet(TEMX);
% 
% TEGUESS = round(varout(1,:));
% 
% TECORR = TELX(1,:);
% 
% mean(TEGUESS == TECORR)
% 
% ageadjust = (SIG(vec2ind(ageout)- (0:4)' )-.5) ./ 2   ;
% 
% mean(  round(varout(1,:)) == TECORR)
% 
% mean(  round(ageadjust(1,:)+varout(1,:)) == TECORR)
% mean(  round(ageadjust(2,:)+varout(1,:)) == TECORR)
% mean(  round(ageadjust(3,:)+varout(1,:)) == TECORR)
% mean(  round(ageadjust(4,:)+varout(1,:)) == TECORR)


% COLORBY = findgroups(TRPHE.APOE);
% COLORBY = findgroups(TRPHE.CONSORTIUM);
% COLORBY = findgroups(TRPHE.Consent);
% COLORBY = findgroups(TRPHE.AD);
% COLORBY = findgroups(TRPHE.COHORTNUM);
% COLORBY = findgroups(TRPHE.RACE);
% COLORBY = findgroups(TRPHE.ETHNIC);
% COLORBY = findgroups(TRPHE.RACE,TRPHE.ETHNIC);
% COLORBY = findgroups(TRPHE.GOODCOH);
% COLORBY = findgroups(TRPHE.SEX);
%---------
% TRPHE.PREVAD(isnan(TRPHE.PREVAD)) = 0;
% TRPHE.INCAD(isnan(TRPHE.INCAD))   = 0;
% [COLORBY,~,~] = findgroups(TRPHE.PREVAD,TRPHE.INCAD);
% TRPHE.INCAD(COLORBY==4) = 0;
% [COLORBY,v,d] = findgroups(TRPHE.PREVAD,TRPHE.INCAD);
%---------
