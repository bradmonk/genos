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


clearvars -except ADSP LOCI CASE CTRL USNP PHEN
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
pause(2)

clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
TRCASE TRCTRL TECASE TECTRL






%% DISCARD VARIANT LOCI WITH KNOWN ULTRA-RARE ALT ALLELE COUNTS

PASS = (LOCI.CASEALT > 5) | (LOCI.CTRLALT > 5);

LOCI  = LOCI(PASS,:);
CASE  = CASE(PASS);
CTRL  = CTRL(PASS);
USNP  = USNP(PASS);






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



%% CREATE HISTOGRAMS OF VARIANT COUNTS

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
disp('Calculating Fishers Exact Test Statistics for Training Group...')
[FISHP, FISHOR] = fishp_mex(LOCI.TRCASEREF,LOCI.TRCASEALT,...
                            LOCI.TRCTRLREF,LOCI.TRCTRLALT);

LOCI.TRFISHP  = FISHP;
LOCI.TRFISHOR = FISHOR;




% COMPUTE FISHERS STATISTICS FOR THE TRAINING GROUP
disp('Calculating Fishers Exact Test Statistics for Testing Group...')
[FISHP, FISHOR] = fishp_mex(LOCI.TECASEREF, LOCI.TECASEALT,...
                            LOCI.TECTRLREF, LOCI.TECTRLALT);

LOCI.TEFISHP  = FISHP;
LOCI.TEFISHOR = FISHOR;





disp('done.'); close all
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
N = 500;

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








%% READ IN ANDREI-GENERATED CSV TABLE
%{
clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL


f=['/Users/bradleymonk/Documents/MATLAB/GIT/genomics/genos/'...
   'genospython/pythondata/train20180403092606_40_1.csv'];

TX = readtable(f);



% GET CASE CONTROL LABELS
txcc = TX.Var42;
TX.Var42 = [];
txcc = string(txcc);
txcc = str2double(txcc);


% GET SRR LISTS
txsrr = TX.Var1(2:end,:);

TRPHE = [TRCASE; TRCTRL];
TRSRR = TRPHE.SRR;

GOODPHE = ADSP.PHEN;
GOODPHE = GOODPHE(GOODPHE.GOODCOH,:);
GOODSRR = GOODPHE.SRR;
numel(GOODSRR)



% GET CHR:POS
txmat = table2array(TX);
txpos = uint64(txmat(1,2:end)');
CHRPOS = VLOCI.CHRPOS;




% BRAD STUFF PREFIXED WITH 'b'; ANDREI STUFF PREFIXED WITH 'a'


% CHECK HOW MANY bSRR IN aSRR
[i,j] = ismember(TRSRR,txsrr);
sum(i)
numel(TRSRR)
numel(txsrr)


% CHECK HOW MANY aSRR ARE FROM ANY MEMBER OF GOOD COHORTS
[i,j] = ismember(txsrr,GOODSRR);
sum(i)
numel(txsrr)
numel(GOODSRR)


% CHECK NUMBER OF aLOCI FOUND IN bLOCI
[i,j] = ismember(txpos,CHRPOS);
sum(i)
numel(txpos)
numel(CHRPOS)





% CREATE PHE MATRIX USING THE SRR VALUES ANDREI USED
[i,j] = ismember(GOODPHE.SRR,txsrr);
sum(i)
numel(GOODPHE.SRR)
numel(txsrr)


PHEX = GOODPHE(i,:);



[i,j] = ismember(ADSP.LOCI.CHRPOS,txpos);
sum(i)
numel(txpos)
numel(GOODSRR)

LOX = ADSP.LOCI(i,:);
CAX = ADSP.CASE(i);
COX = ADSP.CTRL(i);
UNX = ADSP.UNSNP(i);



[FULLTRX, TRX, TRL] = makennet(LOX,CAX,COX,UNX,PHEX,0);
FULLTRX(:,2:4) = [];
BX = FULLTRX;
BX = sortrows(BX);
BX = sortrows(BX')';
BX = uint64(BX);

a = txmat';
a(isnan(a)) = 0;
AX = a';
AX = sortrows(AX);
AX = sortrows(AX')';
AX = uint64(AX);

clc



AX = uint64(AX);
BX = uint64(BX);

disp(' ')
disp( AX(1:5,1:5) )
disp( BX(1:5,1:5) )
disp(' ')
disp( mean(mean(AX == BX)) )


%}



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
%%  ASSIGN NEURAL NET ACTIVATION VALUE TO PARTICIPANT PHE
%################################################################

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



clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE TR TE...
TRX TRL TEX TEL net netz NETi ACTR ACTE




%% ADD THOSE ACTIVATIONS TO PHE TABLES AND WRITE TO CSV FILE


TRPHE.ACTR = ACTR;
TEPHE.ACTE = ACTE;

TRPHE.MUACTR = mean(ACTR,2);
TEPHE.MUACTE = mean(ACTE,2);



BRAAK = TRPHE.BRAAK;

BRAAK(TRPHE.AUTOPSY==0) = NaN;



clc;

nanmean(BRAAK(TRPHE.MUACTR > .90))

nanmean(BRAAK( ((TRPHE.MUACTR<.90)&(TRPHE.MUACTR>.80))  ))

nanmean(BRAAK( ((TRPHE.MUACTR<.80)&(TRPHE.MUACTR>.70))  ))

nanmean(BRAAK( ((TRPHE.MUACTR<.70)&(TRPHE.MUACTR>.60))  ))

nanmean(BRAAK( ((TRPHE.MUACTR<.60)&(TRPHE.MUACTR>.50))  ))

nanmean(BRAAK(TRPHE.MUACTR < .5))



% writetable(TRPHE,'TRPHE.csv')
% writetable(TEPHE,'TEPHE.csv')
% save('GOODBRAAK.mat')


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