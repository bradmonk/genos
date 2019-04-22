%% GENOS: 
% 
%--------------------------------------------------------------------------
% 
% SUMMARY TABLE OF THE 24 COHORTS
% 
% COHID    CONSOR    STUDY        COHORT    CASES    CTRLS    TOTAL    %CASE    EQL  BRAAK ID GOOD
% 01       DGC       Adult_Chng    ACT        323      945     1268       25    323    1   01    1
% 02       ADGC      AD_Centers    ADC       2438      817     3255       75    817    0   02    1
% 03       CHARGE    Athrosclro    ASC         39       18       57       68     18    0   03    0
% 04       CHARGE    Aus_Stroke    SKE        121        5      126       96      5    0   04    0
% 05       ADGC      ChiT_Aging    CHA         27      204      231       12     27    0   05    0
% 06       CHARGE    CardioHlth    CHS        250      583      833       30    250    1   06    1
% 07       ADGC      Hispanic_S    HSP        160      171      331       48    160    0   07    0
% 08       CHARGE    Erasmus_Er    ERF         45        0       45      100      0    0   08    0
% 09       CHARGE    Framingham    FHS        157      424      581       27    157    1   09    1
% 10       ADGC      Gene_Diffs    GDF        111       96      207       54     96    1   10    1
% 11       ADGC      NIA_LOAD      LOD        367      109      476       77    109    1   11    1
% 12       ADGC      Aging_Proj    MAP        138      277      415       33    138    1   12    1
% 13       ADGC      Mayo_Clini    MAY        250       99      349       72     99    1   13    1
% 14       ADGC      Miami_Univ    MIA        186       14      200       93     14    1   14    0
% 15       ADGC      AD_Genetic    MIR        316       15      331       95     15    0   15    0
% 16       ADGC      Mayo_cl_PD    MPD          0       20       20        0      0    1   16    0
% 17       ADGC      NationC_AD    NCA        160        0      160      100      0    1   17    0
% 18       ADGC      Wash_Unive    RAS         46        0       46      100      0    1   18    0
% 19       ADGC      Relig_Ordr    ROS        154      197      351       44    154    1   19    1
% 20       CHARGE    RotterdamS    RDS        276      813     1089       25    276    0   20    1
% 21       ADGC      Texas_AD_S    TAR        132       12      144       92     12    0   21    0
% 22       ADGC      Un_Toronto    TOR          9        0        9      100      0    0   22    0
% 23       ADGC      Vanderbilt    VAN        210       26      236       89     26    1   23    1
% 24       ADGC      WashNY_Age    WCA         34      116      150       23     34    0   24    1
% 
% 
% GOODCOHORTS = [1 2         6 7   9 10 11 12 13               19 20     23 24]
% BRAKCOHORTS = [1           6     9 10 11 12 13 14   16 17 18 19        23   ]
%--------------------------------------------------------------------------
% 
% THE ADSP DATASET - WHAT'S BEING IMPORTED?
%
% 
% The dataset that will be loaded in STEP-1 below will import 5 variables
% and store them into a structural array named 'ADSP'. If you type ADSP
% into the command prompt you will see...
% 
% 
% >> ADSP
% 
% ADSP = 
%   struct with fields:
% 
%     PHEN: [10910×22 table]
%     LOCI: [94483×24 table]
%     CASE: {94483×1 cell}
%     CTRL: {94483×1 cell}
%     USNP: {94483×1 cell}
% 
% 
% ...(maybe not in this specific order) the 5 container variables.
% 
% 
% PHEN    a table containing the phenotype information for each
%         participant. If you type head(ADSP.PHEN) you can see what
%         data each column contains.
% 
% 
% LOCI    a table containing genotype info for each exome variant locus.
%         if you type head(ADSP.LOCI) you can see what data each column
%         contains.
% 
% 
% 
% CASE    the last three are cell arrays, containing a list of 
% CTRL    participant IDs & HET/HOM status. Each have 1 cell per row
% USNP    of LOCI (~94483 cells); they are (at least upon import) 
%         pre-sorted in corresponding order, which allows us to
%         iterate over each loci/cell and tally each HET (+1) or 
%         HOM (+2) that matches subsets of participant IDs. (there
%         is an optimized function specifically designed to
%         perform this task, as you will see below).
%         
% 
%--------------------------------------------------------------------------





%==========================================================================
%% STEP-1: LOAD THE DATASET
%==========================================================================
%{.
close all; clear; clc; rng('shuffle');
genosdir = fileparts(which('GENOS.m'));
matlabdir = fileparts(which('MATLABS.m'));
cd(genosdir);


subfuncpath = [genosdir filesep 'genosfunctions'];
datasetpath = [matlabdir filesep 'genosdata'];
gpath = [genosdir pathsep subfuncpath pathsep datasetpath];
addpath(gpath)


which('GENOSDATA.mat')
ADSP = load('GENOSDATA.mat');



disp('dataset loaded')
clearvars -except IJ ADSP

%}


%==========================================================================
%%   START MAIN LOOP
%==========================================================================

for IJ = 1:50


%==========================================================================
%%   CARBON COPY MAIN VARIABLES FROM ADSP.STRUCT
%==========================================================================
%
% After evaluating this section, each variable will be copied from the
% ADSP structural array to their own base variable. This is done so
% that (1) you can access their data directly (e.g. LOCI.GENE(1:5) instead
% of ADSP.LOCI.GENE(1:5) ) and so that (2) you can always restart fresh
% here, by running this segment of code, rather than having to import the
% data from the .mat file in the section above.


LOCI = ADSP.LOCI;
CASE = ADSP.CASE;
CTRL = ADSP.CTRL;
USNP = ADSP.USNP;
PHEN = ADSP.PHEN;


clc; clearvars -except IJ ADSP PHEN LOCI CASE CTRL USNP
head(PHEN)
head(LOCI)





%==========================================================================
%% IDENTIFY KEEPER PHEN AND MAKE CASE:CTRL EQUAL SIZE
%==========================================================================
% 
% 1. Identify cohorts in PHEN with a reasonable number of CASE and CTRL
% 
% 2. Take a random selection of people from those cohorts so that
%    there are the same number of CASE and CTRL.
% 
% 3. The rest of the people from the keeper cohorts can be used to
%    increase the size of the NN test group.
%  

clc; clearvars -except IJ ADSP PHEN LOCI CASE CTRL USNP




% COPY PHEN TABLE TO A VARIABLE CALLED 'COHSET'
COHSET = PHEN;



%-------------------------------------------------------------
%
%   REMOVE UNWANTED COHORTS FROM THE PHENOTYPE TABLE
%
%-------------------------------------------------------------



% THE COHORTS YOU WANT TO USE THIS RUN
GOODCOHORTS = [1 2 6 7 9 10 11 12 13 19 20 23 24];
BRAKCOHORTS = [1 6 9 10 11 12 13 14 16 17 18 19 23];
USECOHORTS  = unique([GOODCOHORTS BRAKCOHORTS]);

COHSET = COHSET(sum(COHSET.COHORTNUM == USECOHORTS , 2)>0,:);




%-------------------------------------------------------------
%
%         **REMOVE** PEOPLE WITH CERTAIN APOE VARIANTS
%
%-------------------------------------------------------------
% COHSET(COHSET.APOE == 44,:) = [];
% COHSET(COHSET.APOE == 24,:) = [];
% COHSET(COHSET.APOE == 34,:) = [];
COHSET(COHSET.APOE == 33,:) = [];
% COHSET(COHSET.APOE == 23,:) = [];
% COHSET(COHSET.APOE == 22,:) = [];










% REMOVE PEOPLE BASED ON OTHER PHENOTYPE QUALITIES???
% COHSET(COHSET.XXX == YYY,:) = [];



% RANDOMLY SHUFFLE TABLE ROWS TO ENSURE UNIQUE SUBSET EACH RUN
COHSET = COHSET(randperm(size(COHSET,1)),:);    


% close all; cohhist(COHSET); fig = gcf
% cohcounts(PHETRCASE,PHETRCTRL,PHETECASE,PHETECTRL)
% fig2plotly();
% response = fig2plotly(gcf, 'filename', 'matlab-grouped-bar', 'strip', false);
% plotly_url = response.url;

clc; clearvars -except IJ ADSP PHEN LOCI CASE CTRL USNP COHSET


% close all; 
% COHORTS = unique(COHSET.COHORTNUM);
% nCA=sum((COHSET.COHORTNUM==COHORTS') & (COHSET.AD==1));
% nCO=sum((COHSET.COHORTNUM==COHORTS') & (COHSET.AD==0));
% COUNTS = [nCA' nCO'];
% fh=figure; 
% ax=axes;
% ph=bar(COUNTS,'group');
% ax.XTick = 1:numel(COHORTS); 
% ax.XTickLabel=string(COHORTS); 
% ax.XTickLabel=num2str(COHORTS); 
% ax.XTickLabel=cellstr(num2str(COHORTS)); 
% xlabel('Cohort'); 
% ylabel('People');
% legend(ax,'Case','Control')
% PLOTLY
% response = fig2plotly(fh, 'filename', 'matlab-grouped-bar', 'strip', false);
% plotly_url = response.url;





%==========================================================================
%%   GENERATE TRAINING/HOLDOUT SUBSETS - COUNTERBALANCE COHORTS
%==========================================================================
% 
% THIS IS SOMEWHAT CONVOLUTED, BUT IT WORKS. PROBABLY NOT
% WORTH TRYING TO UNDERSTAND EACH LINE; IF ANYTHING, JUST CHECK
% TO SEE THAT THE OUTPUTS ARE BALANCED. IN FACT, IF YOU RUN
% THIS SECTION, THE FINAL TALLY WILL PRINT TO THE CONSOLE
% AND YOU WILL SEE FOR EXAMPLE THAT 
% 
%   Tr_CASE_CTRL_PctRow
% 
% SHOULD EACH BE 50:50 CASE:CTRL. IN THE SECTION THAT FOLLOWS
% WE BREAK SOME OF THAT PERFECT SYMMETRY, OTHERWISE THERE
% JUST ISN'T ENOUGH TRAINING DATA...

clc; clearvars -except IJ ADSP PHEN LOCI CASE CTRL USNP COHSET



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




disp(PHETRCASE(1:9,:))
disp('--------------------'); disp('N TRAINING & TESTING EXAMPLES')
disp(' '); fprintf('PHETRCASE... %.0f \n',size(PHETRCASE,1));
disp(' '); fprintf('PHETRCTRL... %.0f \n',size(PHETRCTRL,1));
disp(' '); fprintf('PHETECASE... %.0f \n',size(PHETECASE,1));
disp(' '); fprintf('PHETECTRL... %.0f \n',size(PHETECTRL,1));
disp('--------------------');
cohcounts(PHETRCASE,PHETRCTRL,PHETECASE,PHETECTRL)




clearvars -except IJ ADSP PHEN LOCI CASE CTRL USNP COHSET...
PHETRCASE PHETRCTRL PHETECASE PHETECTRL


%-------------------------------------------------------------
close all;
[~,COHORTS] = cohhist(COHSET);
% cohcounts(TRCASE,TRCTRL,TECASE,TECTRL,2);
% [G,APOE] = findgroups(COHSET.APOE(COHSET.AD==0));
% N = splitapply(@numel,COHSET.APOE(COHSET.AD==0),G);
% disp([APOE,N]);
% close all; fh = gcf;
% f=fig2plotly(fh, 'filename', 'barplot_test')
% web(f.url);




%==========================================================================
%%   PUT SOME OF THE TEST GROUP BACK INTO THE TRAINING GROUP
%==========================================================================

% PUT 1/3RD OF THE REMAINING DATA INTO THE TRAINING SET
% USE THE OTHER 2/3RDS AS THE HOLDOUT DATASET. THE TRAINING
% DATA SHOULD STILL BE REASONABLY BALANCED - DEPENDING ON
% WHO WAS REMOVED ABOVE (CERTAIN APOE SUBSETS MAY ALTER THIS);
% THE AMOUNT TO RETURN TO THE TRAINING SET CAN EASILY BE CHANGED
% BY MANIPULATING THE NEXT TWO LINES TO BE /3 /4 /5 ETC TO RETURN
% 1/3, 1/4, 1/5 RESPECTIVELY.


szCA = round(size(PHETECASE,1)/2);
szCO = round(size(PHETECTRL,1)/2);

PHETRCASE = [PHETRCASE; PHETECASE(1:szCA,:)];
PHETRCTRL = [PHETRCTRL; PHETECTRL(1:szCO,:)];

PHETECASE(1:szCA,:) = [];
PHETECTRL(1:szCO,:) = [];



disp('--------------------'); disp('N TRAINING & TESTING EXAMPLES')
disp(' '); fprintf('PHETRCASE... %.0f \n',size(PHETRCASE,1));
disp(' '); fprintf('PHETRCTRL... %.0f \n',size(PHETRCTRL,1));
disp(' '); fprintf('PHETECASE... %.0f \n',size(PHETECASE,1));
disp(' '); fprintf('PHETECTRL... %.0f \n',size(PHETECTRL,1));
disp('--------------------');
cohcounts(PHETRCASE,PHETRCTRL,PHETECASE,PHETECTRL)





clearvars -except IJ ADSP PHEN LOCI CASE CTRL USNP COHSET...
PHETRCASE PHETRCTRL PHETECASE PHETECTRL








%==========================================================================
%%   REMOVE CASE & CTRL PARTICIPANTS BASED ON AGE
%==========================================================================
clearvars -except IJ ADSP LOCI CASE CTRL PHEN USNP...
PHETRCASE PHETRCTRL PHETECASE PHETECTRL




% REMOVE CTRL PARTICIPANTS YOUNGER THAN...
PHETRCTRL(PHETRCTRL.AGE < 72 , :) = [];
PHETECTRL(PHETECTRL.AGE < 72 , :) = [];



% REMOVE CASE PARTICIPANTS OLDER THAN...
PHETRCASE(PHETRCASE.AGE > 90 , :) = [];
PHETECASE(PHETECASE.AGE > 90 , :) = [];




disp('--------------------'); disp('N TRAINING & TESTING EXAMPLES')
disp(' '); fprintf('PHETRCASE... %.0f \n',size(PHETRCASE,1));
disp(' '); fprintf('PHETRCTRL... %.0f \n',size(PHETRCTRL,1));
disp(' '); fprintf('PHETECASE... %.0f \n',size(PHETECASE,1));
disp(' '); fprintf('PHETECTRL... %.0f \n',size(PHETECTRL,1));
disp('--------------------');
cohcounts(PHETRCASE,PHETRCTRL,PHETECASE,PHETECTRL)




clearvars -except IJ ADSP PHEN LOCI CASE CTRL USNP COHSET...
PHETRCASE PHETRCTRL PHETECASE PHETECTRL






%==========================================================================
%%  USE PER-PERSON VARIANT COUNTS TO IDENTIFY OUTLIERS
%==========================================================================
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

clearvars -except IJ ADSP PHEN LOCI CASE CTRL USNP COHSET...
PHETRCASE PHETRCTRL PHETECASE PHETECTRL


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




clearvars -except IJ ADSP LOCI CASE CTRL PHEN USNP...
TRCASE TRCTRL TECASE TECTRL





%==========================================================================
%% DISCARD VARIANT LOCI WITH KNOWN ULTRA-RARE ALT ALLELE COUNTS
%==========================================================================

clearvars -except IJ ADSP LOCI CASE CTRL PHEN USNP...
TRCASE TRCTRL TECASE TECTRL




PASS = (LOCI.CASEALT > 20) | (LOCI.CTRLALT > 20);

LOCI  = LOCI(PASS,:);
CASE  = CASE(PASS);
CTRL  = CTRL(PASS);
USNP  = USNP(PASS);


clearvars -except IJ ADSP LOCI CASE CTRL PHEN USNP...
TRCASE TRCTRL TECASE TECTRL





%==========================================================================
%%          COUNT NUMBER OF VARIANTS PER LOCI
%==========================================================================
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


clearvars -except IJ ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL



%==========================================================================
%%               COMPUTE FISHER'S P-VALUE
%==========================================================================
clearvars -except IJ ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL
clc; disp('COMPUTING FISHERS EXACT TEST STATISTICS PER DNA LOCUS')

% AAA = LOCI.TRCASEREF(1:100,:);
% AAB = LOCI.TRCASEALT(1:100,:);
% AAC = LOCI.TRCTRLREF(1:100,:);
% AAD = LOCI.TRCTRLALT(1:100,:);
% [FISHP, FISHOR] = fishp(AAA,AAB,AAC,AAD);







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


clearvars -except IJ ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL


%------------------------------------------%
pause(1)
dt=char(datetime(datetime,'Format','yyyy-MM-dd-HH-mm-ss'));
save(['F:\ML\genosdat\APOE2x3x4x_' dt '.mat'],...
    'LOCI','PHEN','TRCASE','TRCTRL','TECASE','TECTRL');
pause(1)
%------------------------------------------%

end















