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
%     PHEN: [10910�22 table]
%     LOCI: [94483�24 table]
%     CASE: {94483�1 cell}
%     CTRL: {94483�1 cell}
%     USNP: {94483�1 cell}
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

close all; clear; clc; rng('shuffle');

genosdir = fileparts(which('GENOS_RM1.m'));
cd(genosdir);
subfuncpath = [genosdir filesep 'genosfunctions'];
datasetpath = [genosdir filesep 'genosdata'];
gpath = [genosdir pathsep subfuncpath pathsep datasetpath];
addpath(gpath)

diary on;
difi = ['BRADLAB-' char(datetime('now','Format','yyy-MM-dd'))];
diary(difi)
get(0,'Diary')



    which('GENOSDATA.mat')
    disp('Loading GENOSDATA.mat please wait...')
ADSP = load('GENOSDATA.mat');
    disp('GENOSDATA.mat loaded.')






clearvars -except ADSP




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


clc; clearvars -except ADSP PHEN LOCI CASE CTRL USNP
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
% 

clc; clearvars -except ADSP PHEN LOCI CASE CTRL USNP




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
% COHSET(COHSET.APOE == 33,:) = [];
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

clc; clearvars -except ADSP PHEN LOCI CASE CTRL USNP COHSET


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

clc; clearvars -except ADSP PHEN LOCI CASE CTRL USNP COHSET



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




clearvars -except ADSP PHEN LOCI CASE CTRL USNP COHSET...
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





clearvars -except ADSP PHEN LOCI CASE CTRL USNP COHSET...
PHETRCASE PHETRCTRL PHETECASE PHETECTRL








%==========================================================================
%%   REMOVE CASE & CTRL PARTICIPANTS BASED ON AGE
%==========================================================================
clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
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




clearvars -except ADSP PHEN LOCI CASE CTRL USNP COHSET...
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

clearvars -except ADSP PHEN LOCI CASE CTRL USNP COHSET...
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




clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
TRCASE TRCTRL TECASE TECTRL





%==========================================================================
%% DISCARD VARIANT LOCI WITH KNOWN ULTRA-RARE ALT ALLELE COUNTS
%==========================================================================

clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
TRCASE TRCTRL TECASE TECTRL




PASS = (LOCI.CASEALT > 20) | (LOCI.CTRLALT > 20);

LOCI  = LOCI(PASS,:);
CASE  = CASE(PASS);
CTRL  = CTRL(PASS);
USNP  = USNP(PASS);


clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
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


clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
TRCASE TRCTRL TECASE TECTRL








%==========================================================================
%%               COMPUTE FISHER'S P-VALUE
%==========================================================================
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



%==========================================================================
%==========================================================================
%==========================================================================
%%
%
% PREPARE DATA FOR NEURAL NET CLASSIFIER SUPERVISED LEARNING
%
%==========================================================================
%==========================================================================
%==========================================================================
% LOAD THE DATASET
%{
close all; clear; clc; rng('shuffle');

genosdir = fileparts(which('GENOS_RM1.m'));
cd(genosdir);
subfuncpath = [genosdir filesep 'genosfunctions'];
datasetpath = [genosdir filesep 'genosdata'];
gpath = [genosdir pathsep subfuncpath pathsep datasetpath];
addpath(gpath)

diary on;
difi = ['BRADLAB-' char(datetime('now','Format','yyy-MM-dd'))];
diary(difi)
get(0,'Diary')



    which('GENOSDATAFISH.mat')
    disp('Loading GENOSDATA.mat please wait...')
ADSP = load('GENOS_COH-BRAAK_APOE-ALL.mat','LOCI','CASE','CTRL','PHEN','USNP','TRCASE','TRCTRL','TECASE','TECTRL')
    disp('GENOSDATAFISH.mat loaded.')


% save('GOODCOH_BRAAKCOH_APOEALL.mat','LOCI','CASE','CTRL','PHEN','USNP','TRCASE','TRCTRL','TECASE','TECTRL')



clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
TRCASE TRCTRL TECASE TECTRL
%}


%==========================================================================
%% STORE VARIABLES FOR NEURAL NETWORK TRAINING AS 'VLOCI'
%==========================================================================
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
% [~,i]  = sort(VLOCI.TRFISHP);
[~,i]  = sort(VLOCI.CHRPOS);
VLOCI  = VLOCI(i,:);
VCASE  = VCASE(i);
VCTRL  = VCTRL(i);
VUSNP  = VUSNP(i);


clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL





%==========================================================================
% REMOVE SELECT GENES BY NAME
%==========================================================================
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






%==========================================================================
%% LOOP OVER 1:N NUMBER OF VARIANTS
%==========================================================================
clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL


TRMEAN = zeros(1,200);
TRPOPU = zeros(1,200);
TRHIMU = zeros(1,200);
HOMEAN = zeros(1,200);
HOPOPU = zeros(1,200);
HOHIMU = zeros(1,200);



VI = fliplr(1:200);

for vi = 1:200
GRPn = 1;

% REMOVE VARS FROM HIGH-TO-LOW P-VALUE
HI2LOW = 1;
SNPn = VI(vi);
SNPi = 1:SNPn;


% REMOVE VARS FROM LOW-TO-HIGH P-VALUE
% HI2LOW = 0;
% SNPn = vi;
% SNPi = SNPn:200;



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

%==========================================================================
% TAKE THE TOP N or P<x GENES FOR NEURAL NET CLASSIFIER TRAINING
%==========================================================================


% SORT VARIANTS BY EITHER TRFISHP|CHRPOS
[~,j]  = sort(VLOCI.TRFISHP);
VLOCI  = VLOCI(j,:);
VCASE  = VCASE(j);
VCTRL  = VCTRL(j);
VUSNP  = VUSNP(j);


% EXTRACT TOP-N NUMBER OF VARIANTS
VLOCI  = VLOCI(SNPi,:);
VCASE  = VCASE(SNPi);
VCTRL  = VCTRL(SNPi);
VUSNP  = VUSNP(SNPi);




%==========================================================================
%      MAKE  RECTANGLE  NEURAL NET  VARIANT MATRIX
%==========================================================================

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



%==========================================================================
%                  LOGISTIC REGRESSION 
%==========================================================================


TXX = VTRX(:,7:end);
HXX = VTEX(:,7:end);
TX = [ones(size(TXX,1),1) TXX]; % ADD AN INTERCEPT COLUMN
HX = [ones(size(HXX,1),1) HXX]; % ADD AN INTERCEPT COLUMN


TL = VTRX(:,2);
HL = VTEX(:,2);




% PERFORM THE SO-CALLED MACHINE LEARNING STEP
BETA = pinv(TX' * TX) * (TX' * TL);

fprintf('\n Solved GLM OLS for %0.f beta coefficients. \n\n',size(BETA,1));  




% GET TRAINED LINEAR MODEL PREDICTIONS (ACTIVATIONS)
TRfx = nansum( TX .* BETA' ,2);
HOfx = nansum( HX .* BETA' ,2);


% DESCRETIZE PREDICTIONS
TRy = round(TRfx);
HOy = round(HOfx);


% GRADE THE PREDICTIONS
TRmu = nanmean(TRy == TL);
HOmu = nanmean(HOy == HL);




% [TRAINED] HIGH CONFIDENCE PREDICTIONS
TRhi = (TRfx>.8) | (TRfx<.2);
TRhiN  = TRy(TRhi) == TL(TRhi);
TRhiMu = nanmean(TRhiN);
TRhiPop = nanmean(TRhi);


% [HOLDOUT] HIGH CONFIDENCE PREDICTIONS
HOhi = (HOfx>.8) | (HOfx<.2);
HOhiN  = HOy(HOhi) == HL(HOhi);
HOhiMu = nanmean(HOhiN);
HOhiPop = nanmean(HOhi);




TRMEAN(vi) = TRmu;
TRHIMU(vi) = TRhiMu;
TRPOPU(vi) = TRhiPop;
HOMEAN(vi) = HOmu;
HOHIMU(vi) = HOhiMu;
HOPOPU(vi) = HOhiPop;



disp('----------  TRAINING SET  ----------')
fprintf('TRAIN correct %0.0f%% overall(pop:100%%)\n' ,(TRmu .* 100))
fprintf('TRAIN correct %0.0f%% hicon  (pop: %0.0f%%)\n\n',TRhiMu.*100,TRhiPop.*100)

disp('----------  TESTING SET  ----------')
fprintf('TEST correct %0.0f%% overall(pop:100%%)\n' ,(HOmu .* 100))
fprintf('TEST correct %0.0f%% hicon  (pop: %0.0f%%)\n\n',HOhiMu.*100,HOhiPop.*100)
pause(.1)


end


% FILL IN NAN VALUES WIHT NANMEAN

TRMEAN(isnan(TRMEAN)) = nanmean(TRMEAN);
TRHIMU(isnan(TRHIMU)) = nanmean(TRHIMU);
TRPOPU(isnan(TRPOPU)) = nanmean(TRPOPU);
HOMEAN(isnan(HOMEAN)) = nanmean(HOMEAN);
HOHIMU(isnan(HOHIMU)) = nanmean(HOHIMU);
HOPOPU(isnan(HOPOPU)) = nanmean(HOPOPU);


%==========================================================================
%% 4-PACK GRAPHS TRAINING & HOLDOUT PROPORTION CORRECT x N-VARIANT LOCI
%  LOOSE AXES LIMITS
%==========================================================================
clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL HI2LOW...
NVARS TRmu TRhiPop HOmu HOhiPop TRMEAN TRHIMU TRPOPU HOMEAN HOHIMU HOPOPU



%################   FOUR PACK   ################
close all
fh01 = figure('Units','normalized','OuterPosition',[.01 .05 .95 .88],'Color','w');
ax01 = axes('Position',[.06 .56 .4 .4],'Color','none');
ax02 = axes('Position',[.56 .56 .4 .4],'Color','none');
ax03 = axes('Position',[.06 .06 .4 .4],'Color','none');
ax04 = axes('Position',[.56 .06 .4 .4],'Color','none');

pk = pink;
pi = pk(15,:);
PAR = flipud(parula); HOT = hot;
PH = [flipud(HOT(end-5:end,:)); PAR(1:end-5,:)];

NVARS = fliplr(1:200);


axes(ax01);
ph01 = scatter(NVARS,TRMEAN,500,'k.'); hold on;
ph01 = scatter(NVARS,TRMEAN,300,[.05 .2 .6],'.');
% ax01.YLim = [.3 1];
% line([0 200],[.5 .5],'Color','k','LineStyle','--','LineWidth',1);
title('TRAINING SUBSET (ALL CONFIDENCE)')
if HI2LOW == 1
xlabel('N loci of top 200 (removed from HIGH-to-LOW P-value)')
else
xlabel('N loci of top 200 (removed from LOW-to-HIGH P-value)')
end
ylabel('Proportion correct')



axes(ax02);
ph02 = scatter(NVARS,HOMEAN,500,'k.'); hold on;
ph02 = scatter(NVARS,HOMEAN,300,[.05 .2 .6],'.');
% ax02.YLim = [.3 1]; 
% line([0 200],[.5 .5],'Color','k','LineStyle','--','LineWidth',1);
title('HOLDOUT SUBSET (ALL CONFIDENCE)')
if HI2LOW == 1
xlabel('N loci of top 200 (removed from HIGH-to-LOW P-value)')
else
xlabel('N loci of top 200 (removed from LOW-to-HIGH P-value)')
end
ylabel('Proportion correct')



axes(ax03);
ph03 = scatter(NVARS,TRHIMU,100,'k.'); hold on;
ph03 = scatter(NVARS,TRHIMU,400,TRPOPU,'.');
% ax03.YLim = [.3 1];
% line([0 200],[.5 .5],'Color','k','LineStyle','--','LineWidth',1);
title('TRAINED SUBSET (HIGH CONFIDENCE)')
if HI2LOW == 1
xlabel('N loci of top 200 (removed from HIGH-to-LOW P-value)')
else
xlabel('N loci of top 200 (removed from LOW-to-HIGH P-value)')
end
ylabel('Proportion correct')
cb03 = colorbar;
cb03.Label.String = 'Proportion of Sample';
cb03.Label.FontSize = 12;
cb03.Label.HorizontalAlignment = 'center';
cb03.Label.VerticalAlignment = 'bottom';
cb03.Label.Rotation = -90;


axes(ax04);
ph04 = scatter(NVARS,HOHIMU,100,'k.'); hold on;
ph04 = scatter(NVARS,HOHIMU,400,HOPOPU,'.');
% ax03.YLim = [.3 1];
% line([0 200],[.5 .5],'Color','k','LineStyle','--','LineWidth',1);
title('HOLDOUT SUBSET (HIGH CONFIDENCE)')
if HI2LOW == 1
xlabel('N loci of top 200 (removed from HIGH-to-LOW P-value)')
else
xlabel('N loci of top 200 (removed from LOW-to-HIGH P-value)')
end
ylabel('Proportion correct')
cb04 = colorbar;
cb04.Label.String = 'Proportion of Sample';
cb04.Label.FontSize = 12;
cb04.Label.HorizontalAlignment = 'center';
cb04.Label.VerticalAlignment = 'bottom';
cb04.Label.Rotation = -90;



colormap(fh01,PH)

%------------------------------------------%
pause(1)
set(gcf, 'PaperPositionMode', 'auto');
dt=char(datetime(datetime,'Format','yyyy-MM-dd-HH-mm-ss'));
saveas(gcf, ['/Users/bradleymonk/Desktop/APOE33_4G_' dt '.png']);
pause(1)
%------------------------------------------%













%==========================================================================
% 2-PACK GRAPHS TRAINING & HOLDOUT PROPORTION CORRECT x N-VARIANT LOCI
%  STRICT AXES LIMITS
%==========================================================================
clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL HI2LOW...
NVARS TRmu TRhiPop HOmu HOhiPop TRMEAN TRHIMU TRPOPU HOMEAN HOHIMU HOPOPU



%################   TWO PACK   ################
close all
fh01 = figure('Units','normalized','OuterPosition',[.01 .06 .8 .7],'Color','w');
ax01 = axes('Position',[.05 .11 .42 .81],'Color','none');
ax02 = axes('Position',[.55 .11 .42 .81],'Color','none');

pk = pink;
pi = pk(15,:);
PAR = flipud(parula); HOT = hot;
PH = [flipud(HOT(end-5:end,:)); PAR(1:end-5,:)];

NVARS = fliplr(1:200);




axes(ax01);
ph01 = scatter(NVARS,TRMEAN,500,'k.'); hold on;
ph01 = scatter(NVARS,TRMEAN,500,[.01 .1 .5],'.'); hold on;
hold on;
ph03 = scatter(NVARS,TRHIMU,100,'k.'); hold on;
ph03 = scatter(NVARS,TRHIMU,400,TRPOPU,'.');
ax01.YLim = [.3 1];
line([0 200],[.5 .5],'Color','k','LineStyle','--','LineWidth',1);
title('TRAINING SUBSET (HIGH & ALL CONFIDENCE)')
if HI2LOW == 1
xlabel('N loci of top 200 (removed from HIGH-to-LOW P-value)')
else
xlabel('N loci of top 200 (removed from LOW-to-HIGH P-value)')
end
ylabel('Proportion correct')
cb01 = colorbar;
cb01.Label.String = 'Proportion of Sample';
cb01.Label.FontSize = 12;
cb01.Label.HorizontalAlignment = 'center';
cb01.Label.VerticalAlignment = 'bottom';
cb01.Label.Rotation = -90;



axes(ax02);
ph02 = scatter(NVARS,HOMEAN,500,'k.'); hold on;
ph02 = scatter(NVARS,HOMEAN,300,[.01 .1 .5],'.');
ax02.YLim = [.3 1];
hold on;
ph04 = scatter(NVARS,HOHIMU,100,'k.'); hold on;
ph04 = scatter(NVARS,HOHIMU,400,HOPOPU,'.');
ax02.YLim = [.3 1];
line([0 200],[.5 .5],'Color','k','LineStyle','--','LineWidth',1);
title('HOLDOUT SUBSET (HIGH & ALL CONFIDENCE)')
if HI2LOW == 1
xlabel('N loci of top 200 (removed from HIGH-to-LOW P-value)')
else
xlabel('N loci of top 200 (removed from LOW-to-HIGH P-value)')
end
ylabel('Proportion correct')
cb02 = colorbar;
cb02.Label.String = 'Proportion of Sample';
cb02.Label.FontSize = 12;
cb02.Label.HorizontalAlignment = 'center';
cb02.Label.VerticalAlignment = 'bottom';
cb02.Label.Rotation = -90;

colormap(fh01,PH)

%------------------------------------------%
pause(1)
set(gcf, 'PaperPositionMode', 'auto');
dt=char(datetime(datetime,'Format','yyyy-MM-dd-HH-mm-ss'));
saveas(gcf, ['/Users/bradleymonk/Desktop/APOE33_2G_' dt '.png']);
pause(1)
%------------------------------------------%






































return
%==========================================================================
%
%%                ARTIFICIAL NEURAL NETWORKS (MATLAB BUILT-IN)
%
%==========================================================================
rng('shuffle');






%==========================================================================
%% LOOP OVER 1:N NUMBER OF VARIANTS
%==========================================================================
clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL HI2LOW...
NVARS TRmu TRhiPop HOmu HOhiPop TRMEAN TRHIMU TRPOPU HOMEAN HOHIMU HOPOPU


TRMEAN = zeros(1,200);
TRPOPU = zeros(1,200);
TRHIMU = zeros(1,200);
HOMEAN = zeros(1,200);
HOPOPU = zeros(1,200);
HOHIMU = zeros(1,200);



VI = fliplr(1:200);

for vi = 1:1
GRPn = 1;

% REMOVE VARS FROM HIGH-TO-LOW P-VALUE
HI2LOW = 1;
SNPn = VI(vi);
SNPi = 1:SNPn;


% REMOVE VARS FROM LOW-TO-HIGH P-VALUE
% HI2LOW = 0;
% SNPn = vi;
% SNPi = SNPn:200;


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

%==========================================================================
% TAKE THE TOP N or P<x GENES FOR NEURAL NET CLASSIFIER TRAINING
%==========================================================================


% SORT VARIANTS BY EITHER TRFISHP|CHRPOS
[~,j]  = sort(VLOCI.TRFISHP);
VLOCI  = VLOCI(j,:);
VCASE  = VCASE(j);
VCTRL  = VCTRL(j);
VUSNP  = VUSNP(j);


% EXTRACT TOP-N NUMBER OF VARIANTS
VLOCI  = VLOCI(SNPi,:);
VCASE  = VCASE(SNPi);
VCTRL  = VCTRL(SNPi);
VUSNP  = VUSNP(SNPi);


%==========================================================================
%      MAKE  RECTANGLE  NEURAL NET  VARIANT MATRIX
%==========================================================================

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
[VTRX, TRX, TRL] = makeomicsnet(VLOCI,VCASE,VCTRL,VUSNP,TRPHE,[-1 0 1 3]);
[VTEX, TEX, TEL] = makeomicsnet(VLOCI,VCASE,VCTRL,VUSNP,TEPHE,[-1 0 1 3]);


%==========================================================================
%                  TRAIN NEURAL NETS
%==========================================================================


TXX = VTRX(:,7:end);
HXX = VTEX(:,7:end);
TX = [ones(size(TXX,1),1) TXX]; % ADD AN INTERCEPT COLUMN
HX = [ones(size(HXX,1),1) HXX]; % ADD AN INTERCEPT COLUMN


TL = VTRX(:,2);
HL = VTEX(:,2);
TL = dummyvar(categorical(VTRX(:,2)==1))'; 
HL = dummyvar(categorical(VTEX(:,2)==1))';


% ESTABLISH PATTERNNET PARAMETERS
NN = patternnet([200 70 40]);
NN.trainFcn = 'trainscg';
NN.trainParam.max_fail = 40;
NN.divideFcn = 'dividerand';




% MACHINE LEARNING: TRAIN THE NEURAL NET CLASSIFIER
net = train(NN, TX , TL );









end































%% EVALUATE TRAINED NEURAL NETWORK PERFORMANCE


TRAINED_CONF = net(TX);
HOLDOUT_CONF = net(HX);

TRAINED_GUESS = vec2ind(TRAINED_CONF);
HOLDOUT_GUESS = vec2ind(HOLDOUT_CONF);


TRAINED_LABELS = vec2ind(TL);
HOLDOUT_LABELS = vec2ind(HL);


TRAINED_PCTALL = mean(TRAINED_GUESS == TRAINED_LABELS);
HOLDOUT_PCTALL = mean(HOLDOUT_GUESS == HOLDOUT_LABELS);


TRC = abs(TRAINED_CONF-.5); TRC = TRC(1,:);
HOC = abs(HOLDOUT_CONF-.5); HOC = HOC(1,:);

TRHI = TRC > .3;
HOHI = HOC > .3;

TRAINED_PHICON = mean(TRAINED_GUESS(TRHI) == TRAINED_LABELS(TRHI));
HOLDOUT_PHICON = mean(HOLDOUT_GUESS(HOHI) == HOLDOUT_LABELS(HOHI));



clc
fprintf('[TRAINED] PERCENT CORRECT OVERALL: '); disp(TRAINED_PCTALL);

fprintf('[HOLDOUT] PERCENT CORRECT OVERALL: '); disp(HOLDOUT_PCTALL);

fprintf('[TRAINED] PERCENT CORRECT HIGHCON: '); disp(TRAINED_PHICON);

fprintf('[HOLDOUT] PERCENT CORRECT HIGHCON: '); disp(HOLDOUT_PHICON);












