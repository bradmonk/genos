%% GENOS_FISHP
%{
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
%}
%==========================================================================
%% STEP-1: LOAD THE DATASET
%==========================================================================
clc; close all; clear; rng('shuffle'); f = filesep;
P.root  = [f 'Users' f 'bradleymonk' f 'Documents' f 'MATLAB'];
P.home =  [P.root f 'GIT' f 'genomics' f 'genos'];
P.funs  = [P.home f 'genos_functions'];
P.data  = [P.home f 'genos_data'];
P.figs  = [P.home f 'genos_figures'];
P.mat1  = [P.data f 'APOE'];
addpath(join(string(struct2cell(P)),pathsep,1))
cd(P.home); clearvars -except P




which('GENOSDATAFINAL.mat')
ADSP = load('GENOSDATAFINAL.mat');



ADSP.GOODCOHORTS = [1 2 6 7 9 10 11 12 13 19 20 23 24];
ADSP.BRAKCOHORTS = [1 6 9 10 11 12 13 14 16 17 18 19 23];
%ADSP.USE_COHORT = unique([ADSP.GOODCOHORTS ADSP.BRAKCOHORTS]);
ADSP.USE_COHORT = unique([ADSP.GOODCOHORTS]);

% ADSP.USE_APOE = [22 23 24 33 34 44];
% ADSP.USE_APOT = '22_23_24_33_34_44';
% ADSP.USE_APOE = [22 23 24 34 44];
% ADSP.USE_APOT = '22_23_24_34_44';
% ADSP.USE_APOE = [22 23 34 44];
% ADSP.USE_APOT = '22_23_34_44';
ADSP.USE_APOE = [33];
ADSP.USE_APOT = '33';



clearvars -except P ADSP


%==========================================================================
%%   START MAIN LOOP
%==========================================================================
for IJ = 1:50
%% CARBON COPY
LOCI = ADSP.LOCI;
CASE = ADSP.CASE;
CTRL = ADSP.CTRL;
USNP = ADSP.USNP;
PHEN = ADSP.PHEN;
clc; clearvars -except IJ ADSP PHEN LOCI CASE CTRL USNP


%% CHOOSE COHORTS & APOE SUBSET TO USE IN GROUP GENERATOR

COHSET = ADSP.PHEN;

COHSET = COHSET(sum(COHSET.COHORTNUM == ADSP.USE_COHORT , 2)>0,:);

COHSET = COHSET(sum(COHSET.APOE == ADSP.USE_APOE ,2)>0,:);


% RANDOMLY SHUFFLE TABLE ROWS TO ENSURE UNIQUE SUBSET EACH RUN
COHSET = COHSET(randperm(size(COHSET,1)),:);    

clc; clearvars -except IJ ADSP PHEN LOCI CASE CTRL USNP COHSET


%==========================================================================
%%   GENERATE TRAINING/HOLDOUT SUBSETS - COUNTERBALANCE COHORTS
%==========================================================================
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
%-------------------------------------------------------------
rca={};rco={};
for nn = 1:numel(cohNums)
    rca{nn} = randperm( nCACO(nn,1)  , round(minCACO(nn) * .7));
    rco{nn} = randperm( nCACO(nn,2)  , round(minCACO(nn) * .7));
end


% GET ROW INDEX NUMBER FOR PEOPLE NOT CHOSEN ABOVE
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




% disp(PHETRCASE(1:9,:))
% disp('--------------------'); disp('N TRAINING & TESTING EXAMPLES')
% disp(' '); fprintf('PHETRCASE... %.0f \n',size(PHETRCASE,1));
% disp(' '); fprintf('PHETRCTRL... %.0f \n',size(PHETRCTRL,1));
% disp(' '); fprintf('PHETECASE... %.0f \n',size(PHETECASE,1));
% disp(' '); fprintf('PHETECTRL... %.0f \n',size(PHETECTRL,1));
% disp('--------------------');
% cohcounts(PHETRCASE,PHETRCTRL,PHETECASE,PHETECTRL)




clearvars -except IJ ADSP PHEN LOCI CASE CTRL USNP COHSET...
PHETRCASE PHETRCTRL PHETECASE PHETECTRL

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


clearvars -except IJ ADSP PHEN LOCI CASE CTRL USNP COHSET...
PHETRCASE PHETRCTRL PHETECASE PHETECTRL

%==========================================================================
%%  USE PER-PERSON VARIANT COUNTS TO IDENTIFY OUTLIERS
%==========================================================================
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

clearvars -except IJ ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL



%==========================================================================
%%               COMPUTE FISHER'S P-VALUE
%==========================================================================
clearvars -except IJ ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL
clc; disp('COMPUTING FISHERS EXACT TEST STATISTICS PER DNA LOCUS')


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


clearvars -except IJ ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL


%% SAVE MAT FILE
pause(1)

INFO.GOODCOHORTS = ADSP.GOODCOHORTS;
INFO.BRAKCOHORTS = ADSP.BRAKCOHORTS;
INFO.USE_COHORT  = ADSP.USE_COHORT;
INFO.USE_APOE    = ADSP.USE_APOE;
INFO.USE_APOT    = ADSP.USE_APOT;

dt=char(datetime(datetime,'Format','yyyy-MM-dd-HH-mm-ss'));
save(['APOE_' ADSP.USE_APOT '_FISHP_' dt '.mat'],...
    'LOCI','PHEN','TRCASE','TRCTRL','TECASE','TECTRL','INFO');
pause(1)
%------------------------------------------%

end















