
%==========================================================================
%% STEP-1: LOAD THE DATASET
%==========================================================================
close all; clear; clc; rng('shuffle');
P.home = fileparts(which('GENOS.m')); cd(P.home);
P.P1 = [P.home filesep 'genos_functions'];
P.P3 = [P.P1 filesep 'genos_main_functions'];
P.P2 = [P.home filesep 'genosfunctions'];
P.P4 = [P.home filesep 'genos_other'];
addpath(join(string(struct2cell(P)),pathsep,1))
cd(P.home)


which('GENOSDATA.mat')
ADSP = load('GENOSDATA.mat');

which('ADSP_MINI3.mat')
ADSP = load('ADSP_MINI3.mat');

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


clc; clearvars -except P ADSP INFO PHEN LOCI CASE CTRL USNP
head(PHEN)
head(LOCI)


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

%==========================================================================
%%             GET APOE AND TOMM40 LOCI
%==========================================================================
clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL


ALZGENES = string({"APOE";"TOMM40"});

x1 = strcmp(VLOCI.GENE,ALZGENES{1});
sum(x1)
x2 = strcmp(VLOCI.GENE,ALZGENES{2});
sum(x2)

KEEP = x1 | x2;
VLOCI = VLOCI(KEEP,:);
VCASE = VCASE(KEEP);
VCTRL = VCTRL(KEEP);
VUSNP = VUSNP(KEEP);

%%


% ApoE status is defined by a person's combination of these two SNPs:
% 
%     SNP        CHR   POS(GRCh37)   POS(GRCh38)  REF   ALT   MUT
%     rs429358   19    45411941      44908684     T     C     missense
%     rs7412     19    45412079      44908822     C     T     missense
% 
%     APOE 2/2	45411941 ( T , T )	45412079 ( T , T )
%     APOE 2/3	45411941 ( T , T )	45412079 ( T , C )
%     APOE 2/4	45411941 ( T , C )	45412079 ( T , C )
%     APOE 3/3	45411941 ( T , T )	45412079 ( C , C )
%     APOE 3/4	45411941 ( T , C )	45412079 ( C , C )
%     APOE 4/4	45411941 ( C , C )	45412079 ( C , C )
   

% you never two mutations on the same chromosome:
%     APOE 1	45411941 ( C )	45412079 ( T )
   

%==========================================================================
%%      MAKE  RECTANGLE  NEURAL NET  VARIANT MATRIX
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
% [VTRX, TRX, TRL] = mkmx(VLOCI,VCASE,VCTRL,VUSNP,TRPHE,[-1 -0 2 5]);
% [VTEX, TEX, TEL] = mkmx(VLOCI,VCASE,VCTRL,VUSNP,TEPHE,[-1 -0 2 5]);
[VTRX, TRX, TRL] = mkmx(VCASE,VCTRL,VUSNP,TRPHE,[-1 -0 2 3]);
[VTEX, TEX, TEL] = mkmx(VCASE,VCTRL,VUSNP,TEPHE,[-1 -0 2 3]);



clearvars -except P ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL TRPHE TEPHE...
VTRX VTEX TRX TRL TEX TEL










