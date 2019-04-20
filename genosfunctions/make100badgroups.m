function [TRCASE,TRCTRL,TECASE,TECTRL] = make100badgroups(PHEN)
%%          PREPARE TRAINING AND TESTING SETS
%% GENERATE TRAINING AND TESTING GROUPS (CONDITIONAL RANDOM SAMPLE)
% 1. Identify cohorts in PHEN with a reasonable number of CASE and CTRL
% 2. Take a random selection of people from those cohorts so that
%    there are the same number of CASE and CTRL.
% 3. The rest of the people from the keeper cohorts can be used to
%    increase the size of the NN test group.
clearvars -except ADSP LOCI CASE CTRL PHEN USNP
clc; close all;




% GET SET OF GOOD COHORTS
%-------------------------------------------------------------
% COHSET = PHEN(sum(PHEN.COHORTNUM == [1 2 6 7 9 10 11 12 13 19 20] , 2)>0,:);
% COHSET = COHSET(randperm(size(COHSET,1)),:);    % RAND SHUFFLE TABLE
% COHSET = PHEN(PHEN.GOODCOH==1,:);               % GET GOOD COHORTS
% COHSET = COHSET(randperm(size(COHSET,1)),:);    % RAND SHUFFLE TABLE

COHSET = PHEN;
COHSET = COHSET(randperm(size(COHSET,1)),:);    % RAND SHUFFLE TABLE





% disp('Usable cohorts:'); fprintf('%5i%5i%5i\n',unique(COHSET.COHORTNUM));
% fprintf('\nParticipants remaining: %5i \n',height(COHSET));




%% GET RANDOM PARTICIPANT SET FROM EACH COHORT

PHENCASE = COHSET( (COHSET.AD==1) ,:);
PHENCTRL = COHSET( (COHSET.AD==0) ,:);

NCa = size(PHENCASE,1);
NCo = size(PHENCTRL,1);
NTrCa = round(NCa * .7);
NTrCo = round(NCo * .7);

RTrCa = randperm( NCa , NTrCa);
RTrCo = randperm( NCo , NTrCo);

PHENCASE.TRAIN = zeros(NCa,1);
PHENCTRL.TRAIN = zeros(NCo,1);

PHENCASE.TRAIN(RTrCa) = 1;
PHENCTRL.TRAIN(RTrCo) = 1;

PHETRCASE = PHENCASE(PHENCASE.TRAIN==1,:);
PHETRCTRL = PHENCTRL(PHENCTRL.TRAIN==1,:);

PHETECASE = PHENCASE(PHENCASE.TRAIN==0,:);
PHETECTRL = PHENCTRL(PHENCTRL.TRAIN==0,:);




%% RANDOMIZE PHENOTYPE TABLES OF EACH SUBGROUP
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


clc; clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
PHETRCASE PHETRCTRL PHETECASE PHETECTRL




%% PERFORM SOME CHECKS ON THE FINAL TRAINING AND TESTING GROUPS
%{.
clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
PHETRCASE PHETRCTRL PHETECASE PHETECTRL




% COMPUTE Z-SCORE FOR CASE-AGE & CTRL-AGE, SEPARATELY
%------------------------------------------------
PHECASE = [PHETRCASE; PHETECASE];
PHECASE.AGE(PHECASE.AGE<55) = 55;
PHECASE.AGE = round(PHECASE.AGE);
PHECASE.AGEZ = zeros(size(PHECASE,1),1);
PHECASE.AGEZ = zscore(PHECASE.AGE);       % KEEP SIGN OF CASE GROUP Z

PHECTRL = [PHETRCTRL; PHETECTRL];
PHECTRL.AGE(PHECTRL.AGE<60) = 60;
PHECTRL.AGE = round(PHECTRL.AGE);
PHECTRL.AGEZ = zeros(size(PHECTRL,1),1);
PHECTRL.AGEZ = zscore(PHECTRL.AGE);       % REMEMBER TO FLIP SIGN OF CTRL GROUP Z
% PHECTRL.AGEZ = zscore(PHECTRL.AGE) .* -1; % REMEMBER TO FLIP SIGN OF CTRL GROUP Z
%------------------------------------------------



% CREATE A BRAAK-BASED AGE: BRAGE = AGE - BRAAK
%------------------------------------------------
CASENAN = isnan(PHECASE.BRAAK);
CTRLNAN = isnan(PHECTRL.BRAAK);

PHECASE.BRAK = CASENAN .* -1;      % MAKE NAN EQUAL -1 FOR NOW, ZERO LATER
PHECTRL.BRAK = CTRLNAN .* -1;      % MAKE NAN EQUAL -1 FOR NOW, ZERO LATER

PHECASE.BRAK(~CASENAN) = PHECASE.BRAAK(~CASENAN);
PHECTRL.BRAK(~CTRLNAN) = PHECTRL.BRAAK(~CTRLNAN);

PHECASE.BRAK(PHECASE.BRAK==0) = -6.4;
PHECTRL.BRAK(PHECTRL.BRAK==0) = -6.4;

PHECASE.BRAK(PHECASE.BRAK==1) = -4.4;
PHECTRL.BRAK(PHECTRL.BRAK==1) = -4.4;

PHECASE.BRAK(PHECASE.BRAK==2) = -2.4;
PHECTRL.BRAK(PHECTRL.BRAK==2) = -2.4;

PHECASE.BRAK(PHECASE.BRAK==3) = 1.5;
PHECTRL.BRAK(PHECTRL.BRAK==3) = 1.5;

PHECASE.BRAK(PHECASE.BRAK==4) = 4.5;
PHECTRL.BRAK(PHECTRL.BRAK==4) = 4.5;

PHECASE.BRAK(PHECASE.BRAK==5) = 6.5;
PHECTRL.BRAK(PHECTRL.BRAK==5) = 6.5;

PHECASE.BRAK(PHECASE.BRAK>=6) = 8.5;
PHECTRL.BRAK(PHECTRL.BRAK>=6) = 8.5;


PHECASE.BRAK = PHECASE.BRAK+1;      % ADD +1 BECAUSE NAN WAS MADE -1 ABOVE
PHECTRL.BRAK = PHECTRL.BRAK+1;      % ADD +1 BECAUSE NAN WAS MADE -1 ABOVE

PHECASE.BRAGE = PHECASE.AGE - PHECASE.BRAK;
PHECTRL.BRAGE = PHECTRL.AGE - PHECTRL.BRAK;
%------------------------------------------------


% COMPUTE Z-SCORE FOR CASE-AGE & CTRL-AGE, SEPARATELY
%------------------------------------------------
PHECASE.BRAGEZ = zeros(size(PHECASE,1),1);
PHECASE.BRAGEZ = zscore(PHECASE.BRAGE);       % FLIP SIGN OF CASE Z BELOW

PHECTRL.BRAGEZ = zeros(size(PHECTRL,1),1);
PHECTRL.BRAGEZ = zscore(PHECTRL.BRAGE);       % KEEP SIGN OF CTRL Z BELOW
%------------------------------------------------


clearvars -except ADSP LOCI CASE CTRL PHEN USNP PHECASE PHECTRL


%###############################################################
%%   FLIP SIGN OF CONTROL BRAGEZ & AGEZ VARIABLES
%###############################################################


PHECASE.AGEZ   = zscore(PHECASE.AGE) .* -1; 

PHECASE.BRAGEZ = zscore(PHECASE.BRAGE) .* -1;




%###############################################################
%%   SAVE BACK INTO TRAINING AND TESTING SETS
%###############################################################


TRCASE = PHECASE(PHECASE.TRAIN==1 ,:);
TRCTRL = PHECTRL(PHECTRL.TRAIN==1 ,:);
TECASE = PHECASE(PHECASE.TRAIN==0 ,:);
TECTRL = PHECTRL(PHECTRL.TRAIN==0 ,:);

% PHETRCASE = PHENCASE(PHENCASE.TRAIN==1,:);
% PHETRCTRL = PHENCTRL(PHENCTRL.TRAIN==1,:);
% PHETECASE = PHENCASE(PHENCASE.TRAIN==0,:);
% PHETECTRL = PHENCTRL(PHENCTRL.TRAIN==0,:);


clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
TRCASE TRCTRL TECASE TECTRL


%###############################################################
%%  USE PER-PERSON VARIANT COUNTS TO IDENTIFY OUTLIERS
%###############################################################

% Qlo = quantile(PHETRCASE.TOTvars,.01);
% Qhi = quantile(PHETRCASE.TOTvars,.99);
% TRCASE  = PHETRCASE(((PHETRCASE.TOTvars > Qlo) & (PHETRCASE.TOTvars < Qhi)),:);
% 
% 
% Qlo = quantile(PHETRCTRL.TOTvars,.01);
% Qhi = quantile(PHETRCTRL.TOTvars,.99);
% TRCTRL  = PHETRCTRL(((PHETRCTRL.TOTvars > Qlo) & (PHETRCTRL.TOTvars < Qhi)),:);
% 
% 
% Qlo = quantile(PHETECASE.TOTvars,.01);
% Qhi = quantile(PHETECASE.TOTvars,.99);
% TECASE  = PHETECASE(((PHETECASE.TOTvars > Qlo) & (PHETECASE.TOTvars < Qhi)),:);
% 
% 
% Qlo = quantile(PHETECTRL.TOTvars,.01);
% Qhi = quantile(PHETECTRL.TOTvars,.99);
% TECTRL  = PHETECTRL(((PHETECTRL.TOTvars > Qlo) & (PHETECTRL.TOTvars < Qhi)),:);


% clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
% TRCASE TRCTRL TECASE TECTRL

%%
end