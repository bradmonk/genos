%% MACHINE LEARNING CLASSIFIER FOR GENOMICS ALZHEIMERS DIAGNOSIS
% The following script will import a genomics dataset .mat file
% that contains four variables... (expando below see extensive notes)
%{

1. DATATABLE
2. SNPCASE
3. SNPCTRL
4. PHE


Below, I will clarify what each of those variables contains, but
First, let me define the term 'VARIANT', since I'll be using it a lot to
describe the organization of the genomics data. 

For the purposes of genomics studies, and as an approximate truth 
particularly within ethnic populations...

    **humans are genetically identical**

...Given this approximate truth, we as a research community have 
synthesized a prototypical human genome. This prototype reference genome
has undergone several revisions. The most current rev is: GRCh38.p11
However, the genomes here were compared to a prior reference assembly. It
happens to be the second latest assembly: GRCh37



-------------
NB
For more info on these revisions see:
https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/

Also if you want to know what genetic ethnic ancestry looks like,
visit this link and view the video halfway down the page...

https://macarthurlab.org/2017/02/27/the-genome-aggregation-database-gnomad

There are other good reasons to visit that link, if you are interested in
random forest filtering, or standard genomics data QC pipelines, or other
things that can fall out of PCA-related methods.
-------------


With this reference genome in-hand, one can then compare any individual
person's genome, nucleotide-by-nucleotide to see if any of these
nucleotide base-pairs differ...

REFERENCE:  ATCGAAGTCC
INDIVIDUAL: ATCGAAGACC
                   |
                variant


Whenever we find a variant, we create a row in a table of variants,
and note it's genomic location (CHR, POS). If that row already exists
in the table, we add it to the tally. If more than one group is being
assessed (as in most CASE/CTRL studies), we add that variant to the
tally for the appropriate group. We've created such a table for
~5000 Alzheimers CASE and ~5000 age-matched CTRL people.

After importting the .MAT file, the first thing that happens is a few 
ancillary columns are removed from DATATABLE. After this though, 
if we display a few rows from DATATABLE, this is what it should look like:

>> disp(DATATABLE(1:8 , :))

    VID    CHR     POS      REF    ALT    EFFECT       GENE       CASEREF    CTRLREF    CASEALT    CTRLALT     FISHP      FISHOR
    ___    ___    ______    ___    ___    ______    __________    _______    _______    _______    _______    ________    ______
    1      1      861315    G      A      UTR       SAMD11        11858      9876       0          1           0.45443       Inf
    2      1      861333    G      T      SYN       SAMD11        11870      9894       0          1           0.45463       Inf
    3      1      861343    G      C      MIS       SAMD11        11866      9904       1          0                 1         0
    4      1      861349    C      T      MIS       SAMD11        11870      9894       0          3          0.093982       Inf
    5      1      861356    T      C      MIS       SAMD11        11874      9894       0          1           0.45455       Inf
    6      1      861363    C      T      SYN       SAMD11        11866      9892       1          0                 1         0
    7      1      861368    C      T      MIS       SAMD11        11862      9878       1          1                 1    1.2009
    8      1      861384    C      G      SYN       SAMD11        11836      9882       1          0                 1         0


VID is just counting from 1:Nrows
CHR is the chromosome of the given variant
POS is the position on that chromosome
REF is the nucleotide we expected based on the reference prototype
ALT is the nucleotide that was actually found there


EFFECT    is the consequence of that particular switch. 
          for this, just know that 'SYN' isn't quite as consequential as the others EFFECTS


GENE     is the name of the gene where the variant was found


CASEREF  is the total number of CASE people x2 (because each person has 2 chromosomes)
         that have the same nucleotide as the prototype


CTRLREF  is the total number of CTRL people x2 (because each person has 2 chromosomes)
         that have the same nucleotide as the prototype


CASEALT  is the total number of CASE variants with a nucleotide at that position 
         that differs from the prototype

CTRLALT  is the total number of CTRL variants with a nucleotide at that position 
         that differs from the prototype





1. DATATABLE (check)
2. SNPCASE
3. SNPCTRL
4. PHE


SNPCASE and SNPCTRL simply contain ID numbers for each person who had
a variant at each row in the DATATABLE. So you will notice that all
three of these variables has exactly the same number of rows.


1. DATATABLE (check)
2. SNPCASE (check)
3. SNPCTRL (check)
4. PHE


PHE contains the phenotype information for each participant. The ID numbers
found in SNPCASE and SNPCTRL, every single possible one of those IDs should
be listed in the first column of the PHE table. If we display a few rows of
this table...

>> disp(PHE(1:6 , :))

       SRR        SEX     AGE      APOE    AUTOPSY    BRAAK    PREVAD    INCAD    AD
    __________    ___    ______    ____    _______    _____    ______    _____    __
    1.3487e+06    1          88    33      0          NaN      0         0        0 
    1.2681e+06    1       89.85    33      0          NaN      0         1        1 
    1.4504e+06    1      87.488    23      0          NaN      0         0        0 
    1.1253e+06    1          90    23      0          NaN      0         0        0 
     1.304e+06    0          73    33      1            4      1         0        1 
    1.2797e+06    1          87    33      0          NaN      0         0        0 


...we can see there are other phenotypic features that may be useful. The
major one to keep an eye on is "AD". This indicates whether the person was
diagnosed with Alzheimers (is a CASE), or was an age-matched non-diagnosed
person (is a CTRL). The other columns are either self-explainatory, or
too complicated to explain given their very low utility in any sort of 
statistical analysis; other than 'APOE', which is a known risk factor 
for Alzheimers. People with higher APOE scores are more succeptible.
Particularly if they have 34 or 44 APOE status (all possible APOE combos
are: 22 23 33 34 44 24).


1. DATATABLE (check)
2. SNPCASE (check)
3. SNPCTRL (check)
4. PHE  (check)


READY!......

%}


%% IMPORT TRAINING DATA
clc; close all; clear
% system('sudo purge')
cd(fileparts(which('EDNN.m')));



% LOAD MAIN .MAT FILE
load('GENOMICSDATA_MAINN.mat')


% REMOVE ANCILLARY COLUMNS FROM DATATABLE
DATATABLE.GCASEREF  = [];
DATATABLE.GCTRLREF  = [];
DATATABLE.GCASEALT  = [];
DATATABLE.GCTRLALT  = [];
DATATABLE.GENEi     = [];
DATATABLE.A2A       = [];
DATATABLE.AAN       = [];
DATATABLE.PPHEN     = [];


disp('DATATABLE:')
disp(DATATABLE(1:8 , :))





% REMOVE INSTANCES WHERE ALT ALLELE REALLY SHOULD BE THE REFERENCE
ALTDOM = (DATATABLE.CASEALT > DATATABLE.CASEREF) | ...
         (DATATABLE.CTRLALT > DATATABLE.CTRLREF);

disp(sum(ALTDOM));                 % how many ALT alleles removed?
disp(numel(ALTDOM));               % number of total variants?
disp(sum(ALTDOM)/numel(ALTDOM));   % what is the percent removed

DATATABLE(ALTDOM,:)  = [];
SNPCASE(ALTDOM)      = [];
SNPCTRL(ALTDOM)      = [];


% WHENEVER *ROWS* ARE REMOVED FROM DATATABLE REASSIGN VID NUMBERS
DATATABLE.VID = ( 1:size(DATATABLE,1) )';




% MAKE GENE NAMES READABLE IN COMMAND WINDOW DISPLAY
DATATABLE.GENE = DATATABLE.GENE(:,1:10);




clc; clearvars -except DATATABLE SNPCASE SNPCTRL PHE
disp('DATATABLE:')
disp(DATATABLE(1:8 , :))





%% GET SOME SUBSET OF THE FULL DATA TABLE FOR MACHINE LEARNING

% Data must be converted from sparse format to full matrix format.
%   that means each person gets their own row, and each variant
%   gets its own column. Because of this, the representation in
%   RAM will baloon significantly. So we can only test about
%   500 variants on a typical personal computer setup.
%   With this in mind, the variant subset selected for training
%   need not be random. However, the ideal set is an empirical
%   question. Here I'm selecting variants with high assymetry
%   (low Fisher's P-values), and a few other criteria...


ASYMX = DATATABLE( (((DATATABLE.CASEALT + DATATABLE.CTRLALT)>30) & ...
                   ((DATATABLE.CASEALT + DATATABLE.CTRLALT)<3000))  & ...
                   (DATATABLE.FISHP < .2) & ...
                   (DATATABLE.FISHOR>1.2 | DATATABLE.FISHOR<0.8) , :);


[J , I] = sortrows(ASYMX.FISHP);

ASYMX = ASYMX(I,:);


%-----------------------------------------------
% CREATE TABLE GUI TO DISPLAY TOP VARIANT SITES
clc; close all
fh1=figure('Units','normalized','OuterPosition',[.01 .05 .95 .6],'Color','w','MenuBar','none');
t = uitable(fh1,'Data',[ASYMX.Properties.VariableNames; table2cell(ASYMX(1:25,:))],...
'Units','normalized','Position',[.01 .01 .98 .95]);
%-----------------------------------------------





%% MAKE NEURAL NET TRAINING MATRIX

CASEID = PHE.SRR(PHE.AD==1);
CTRLID = PHE.SRR(PHE.AD==0);

VID   = ASYMX.VID;
caSNP = SNPCASE(VID);
coSNP = SNPCTRL(VID);
caMX  = zeros( size(CASEID,1) , size(VID,1) );
coMX  = zeros( size(CTRLID,1) , size(VID,1) );


for nn = 1:size(caSNP,1)
    
    caMX(:,nn) = ismember(CASEID , caSNP{nn} );
    coMX(:,nn) = ismember(CTRLID , coSNP{nn} );
    
    if ~mod(nn,1000); disp(nn/size(caSNP,1)); end
end

cacoMX = [caMX;coMX];

ADNN = padarray(cacoMX,[1 3],0,'pre');

nc=length(CASEID)+1;

ADNN(2:end , 3)  =  1;                  % COL#3: bias
ADNN(2:nc ,  2)  =  1;                  % COL#2: case/ctrl 1/0
ADNN(2:end , 1)  =  [CASEID; CTRLID];   % COL#1: subject IDs
ADNN(1 , 4:end)  =  ASYMX.VID';         % ROW#1: VID


disp('COVERAGE...')
disp(   sum(sum( ADNN(2:end , 4:end)  ,  2) > 0) / (size(CTRLID,1)+size(CASEID,1)) )


% CLEAN UP WORKSPACE
clc; ADNN(1:8,1:8)
clearvars -except DATATABLE SNPCASE SNPCTRL PHE CASEID CTRLID ASYMX ADNN





%% PREP DATA FOR NN CLASSIFIER TRAINING


% SETABLE PARAMETERS
%-----------------------------------
UseNgenes = 70;   % 140 genes work best
maxiters  = 21;    % 21 iterations is perfect
lambda = 0.006;    % .004 - .008 is best window
L2neurons = 23;    % 17 - 30 neurons works best
useBias = 1;       % always use bias
%-----------------------------------


INPUTMX       = ADNN(2:end,3:end);  % col 3 of DATATABLE is bias
INPUTYN       = ADNN(2:end,2);

ADNN(1:9,1:9)
INPUTMX(1:9,1:9)


% DATATABLE
%     1       2       3      4      5      6 ...
%------------------------------------------------
%1|   0       0       0  SNPid  SNPid  SNPid
%2| SRR    CaCo    bias    1/0    1/0    1/0
%3| SRR    CaCo    bias    1/0    1/0    1/0
%4| SRR    CaCo    bias    1/0    1/0    1/0


% INPUTMX
%      1      2      3      4
%------------------------------------------------
%1| bias    1/0    1/0    1/0
%2| bias    1/0    1/0    1/0
%3| bias    1/0    1/0    1/0


INPUTMXdims   = size(INPUTMX);
INPUTYNdims   = size(INPUTYN);

% RANDOMIZE ROWS
randp         = randperm(INPUTMXdims(1));
INPUTMX       = INPUTMX(randp,:);
INPUTYN       = INPUTYN(randp,:);


disp('INPUTMX:')
disp(INPUTMX(1:9,1:8));    % now col 1 of INPUTMX is bias


clearvars -except ...
NNET ADGENES INPUTMX INPUTYN INPUTMXbias ...
UseNgenes UseGenes maxiters lambda L2neurons useBias xlsN xlsG ...
initTheta1 initTheta2 DATATABLE SNPCASE SNPCTRL PHE CASEID CTRLID ASYMX ADNN





%% LOGISTIC REGRESSION NEURAL NET WITH GRADIENT DECENT

sig    = @(z) 1 ./ (1 + exp(-z) );     % sigmoid activation function
sigrad = @(g) sig(g) .* (1-sig(g));     % sigmoid gradient function


UseGenes = [1:65];

X  = INPUTMX( : , UseGenes );         % INCLUDE BIAS COL

y = INPUTYN + 1;    % ADD 1 SO LABELS ARE {1,2}


m = size(X,1);      % ROWS: NUMBER OF TRAINING EXAMPLES (PEOPLE)
n = size(X,2);      % COLS: NUMBER OF FEATURES (SNPs)


L1neurons  = n;
num_labels = 2;

% RANDOMIZE INITIAL THETA WEIGHTS
eps_init = 0.12;

initTheta1 = rand(L2neurons, L1neurons+1) * 2 * eps_init - eps_init;
initTheta2 = rand(num_labels, L2neurons+1) * 2 * eps_init - eps_init;


% UNROLL THETA PARAMETERS 
initial_Thetas = [initTheta1(:) ; initTheta2(:)];
nn_Thetas = initial_Thetas;


% EXECUTE COST FUNCTION

J = GXNNCostF(nn_Thetas, L1neurons, L2neurons, num_labels, X, y, lambda);





%% TRAIN NEURAL NETWORK
disp('Training Neural Network...')


options = optimset('MaxIter', maxiters);

GcostFun = @(T) GXNNCostF(T, L1neurons, L2neurons, num_labels, X, y, lambda);

[nn_Thetas, cost] = GXfmincg(GcostFun, initial_Thetas, options);


% REROLL Theta1 AND Theta2
Theta1 = reshape(nn_Thetas(1:L2neurons * (L1neurons + 1)), L2neurons, (L1neurons + 1));
Theta2 = reshape(nn_Thetas((1 + (L2neurons * (L1neurons + 1))):end), num_labels, (L2neurons + 1));


%% HOW DOES NN PERFORM ON **TRAINING** SET?

[p , a , h] = GXNNpredict(Theta1, Theta2, X);

TRAINPCTCORRECT = mean(p == y);

disp('Percent accuracy on training data:')
disp( TRAINPCTCORRECT )


%% EVALUATE NEURAL NETWORK ON **TEST** SET


INPUTMX       = ADNN(2:end,3:end);  % col 3 of DATATABLE is bias
INPUTYN       = ADNN(2:end,2);

ADNN(1:9,1:9)
INPUTMX(1:9,1:9)

INPUTMXdims   = size(INPUTMX);
INPUTYNdims   = size(INPUTYN);


% RANDOMIZE ROWS
randp         = randperm(INPUTMXdims(1));
INPUTMX       = INPUTMX(randp,:);
INPUTYN       = INPUTYN(randp,:);

X  = INPUTMX( : , UseGenes );         % INCLUDE BIAS COL

y = INPUTYN + 1;    % ADD 1 SO LABELS ARE {1,2}

[p , a , h] = GXNNpredict(Theta1, Theta2, X);

TESTPCTCORRECT = mean(p == y);

% REPORT ACCURACY FOR FULL TESTING DATASET
disp('Percent accuracy on full testing data:')
disp( TESTPCTCORRECT )




%% EXAMINE & REPORT ACCURACY FOR HIGH CONFIDENCE DATA
LowConfThresh = .51;
MidConfThresh = .80;
HiConfThresh  = .85;


numCorrect = p == y;


LowConf      = a >= LowConfThresh & a < MidConfThresh;
LowConfPCT   = (mean(p(LowConf) == y(LowConf)))*100;
LowConfNp    = numel(p(LowConf));
LowConfPCTp  = (numel(p(LowConf)) / numel(p))*100;

MidConf      = a >= MidConfThresh & a < HiConfThresh;
MidConfPCT   = (mean(p(MidConf) == y(MidConf)))*100;
MidConfNp    = numel(p(MidConf));
MidConfPCTp  = (numel(p(MidConf)) / numel(p))*100;

HiConf       = a >= HiConfThresh;                   % Logical index of high confidence predictions
HiConfPCT    = (mean(p(HiConf) == y(HiConf)))*100;  % Percent correct hi conf predictions
HiConfNp     = numel(p(HiConf));                    % Total number of hi conf predictions
HiConfPCTp   = (numel(p(HiConf)) / numel(p))*100;   % Percent of hi conf predictions
HiConfNCorr  = HiConfNp * (HiConfPCT / 100);        % Total number of correct hi conf predictions


fprintf('\nPercent accuracy on mild conf (a > %.2f) testing data: >>> %.1f%% <<<\n',...
    LowConfThresh, LowConfPCT)
fprintf('    (includes %.0f of %.0f [%.1f%%] HighConf:Total patients) \n\n',...
    LowConfNp , numel(p) , LowConfPCTp )

fprintf('\nPercent accuracy on mid conf (a > %.2f) testing data: >>> %.1f%% <<<\n',...
    MidConfThresh, MidConfPCT)
fprintf('    (includes %.0f of %.0f [%.1f%%] HighConf:Total patients) \n\n',...
    MidConfNp , numel(p) , MidConfPCTp )

fprintf('\nPercent accuracy on high conf (a > %.2f) testing data: >>> %.1f%% <<<\n',...
    HiConfThresh, HiConfPCT)
fprintf('    (includes %.0f of %.0f [%.1f%%] HighConf:Total patients) \n\n',...
    HiConfNp , numel(p) , HiConfPCTp )



%% WHICH GENES WERE USED? WHICH WERE IMPORTANT?


% nprtool
% load('net1.mat')
% yguess = net(X');
% mean(INPUTYN' == round(yguess))
% p = round(yguess)' + 1;
