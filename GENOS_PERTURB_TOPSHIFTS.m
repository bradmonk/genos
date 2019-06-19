%% GENOS: 
%{
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
% APOE4 CHRPOS: 190045411941
% APOE2 CHRPOS: 190045412079
%
% APOE  190045411941  190045412079
% e22:      -1            5
% e23:      -1            2
% e24:       2            2
% e33:      -1           -1
% e34:       2           -1
% e44:       5           -1
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
%}
%--------------------------------------------------------------------------

%==========================================================================
%% STEP-1: LOAD THE DATASET
%==========================================================================
close all; clear; clc; rng('shuffle');
P.home = fileparts(which('GENOS.m')); cd(P.home);
P.P1 = [P.home filesep 'genos_functions'];
P.P3 = [P.P1 filesep 'genos_main_functions'];
P.P4 = [P.home filesep 'genos_other'];
addpath(join(string(struct2cell(P)),pathsep,1))
cd(P.home)


ADSP = load('GENOSDATAFINAL.mat');


P.mainmatfile = which('GENOSDATAFINAL.mat');
disp('dataset loaded')
clearvars -except P ADSP INFO



%==========================================================================
%% APOE related options
%==========================================================================
clc; clearvars -except P ADSP INFO

P.basedir = '/Users/bradleymonk/Documents/MATLAB/GIT/genomics/genos/genos_data/APOE';
f = filesep;


P.APOES = '22_23_24_33_34_44'; INFO.APOE = [22 23 24 33 34 44];
P.importdir = [P.basedir f 'APOE_22_23_24_33_34_44' f 'APOE_22_23_24_33_34_44_FISHP'];

% P.APOES = '22_23_24_34_44'; INFO.APOE = [22 23 24 34 44];
% P.importdir = [P.basedir f 'APOE_22_23_24_34_44' f 'APOE_22_23_24_34_44_FISHP'];

% P.APOES = '33'; INFO.APOE = [33];
% P.importdir = [P.basedir f 'APOE_33' f 'APOE_33_FISHP'];

% P.APOES = '34_44'; INFO.APOE = [34 44];
% P.importdir = [P.basedir f 'APOE_34_44' f 'APOE_34_44_FISHP'];

% P.APOES = '33_44'; INFO.APOE = [34 44];
% P.importdir = [P.basedir f 'APOE_34_44' f 'APOE_34_44_FISHP'];


clearvars -except P ADSP INFO





%==========================================================================
%% (GENOX) PERTURB TEST OPTIONS & FOLDERS 
%==========================================================================
clc; clearvars -except P ADSP INFO


% SET PERTURB TEST OPTIONS
%------------------------------------------------------
P.f = filesep;
P.NGenes = 100;
P.Nloops = 10;
P.FileStart = 1;
P.Nvars = 200;
P.windowSize = 50;
P.Ndots = 5;
P.Lo2Hi = true;
P.RemoveGenesByName = false;
P.RemoveBadGenes = false;



% SET DEFAULT FOLDERS 
%------------------------------------------------------
P.genox.DIRroot = [P.home P.f 'genos_data' P.f 'GENOX'];


% P.genox.DIRmat  = [P.genox.DIRroot P.f 'PLO_MAT'];
% P.genox.DIRimg  = [P.genox.DIRroot P.f 'PLO_IMG'];
% P.genox.DIRstatsmat = [P.genox.DIRroot P.f 'PLO_STATS_MAT'];
% P.genox.DIRstatsimg = [P.genox.DIRroot P.f 'PLO_STATS_IMG'];

% P.genox.DIRimg  = [P.genox.DIRroot P.f 'ORHI_MAT'];
% P.genox.DIRimg  = [P.genox.DIRroot P.f 'ORHI_IMG'];
% P.genox.DIRstatsmat = [P.genox.DIRroot P.f 'ORHI_STATS_MAT'];
% P.genox.DIRstatsimg = [P.genox.DIRroot P.f 'ORHI_STATS_IMG'];

P.genox.DIRimg  = [P.genox.DIRroot P.f 'ORLO_MAT'];
P.genox.DIRimg  = [P.genox.DIRroot P.f 'ORLO_IMG'];
P.genox.DIRstatsmat = [P.genox.DIRroot P.f 'ORLO_STATS_MAT'];
P.genox.DIRstatsimg = [P.genox.DIRroot P.f 'ORLO_STATS_IMG'];




% GET PATHS TO ALL MAT FILES IN SELECTED FOLDER
%------------------------------------------------------
P.PERT.w = what(P.genox.DIRstatsmat);
P.PERT.finfo = dir(P.PERT.w.path);
P.PERT.finames = {P.PERT.finfo.name};
c=~cellfun(@isempty,regexp(P.PERT.finames,'((\S)+(\.mat+))'));
P.PERT.finames = string(P.PERT.finames(c)');
P.PERT.folder = P.PERT.finfo.folder;
P.PERT.fipaths = fullfile(P.PERT.folder,P.PERT.finames);
disp(P.PERT.fipaths); disp(P.PERT.finames);


GENOS = load(P.PERT.fipaths{1});
%---
STATS(1).GENE   = GENOS.PER.GENE;
STATS(1).CHRPOS = GENOS.PER.CHRPOS;
STATS(1).AREA   = GENOS.PER.area;
STATS(1).CASEMU = GENOS.PER.caseMu;
STATS(1).CTRLMU = GENOS.PER.ctrlMu;
STATS(1).CACOMU = GENOS.PER.cacoZd;





GENOS = load(P.PERT.fipaths{1});

TAB = struct2table(GENOS.PER);
TAB.GENE(2) = TAB.GENE(1);
TAB.GENE(3) = TAB.GENE(1);
TAB.GENE(4) = TAB.GENE(1);

TAB.CHRPOS(2) = TAB.CHRPOS(1);
TAB.CHRPOS(3) = TAB.CHRPOS(1);
TAB.CHRPOS(4) = TAB.CHRPOS(1);

TAB.TEST = [1;2;3;4];

STATS = TAB;



for ii = 2:numel(P.PERT.fipaths)

    GENOS = load(P.PERT.fipaths{ii});

    TAB = struct2table(GENOS.PER);
    TAB.GENE(2) = TAB.GENE(1);
    TAB.GENE(3) = TAB.GENE(1);
    TAB.GENE(4) = TAB.GENE(1);

    TAB.CHRPOS(2) = TAB.CHRPOS(1);
    TAB.CHRPOS(3) = TAB.CHRPOS(1);
    TAB.CHRPOS(4) = TAB.CHRPOS(1);

    TAB.TEST = [1;2;3;4];

    STATS = [STATS; TAB];

end


clearvars -except P ADSP INFO PERT STATS


%==========================================================================
%% PLOT AREA HISTOGRAM
%==========================================================================

STATISTICS = STATS(:,[1 2 4 5 6 7]);






