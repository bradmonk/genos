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


GENO = load('GENOSDATA.mat');

% P.datapath1 = ['/Users/bradleymonk/Documents/MATLAB/GIT/genomics/genos/genos_data/APOE/'...
% 'APOE_22_23_24_33_34_44/APOE_22_23_24_33_34_44_NNWIN.mat'];
% P.alleles1 = 'APOE_22_23_24_33_34_44';



P.datapath1 = ['/Users/bradleymonk/Documents/MATLAB/GIT/genomics/genos/genos_data/APOE/'...
'APOE_22_23_24_33_34_44/APOE_22_23_24_33_34_44_FISHP_V2'];
P.alleles1 = 'APOE_22_23_24_33_34_44';


P.datapath2 = ['/Users/bradleymonk/Documents/MATLAB/GIT/genomics/genos/genos_data/APOE/'...
'APOE_22_23_24_33_34_44/APOE_22_23_24_33_34_44_FISHP_V2'];
P.alleles2 = 'APOE_22_23_24_33_34_44';


P.datapath3 = ['/Users/bradleymonk/Documents/MATLAB/GIT/genomics/genos/genos_data/APOE/'...
'APOE_22_23_24_33_34_44/APOE_22_23_24_33_34_44_FISHP_V3'];
P.alleles3 = 'APOE_22_23_24_33_34_44';



clc; disp('DATASETS LOADED');
who; clearvars -except P ADSP GENO



%==========================================================================
%%
clearvars -except P ADSP GENO



P.FILES.w1 = what(P.datapath1);
P.FILES.w2 = what(P.datapath2);
P.FILES.w3 = what(P.datapath3);


P.Nmats1 = numel(P.FILES.w1.mat);
P.Nmats2 = numel(P.FILES.w2.mat);
P.Nmats3 = numel(P.FILES.w3.mat);



%=====================================================

MAT.DAT1 = load([P.FILES.w1.path filesep P.FILES.w1.mat{1}]);
MAT.DAT2 = load([P.FILES.w2.path filesep P.FILES.w2.mat{1}]);
MAT.DAT3 = load([P.FILES.w3.path filesep P.FILES.w3.mat{1}]);


LOCI1 = MAT.DAT1.LOCI;
LOCI2 = MAT.DAT2.LOCI;
LOCI3 = MAT.DAT3.LOCI;

[~,i]  = sort(LOCI1.TRFISHP);
LOCI1  = LOCI1(i,:);

[~,i]  = sort(LOCI2.TRFISHP);
LOCI2  = LOCI2(i,:);

[~,i]  = sort(LOCI3.TRFISHP);
LOCI3  = LOCI3(i,:);


TAB1 = LOCI1(1:50,1:8);
TAB1.VID = (1:50)';

TAB2 = LOCI2(1:50,1:8);
TAB2.VID = (1:50)';

TAB3 = LOCI3(1:50,1:8);
TAB3.VID = (1:50)';


LOCITAB = [TAB1; TAB2; TAB3];

%=====================================================
% clearvars -except P ADSP GENO 

for IJ = 2:P.Nmatfiles

MAT.DAT1 = load([P.FILES.w1.path filesep P.FILES.w1.mat{IJ}]);
MAT.DAT2 = load([P.FILES.w2.path filesep P.FILES.w2.mat{IJ}]);
MAT.DAT3 = load([P.FILES.w3.path filesep P.FILES.w3.mat{IJ}]);


LOCI1 = MAT.DAT1.LOCI;
LOCI2 = MAT.DAT2.LOCI;
LOCI3 = MAT.DAT3.LOCI;

[~,i]  = sort(LOCI1.TRFISHP);
LOCI1  = LOCI1(i,:);

[~,i]  = sort(LOCI2.TRFISHP);
LOCI2  = LOCI2(i,:);

[~,i]  = sort(LOCI3.TRFISHP);
LOCI3  = LOCI3(i,:);


TAB1 = LOCI1(1:50,1:8);
TAB1.VID = (1:50)';

TAB2 = LOCI2(1:50,1:8);
TAB2.VID = (1:50)';

TAB3 = LOCI3(1:50,1:8);
TAB3.VID = (1:50)';


LOCITAB = [LOCITAB; TAB1; TAB2; TAB3];


end



%==========================================================================
%%
clearvars -except P ADSP GENO LOCITAB


[C,ia,ic] = unique(LOCITAB.CHRPOS,'stable');

LOCT = LOCITAB(ia,:);

LOCT.COUNT = accumarray(ic,1);

LOCT.SUMPOS = zeros(size(LOCT,1),1);


for i = 1:size(LOCT,1)

    LOCT.SUMPOS(i) = sum( LOCITAB.VID(ic==i) );

    %LOCT.SUMPOS(i) = sum( LOCITAB.VID(ic==i) );

end

LOCT.SUMPOSADJ = ((60-LOCT.COUNT).*55) + LOCT.SUMPOS;

LOCT.MEANPOS = LOCT.SUMPOSADJ ./ 60;


[~,i]  = sort(LOCT.MEANPOS,'ascend');
LOCT  = LOCT(i,:);


LOCT.MEANCOUNT = LOCT.COUNT ./ 60;




