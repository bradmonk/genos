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






% APOE_22_23_24_33_34_44
%==================================
%{
P.datapath1 = ['/Users/bradleymonk/Documents/MATLAB/GIT/genomics/genos/genos_data/APOE'...
'/APOE_22_23_24_33_34_44/APOE_22_23_24_33_34_44_FISHP_V2'];
P.alleles1 = 'APOE_22_23_24_33_34_44';


P.datapath2 = ['/Users/bradleymonk/Documents/MATLAB/GIT/genomics/genos/genos_data/APOE'...
'/APOE_22_23_24_33_34_44/APOE_22_23_24_33_34_44_FISHP_V2'];
P.alleles2 = 'APOE_22_23_24_33_34_44';


P.datapath3 = ['/Users/bradleymonk/Documents/MATLAB/GIT/genomics/genos/genos_data/APOE'...
'/APOE_22_23_24_33_34_44/APOE_22_23_24_33_34_44_FISHP_V3'];
P.alleles3 = 'APOE_22_23_24_33_34_44';
%}







% APOE_22_23_24_34_44
%==================================
%{

P.datapath{1} = ['/Users/bradleymonk/Documents/MATLAB/GIT/genomics/genos/genos_data/APOE'...
'/APOE_22_23_24_34_44/APOE_22_23_24_34_44_FISHP_V2'];
P.alleles{1} = 'APOE_22_23_24_34_44';

P.datapath{2} = ['/Users/bradleymonk/Documents/MATLAB/GIT/genomics/genos/genos_data/APOE'...
'/APOE_22_23_24_34_44/APOE_22_23_24_34_44_FISHP_V3'];
P.alleles{2} = 'APOE_22_23_24_34_44';

%}



% APOE_33
%==================================
%{.

P.datapath{1} = ['/Users/bradleymonk/Documents/MATLAB/GIT/genomics/genos/genos_data/APOE'...
'/APOE_33/APOE_33_FISHP_P01-50_V2'];
P.alleles{1} = 'APOE_33';

P.datapath{2} = ['/Users/bradleymonk/Documents/MATLAB/GIT/genomics/genos/genos_data/APOE'...
'/APOE_33/APOE_33_FISHP_V4'];
P.alleles{2} = 'APOE_33';

%}






clc; disp('DATASETS LOADED');
who; clearvars -except P ADSP GENO



%==========================================================================
%%
clearvars -except P ADSP GENO




for i = 1:numel(P.datapath)

P.FILES.w{i} = what(P.datapath{i});

P.Nmats{i} = numel(P.FILES.w{i}.mat);

end


%=====================================================

MAT.DAT1 = load([P.FILES.w{1}.path filesep P.FILES.w{1}.mat{1}]);
MAT.DAT2 = load([P.FILES.w{2}.path filesep P.FILES.w{2}.mat{1}]);


LOCI1 = MAT.DAT1.LOCI;
LOCI2 = MAT.DAT2.LOCI;

[~,i]  = sort(LOCI1.TRFISHP);
LOCI1  = LOCI1(i,:);

[~,i]  = sort(LOCI2.TRFISHP);
LOCI2  = LOCI2(i,:);


TAB1 = LOCI1(1:50,1:8);
TAB1.VID = (1:50)';

TAB2 = LOCI2(1:50,1:8);
TAB2.VID = (1:50)';


LOCITAB = [TAB1; TAB2];
LOCITAB(1:end,:) = [];

%=====================================================
% clearvars -except P ADSP GENO 

for IJ = 1:numel(P.FILES.w)
for ff = 1:numel(P.FILES.w{IJ}.mat)

MAT.DAT = load([P.FILES.w{IJ}.path filesep P.FILES.w{IJ}.mat{ff}]);


LOCIS = MAT.DAT.LOCI;


[~,i]  = sort(LOCIS.TRFISHP);
LOCIS  = LOCIS(i,:);

TABI = LOCIS(1:50,1:8);
TABI.VID = (1:50)';

LOCITAB = [LOCITAB; TABI];


end
end



%==========================================================================
%%
clearvars -except P ADSP GENO LOCITAB



nFiles = size(LOCITAB,1) ./ 50;


[C,ia,ic] = unique(LOCITAB.CHRPOS,'stable');

LOCT = LOCITAB(ia,:);

LOCT.COUNT = accumarray(ic,1);

LOCT.SUMPOS = zeros(size(LOCT,1),1);


for i = 1:size(LOCT,1)

    LOCT.SUMPOS(i) = sum( LOCITAB.VID(ic==i) );

    %LOCT.SUMPOS(i) = sum( LOCITAB.VID(ic==i) );

end

LOCT.SUMPOSADJ = ((nFiles-LOCT.COUNT).*(nFiles+5)) + LOCT.SUMPOS;

LOCT.MEANPOS = LOCT.SUMPOSADJ ./ nFiles;


[~,i]  = sort(LOCT.MEANPOS,'ascend');
LOCT  = LOCT(i,:);


LOCT.MEANCOUNT = LOCT.COUNT ./ nFiles;




