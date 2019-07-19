%% GENOS: 
%==========================================================================
%% STEP-1: LOAD THE DATASET
%==========================================================================
close all; clear; clc; rng('shuffle');
P.home = fileparts(which('GENOS.m')); cd(P.home);
P.funs = [P.home filesep 'genos_functions'];
P.mfuns = [P.funs filesep 'genos_main_functions'];
P.other = [P.home filesep 'genos_other'];
P.data = [P.home filesep 'genos_data'];
% P.data = 'F:\GENOSDATA';
addpath(join(string(struct2cell(P)),pathsep,1))
cd(P.home); P.f = filesep;


ADSP = load('GENOSDATAFINAL.mat');


clc; clearvars -except P ADSP
%==========================================================================
%% (GENOX) PERTURB TEST OPTIONS & FOLDERS 
%==========================================================================
clc; clearvars -except P ADSP


P.serial.DIRroot = [P.home P.f 'genos_data' P.f 'SERIAL' P.f 'RUN3'];
P.serial.DIRmat  = [P.serial.DIRroot P.f 'MAT'];
P.serial.DIRimg  = [P.serial.DIRroot P.f 'IMG'];
% P.serial.DIRstatsmat = [P.serial.DIRroot P.f 'STATS_MAT'];
% P.serial.DIRstatsimg = [P.serial.DIRroot P.f 'STATS_IMG'];



% IMPORT MAT FILES
%------------------------------------------------------
P.LDATA.w = what(P.serial.DIRmat);
P.LDATA.finfo = dir(P.LDATA.w.path);
P.LDATA.finames = {P.LDATA.finfo.name};
c=~cellfun(@isempty,regexp(P.LDATA.finames,'((\S)+(\.mat+))'));
P.LDATA.finames = string(P.LDATA.finames(c)');
P.LDATA.folder = P.LDATA.finfo.folder;
P.LDATA.fipaths = fullfile(P.LDATA.folder,P.LDATA.finames);
disp(P.LDATA.fipaths); disp(P.LDATA.finames);
disp('TOTAL MAT FILES:'); disp(numel(P.LDATA.finames));



% LDATA = load(P.LDATA.fipaths{1});



clearvars -except P ADSP
%==========================================================================
%% GET SERIAL STATS
%==========================================================================
clc; clearvars -except P ADSP


DATA = load(P.LDATA.fipaths{1});

LOOP.GENECHRPOS  = DATA.LOOPDATA.GENECHRPOS;
LOOP.XLOCI       = DATA.LOOPDATA.XLOCI;
LOOP.NETS        = DATA.LOOPDATA.NETS;

LOOP.PVXt   = DATA.LOOPDATA.PVXt;
LOOP.PVXh   = DATA.LOOPDATA.PVXh;

LOOP.NNREFt = DATA.LOOPDATA.NNREFt;
LOOP.NNREFh = DATA.LOOPDATA.NNREFh;
LOOP.NNALTt = DATA.LOOPDATA.NNALTt;
LOOP.NNALTh = DATA.LOOPDATA.NNALTh;


clc; clearvars -except P ADSP LOOP
%==========================================================================


























%==========================================================================
%% BUILD LOCI TABLE
%==========================================================================
clc; clearvars -except P ADSP LOOP


% GET LOOP PARAMETERS FROM DATA
%----------------------------------------
XLOCIS = LOOP.XLOCI';
XLOCI = XLOCIS{1};

P.SnpSets     = numel(P.LDATA.fipaths);
P.LoopsPerSet = size(XLOCIS,1);
P.nTopTarg    = size(XLOCI,1);
P.nTarg       = size(XLOCI,1) - 50;
P.nSNP        = P.SnpSets * P.nTarg;

fprintf('N LOOPS PER TARGET SET: %3.f \n',P.LoopsPerSet);
fprintf('N TARGETS + TOPSNPS: %6.f \n',P.nTopTarg);
fprintf('N TARGETS: %16.f \n',P.nTarg);



%==========================================================================
for i = 1:P.SnpSets
%==========================================================================
disp(i);


% LOAD DATASET FROM MAT FILE
DATA    = load(P.LDATA.fipaths{i});
XLOCIS  = DATA.LOOPDATA.XLOCI';




% EXTRACT VALUES FROM EACH LOOP
%----------------------------------------
for j = 1:P.LoopsPerSet

    XLOCI         = XLOCIS{j};

    % SET TRAINING FISHP BACK TO REAL VALUE FROM FISHP
    XLOCI.TRFISHP    = XLOCI.FISHP;
    XLOCI.FISHP      = XLOCI.OKFISHP;
    XLOCI.FISHOR     = XLOCI.OKFISHOR;
    XLOCI.CASEREF    = XLOCI.OKCASEREF;
    XLOCI.CASEALT    = XLOCI.OKCASEALT;
    XLOCI.CTRLREF    = XLOCI.OKCTRLREF;
    XLOCI.CTRLALT    = XLOCI.OKCTRLALT;

    XLOCI         = XLOCI(1:P.nTarg,[1 3:8 13:18 25:32 41:44]);
    XLOCI.TARGET  = (1:size(XLOCI,1))';
    XLOCI.SET     = (zeros(size(XLOCI,1),1) + i);
    XLOCI.srtVAL  = log(abs(-log(XLOCI.FISHP)))+abs(-log(XLOCI.FISHOR));

    if j == 1
        jLOC = XLOCI;
    else
        if any(~strcmp(jLOC.GENE,XLOCI.GENE)); keyboard; end
        jLOC.TRCASEREF = jLOC.TRCASEREF + XLOCI.TRCASEREF;
        jLOC.TRCTRLREF = jLOC.TRCTRLREF + XLOCI.TRCTRLREF;
        jLOC.TRCASEALT = jLOC.TRCASEALT + XLOCI.TRCASEALT;
        jLOC.TRCASEALT = jLOC.TRCASEALT + XLOCI.TRCASEALT;
        jLOC.TECASEREF = jLOC.TECASEREF + XLOCI.TECASEREF;
        jLOC.TECTRLREF = jLOC.TECTRLREF + XLOCI.TECTRLREF;
        jLOC.TECASEALT = jLOC.TECASEALT + XLOCI.TECASEALT;
        jLOC.TECASEALT = jLOC.TECASEALT + XLOCI.TECASEALT;
        jLOC.TRFISHP   = jLOC.TRFISHP   + XLOCI.TRFISHP;
        jLOC.TEFISHP   = jLOC.TEFISHP   + XLOCI.TEFISHP;
        jLOC.TRFISHOR  = jLOC.TRFISHOR  + XLOCI.TRFISHOR;
        jLOC.TEFISHOR  = jLOC.TEFISHOR  + XLOCI.TEFISHOR;
    end

end

    jLOC.TRCASEREF = round(jLOC.TRCASEREF ./ P.LoopsPerSet);
    jLOC.TRCTRLREF = round(jLOC.TRCTRLREF ./ P.LoopsPerSet);
    jLOC.TRCASEALT = round(jLOC.TRCASEALT ./ P.LoopsPerSet);
    jLOC.TRCASEALT = round(jLOC.TRCASEALT ./ P.LoopsPerSet);
    jLOC.TECASEREF = round(jLOC.TECASEREF ./ P.LoopsPerSet);
    jLOC.TECTRLREF = round(jLOC.TECTRLREF ./ P.LoopsPerSet);
    jLOC.TECASEALT = round(jLOC.TECASEALT ./ P.LoopsPerSet);
    jLOC.TECASEALT = round(jLOC.TECASEALT ./ P.LoopsPerSet);
    jLOC.TRFISHP   = jLOC.TRFISHP   ./ P.LoopsPerSet;
    jLOC.TEFISHP   = jLOC.TEFISHP   ./ P.LoopsPerSet;
    jLOC.TRFISHOR  = jLOC.TRFISHOR  ./ P.LoopsPerSet;
    jLOC.TEFISHOR  = jLOC.TEFISHOR  ./ P.LoopsPerSet;

    if i == 1
    iLOC = jLOC;
    else
    iLOC = [iLOC; jLOC];
    end
end
%==========================================================================

writetable(iLOC,'/Users/bradleymonk/Desktop/iLOC.csv');

return






























%==========================================================================
%% GET SERIAL STATS
%==========================================================================
clc; clearvars -except P ADSP


DATA = load(P.LDATA.fipaths{1});

LOOP.GENECHRPOS  = DATA.LOOPDATA.GENECHRPOS;
LOOP.XLOCI       = DATA.LOOPDATA.XLOCI;
LOOP.NETS        = DATA.LOOPDATA.NETS;

LOOP.PVXt   = DATA.LOOPDATA.PVXt;
LOOP.PVXh   = DATA.LOOPDATA.PVXh;

LOOP.NNREFt = DATA.LOOPDATA.NNREFt;
LOOP.NNREFh = DATA.LOOPDATA.NNREFh;
LOOP.NNALTt = DATA.LOOPDATA.NNALTt;
LOOP.NNALTh = DATA.LOOPDATA.NNALTh;


clc; clearvars -except P ADSP LOOP
%==========================================================================




%==========================================================================
%% GET SERIAL STATS
%==========================================================================
clc; clearvars -except P ADSP LOOP


% GET LOOP PARAMETERS FROM DATA
%----------------------------------------
XLOCIS = LOOP.XLOCI';
XLOCI = XLOCIS{1};

P.SnpSets     = numel(P.LDATA.fipaths);
P.LoopsPerSet = size(XLOCIS,1);
P.nTopTarg    = size(XLOCI,1);
P.nTarg       = size(XLOCI,1) - 50;
P.nSNP        = P.SnpSets * P.nTarg;

fprintf('N LOOPS PER TARGET SET: %3.f \n',P.LoopsPerSet);
fprintf('N TARGETS + TOPSNPS: %6.f \n',P.nTopTarg);
fprintf('N TARGETS: %16.f \n',P.nTarg);




















%==========================================================================
% GET FULL PHENOTYPE TO CAT WITH REF&ALT EDITED TARGET NN ACTIVITY
% DUPLICATE FOR REF & ALT
% ADD 1001 COLUMNS
%----------------------------------------
PHEN  = ADSP.PHEN;
PHEN(:,11:17) = [];
PHEN  = table2array(PHEN);
PHREF = [PHEN  zeros(size(PHEN,1),P.nSNP)];
PHALT = [PHEN  zeros(size(PHEN,1),P.nSNP)];
PHNAT = [PHEN  zeros(size(PHEN,1),P.nSNP)];
PHREF = repmat(PHREF,1,1,5);
PHALT = repmat(PHALT,1,1,5);
PHNAT = repmat(PHNAT,1,1,5);






% GET FULL PHENOTYPE TABLE TO CAT WITH TARGET GENOTYPES
% ADD 1001 COLUMNS
%----------------------------------------
PVMX  = [PHEN  zeros(size(PHEN,1),P.nSNP)];
PVMX  = repmat(PVMX,1,1,5);


% GET PARAMETERS
%----------------------------------------
S     = size(PHEN,2)+1;
N     = P.nTarg;
i=1;
j=1;


%%
%==========================================================================
for i = 1:P.SnpSets
%==========================================================================
clc; clearvars -except i j S N P ADSP LOOP PHEN PHREF PHALT PHNAT PVMX LOXI
disp(P.LDATA.finames{i}); 
disp(LOOP.GENECHRPOS{1}(1:P.nTarg)')



% GET DATA FROM MAT FILE
%----------------------------------------
DATA = load(P.LDATA.fipaths{i});

LOOP.GENECHRPOS = DATA.LOOPDATA.GENECHRPOS;
LOOP.XLOCI      = DATA.LOOPDATA.XLOCI;
LOOP.NETS       = DATA.LOOPDATA.NETS;

LOOP.PVXt       = DATA.LOOPDATA.PVXt;
LOOP.PVXh       = DATA.LOOPDATA.PVXh;

LOOP.NNREFt     = DATA.LOOPDATA.NNREFt;
LOOP.NNREFh     = DATA.LOOPDATA.NNREFh;
LOOP.NNALTt     = DATA.LOOPDATA.NNALTt;
LOOP.NNALTh     = DATA.LOOPDATA.NNALTh;



% GET SNP LOCI TABLE; ADD TARGET & SET COLUMNS
LOX           = LOOP.XLOCI{1};
LOX.TARGET    = 1:size(LOX,1);
LOX.SET       = zeros(size(LOX,1),1) + i;
% LOX(11:end,:) = [];
% LOX(:,2:12)   = [];
% LOX.CHRPOS    = double(LOX.CHRPOS);
% LOX           = table2array(LOX);
% disp(LOX(1:5,1:5)); disp(' '); pause(2)





% EXTRACT VALUES FROM 1ST LOOP
%----------------------------------------
for j = 1:P.LoopsPerSet



% GET NN MATRICES & NN OUTPUTS
PVX   = [LOOP.PVXt(:,:,j); LOOP.PVXh(:,:,j)];
inan  = isnan(PVX(:,1)); PVX(inan,:)   = [];

REF   = [LOOP.NNREFt(:,:,j); LOOP.NNREFh(:,:,j)];
inan  = isnan(REF(:,1)); REF(inan,:)   = [];

ALT   = [LOOP.NNALTt(:,:,j); LOOP.NNALTh(:,:,j)];
inan  = isnan(ALT(:,1)); ALT(inan,:)   = [];



% DETERMINE ORDER OF PARTICIPANTS IN PVX REF ALT
[~,ai] = ismember( PHEN(:,1,1), PVX(:,1) );  aj = ai(ai>0);
[~,bi] = ismember( PHEN(:,1,1), REF(:,1) );  bj = bi(bi>0);
[~,ci] = ismember( PHEN(:,1,1), ALT(:,1) );  cj = ci(ci>0);


% ADD SNPS TO THE GENOTYPE & PHENOTYPE TABLES
kPVX(ai>0, (i*N-N+S):(i*N+S-1)  ,j) = PVX(aj,10:(10+N-1));
kREF(bi>0, (i*N-N+S):(i*N+S-1)  ,j) = REF(bj,10:(10+N-1));
kALT(ci>0, (i*N-N+S):(i*N+S-1)  ,j) = ALT(cj,10:(10+N-1));


% ADD NN OUTPUT FOR THE NON-EDIT-TARGETS (1 COL) TO END OF PHENOTYPE TABLE
NX = REF(bj,end);
NX = repmat(NX,1,N);
kNAT(bi>0, (i*N-N+S):(i*N+S-1)  ,j) = NX;



% ADD SNPS TO LOXI TABLE
kLOX           = LOOP.XLOCI{j};
kLOX.TARGET    = 1:size(LOX,1);
kLOX.SET       = zeros(size(LOX,1),1) + i;



jLOX(:,:,j) = kLOX;
jPVX(:,:,j) = kPVX;
jREF(:,:,j) = kREF;
jALT(:,:,j) = kALT;
jNAT(:,:,j) = kNAT;
end



r=(i*N-N+1):(i*N);
iLOX(r,:,j) = mean(jLOX,3);
iPVX(r,:,j) = mean(jPVX,3);
iREF(r,:,j) = mean(jREF,3);
iALT(r,:,j) = mean(jALT,3);
iNAT(r,:,j) = mean(jNAT,3);

end
%==========================================================================
%% POST-PROCESS AVERAGED DATA
%==========================================================================
clc; clearvars -except i j S N P ADSP LOOP PHEN PHREF PHALT PHNAT PVMX LOXI


% GET MEAN OVER THE 5 LOOPS FOR EACH SNP
LOXm = mean(LOXI,3);
PVXm = mean(PVMX,3);
REFm = mean(PHREF,3);
ALTm = mean(PHALT,3);
NATm = mean(PHNAT,3);


%==========================================================================
PHEN = ADSP.PHEN;
PHEN(2:end,:)   = [];
PHEN(:,11:17)   = [];

REFmu = [PHEN  array2table(zeros(size(PHEN,1),1000))];
ALTmu = [PHEN  array2table(zeros(size(PHEN,1),1000))];
NATmu = [PHEN  array2table(zeros(size(PHEN,1),1000))];
DIFmu = [PHEN  array2table(zeros(size(PHEN,1),1000))];

REFmu(1:size(REFm,1),:) = array2table(REFm);
ALTmu(1:size(ALTm,1),:) = array2table(ALTm);
NATmu(1:size(NATm,1),:) = array2table(NATm);
DIFm = ALTm;
DIFm(:,16:end) = ALTm(:,16:end) - REFm(:,16:end);
DIFmu(1:size(ALTm,1),:) = array2table( DIFm );

writetable(REFmu,'/Users/bradleymonk/Desktop/SERIALSNPS_NETREF.csv');
writetable(ALTmu,'/Users/bradleymonk/Desktop/SERIALSNPS_NETALT.csv');
writetable(DIFmu,'/Users/bradleymonk/Desktop/SERIALSNPS_NETDIF.csv');
writetable(NATmu,'/Users/bradleymonk/Desktop/SERIALSNPS_NETNAT.csv');
%==========================================================================







LOXm = mean(LOXI,3);
PVXm = mean(PVMX,3);

DNAmu = [PHEN  array2table(zeros(size(PHEN,1),1000))];
PVXmu = [PHEN  array2table(zeros(size(PHEN,1),1000))];


%==========================================================================

LOCI  = ADSP.LOCI;


LOXmu           = LOOP.XLOCI{1};
LOXmu.TARGET    = (1:size(LOXmu,1))';
LOXmu.SET       = zeros(size(LOXmu,1),1);
LOXmu(2:end,:)  = [];
LOXmu(:,2:12)   = [];
LOXmu(1:size(LOXm,1),:) = array2table(LOXm);





[C,ia,ic] = unique(LOXmu.CHRPOS,'stable');
LOXmu(ia,:);

LOCI  = ADSP.LOCI;
[Ai,Bi] = ismember( LOXmu.CHRPOS  ,  LOCI.CHRPOS);
LOCI  = LOCI(a,:);

a
numel(a)
sum(a)




SNPmu = [LOXmu(:,1) LOCI(:,[3:8])  LOXmu(:,[8:13 2:7 18:21 31:35])];


PVXmu(1:size(PVXm,1),:) = array2table(PVXm);

% ij=sum((LOCI.CHRPOS == LOXmu.CHRPOS')>0,2)>0;
% LOCI = LOCI(ij,:);
% sum(ij)
% writetable(LOCI ,'/Users/bradleymonk/Desktop/SERIALSNPS_LOCI.csv');
% writetable(LOXmu,'/Users/bradleymonk/Desktop/SERIALSNPS_LOXmu.csv');








% Ad = abs(ALTm(:,16:end) - NATm(:,16:end));
% Rd = abs(REFm(:,16:end) - NATm(:,16:end));
% AdRd = Ad;
% AdRd(:,:,2) = Rd;
% [MaxAdRd, Si] = max(abs(AdRd),[],3);
% DNAm = NATmu;
% DNAm(1:size(MaxAdRd,1),:) = array2table( MaxAdRd );







return
%==========================================================================
% GET GENE-CHRLOC TO PERTURB AND ENSURE IT'S #1 AMONG SNPS
%==========================================================================
% GET LOCI TABLE 


LOXmu.TARGET = repmat( [1:10]'  ,100 ,1);
set = repmat( [1:10]  ,100 ,1);
LOXmu.SET = set(:);
LOXmu.ORDER = (1:1000)';



% LOCI = ADSP.LOCI;
% ij=sum((LOCI.CHRPOS == LOXmu.CHRPOS')>0,2)>0;
% LOCI = LOCI(ij,:);
% sum(ij)

% LOCI  = sortrows(LOCI,'CHRPOS');
% LOXmu = sortrows(LOXmu,'CHRPOS');
% SNPmu = sortrows(SNPmu,'ORDER');

SNPmu = [LOXmu(:,1) LOCI(:,[3:8])  LOXmu(:,[8:13 2:7 18:21 32:36])];


%==========================================================================
%% EXPORT VALUES TO TABLES
%==========================================================================


writetable(PVXmu,'/Users/bradleymonk/Desktop/SERIALSNPS_GENO.csv');
writetable(SNPmu,'/Users/bradleymonk/Desktop/SERIALSNPS_LOCI.csv');

writetable(REFmu,'/Users/bradleymonk/Desktop/SERIALSNPS_NETREF.csv');
writetable(ALTmu,'/Users/bradleymonk/Desktop/SERIALSNPS_NETALT.csv');
writetable(DIFmu,'/Users/bradleymonk/Desktop/SERIALSNPS_NETDIF.csv');
writetable(NATmu,'/Users/bradleymonk/Desktop/SERIALSNPS_NETNAT.csv');


%==========================================================================
%%
