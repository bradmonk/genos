
%==========================================================================
%% STEP-1: LOAD THE DATASET
%==========================================================================
close all; clear; clc; rng('shuffle');
P.home = fileparts(which('GENOS.m')); cd(P.home);
P.P1 = [P.home filesep 'genos_functions'];
P.P2 = [P.P1 filesep 'genos_main_functions'];
addpath(join(string(struct2cell(P)),pathsep,1))
cd(P.home)


% which('GENOSDATA.mat')
% ADSP = load('GENOSDATA.mat');

which('ADSP_MINI3.mat')
ADSP = load('ADSP_MINI3.mat');


clearvars -except P ADSP

%==========================================================================
%%   PUT ADSP_MINI3 INTO ADSP STRUCT
%==========================================================================


ADSP.LOCI.VID = (1:numel(ADSP.LOCI.VID))';

ADSP.GOODCOHORTS = [1 2 6 7 9 10 11 12 13 19 20 23 24];
ADSP.USE_APOE = [22 23 24 33 34 44];



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



%==========================================================================
%%             GET APOE AND TOMM40 LOCI
%==========================================================================
clearvars -except IJ ADSP PHEN LOCI CASE CTRL USNP COHSET...
TRCASE TRCTRL TECASE TECTRL

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
%  
% 
% you never two mutations on the same chromosome:
%     APOE 1	45411941 ( C )	45412079 ( T )




ALZGENES = string(["APOE";"TOMM40"]);

x1 = strcmp(LOCI.GENE,ALZGENES{1});
sum(x1)
x2 = strcmp(LOCI.GENE,ALZGENES{2});
sum(x2)

KEEP = x1 | x2;
VLOCI = LOCI(KEEP,:);
VCASE = CASE(KEEP);
VCTRL = CTRL(KEEP);
VUSNP = USNP(KEEP);






clearvars -except ADSP LOCI CASE CTRL PHEN USNP VLOCI VCASE VCTRL VUSNP


%==========================================================================
%% CREATE PERSON X VARIANT MATRIX
%==========================================================================
clearvars -except ADSP LOCI CASE CTRL PHEN USNP VLOCI VCASE VCTRL VUSNP


SRR = PHEN.SRR;
vMX  = zeros( size(SRR,1) , size(VCASE,1) );

for nn = 1:size(VCASE,1)

    CASES = VCASE{nn};
    CTRLS = VCTRL{nn};

    UNCAS = [VUSNP{nn} (-ones(numel(VUSNP{nn}),1))];
    
    CACOUN = [CASES; CTRLS; UNCAS];
    
    if any(any(CACOUN))
        CACOSRR = CACOUN(:,1);
        CACOHH  = CACOUN(:,2);
        [~,Aj] = ismember(SRR , CACOSRR );
        Af = Aj(Aj>0);
        % UNCALL:-1  HOMREF:0  HETALT:1  HOMALT:2  
        vMX(Aj>0,nn) = CACOHH(Af); 
    end

end


vMX = vMX + 20;   % UNCALL:19  HOMREF:20  HETALT:21  HOMALT:22
vMX(vMX==19) =  0;    % UNCALL
vMX(vMX==20) = -2;    % REF/REF
vMX(vMX==21) =  1;    % REF/ALT
vMX(vMX==22) =  2;    % ALT/ALT


PVTAB = PHEN(:,[1 9 18 3 4 2 6 4 4]);
PVTAB.VARS = vMX;


%% ApoE status is defined by a person's combination of these two SNPs:
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
% 
% 
%     APOE 2/2	45411941 ( REF , REF )	45412079 ( ALT , ALT )
%     APOE 2/3	45411941 ( REF , REF )	45412079 ( ALT , REF )
%     APOE 2/4	45411941 ( REF , ALT )	45412079 ( ALT , REF )
%     APOE 3/3	45411941 ( REF , REF )	45412079 ( REF , REF )
%     APOE 3/4	45411941 ( REF , ALT )	45412079 ( REF , REF )
%     APOE 4/4	45411941 ( ALT , ALT )	45412079 ( REF , REF )
% 
% 
% 
%     APOE 2/2	45411941 ( -2 )	45412079 (  2 )
%     APOE 2/3	45411941 ( -2 )	45412079 (  1 )
%     APOE 2/4	45411941 (  1 )	45412079 (  1 )
%     APOE 3/3	45411941 ( -2 )	45412079 ( -2 )
%     APOE 3/4	45411941 (  1 )	45412079 ( -2 )
%     APOE 4/4	45411941 (  2 )	45412079 ( -2 )
% 
% THE VARIANTS IN THE TABLE ABOVE ARE IN THIS ORDER...
% 190045395714	1366260	19	45395714	T	C	SYN	"TOMM40"
% 190045396144	1366264	19	45396144	C	T	SYN	"TOMM40"
% 190045397229	1366271	19	45397229	G	A	SYN	"TOMM40"
% 190045409167	1366301	19	45409167	C	G	MIS	"APOE"
% 190045411110	1366313	19	45411110	T	C	MIS	"APOE"
% *190045411941	1366323	19	45411941	T	C	MIS	"APOE"
% *190045412079	1366323	19	45412079	C	T	MIS	"APOE"
clc; clearvars -except ADSP LOCI CASE CTRL PHEN USNP VLOCI VCASE VCTRL VUSNP PVTAB




APOE22 = (PVTAB.VARS(:,6) == -2) & (PVTAB.VARS(:,7) ==  2);
APOE23 = (PVTAB.VARS(:,6) == -2) & (PVTAB.VARS(:,7) ==  1);
APOE24 = (PVTAB.VARS(:,6) ==  1) & (PVTAB.VARS(:,7) ==  1);
APOE33 = (PVTAB.VARS(:,6) == -2) & (PVTAB.VARS(:,7) == -2);
APOE34 = (PVTAB.VARS(:,6) ==  1) & (PVTAB.VARS(:,7) == -2);
APOE44 = (PVTAB.VARS(:,6) ==  2) & (PVTAB.VARS(:,7) == -2);


APOES = zeros(size(APOE22));
APOES(APOE22) = 22;
APOES(APOE23) = 23;
APOES(APOE24) = 24;
APOES(APOE33) = 33;
APOES(APOE34) = 34;
APOES(APOE44) = 44;

PVTAB.APOE_2 = APOES;

GOODCOH = [1 2 6 7 9 10 11 12 13 19 20 23 24];

clc
MISMATCH = PVTAB.APOE_1 ~= PVTAB.APOE_2;
sum(MISMATCH)

MMGC = PVTAB(MISMATCH & ismember(PVTAB.COHORTNUM,GOODCOH),:);
MMBC = PVTAB(MISMATCH & ~ismember(PVTAB.COHORTNUM,GOODCOH),:);

MISSING = MMGC.APOE_2 == 0;
sum(MISSING)

MISSING = MMBC.APOE_2 == 0;
sum(MISSING)



numel(MMGC.APOE_1( (MMGC.APOE_2 == 0) & (MMGC.AD==1)))
numel(MMBC.APOE_1( (MMBC.APOE_2 == 0) & (MMBC.AD==1)))
numel(MMGC.APOE_1( (MMGC.APOE_2 == 0) & (MMGC.AD==0)))
numel(MMBC.APOE_1( (MMBC.APOE_2 == 0) & (MMBC.AD==0)))

numel(MMGC.APOE_1( (MMGC.APOE_2 ~= 0) & (MMGC.AD==1)))
numel(MMBC.APOE_1( (MMBC.APOE_2 ~= 0) & (MMBC.AD==1)))
numel(MMGC.APOE_1( (MMGC.APOE_2 ~= 0) & (MMGC.AD==0)))
numel(MMBC.APOE_1( (MMBC.APOE_2 ~= 0) & (MMBC.AD==0)))



subplot(2,2,1); histogram(MMGC.APOE_1( (MMGC.APOE_2 == 0) & (MMGC.AD==1)))
title('MISSING & CASE & GOOD COHORT'); xlim([21 45]);

subplot(2,2,2); histogram(MMBC.APOE_1( (MMBC.APOE_2 == 0) & (MMBC.AD==1)))
title('MISSING & CASE & BAD COHORT'); xlim([21 45]);


subplot(2,2,3); histogram(MMGC.APOE_1( (MMGC.APOE_2 == 0) & (MMGC.AD==0)))
title('MISSING & CTRL & GOOD COHORT'); xlim([21 45]);

subplot(2,2,4); histogram(MMBC.APOE_1( (MMBC.APOE_2 == 0) & (MMBC.AD==0)))
title('MISSING & CTRL & BAD COHORT'); xlim([21 45]);



subplot(2,2,1); histogram(MMGC.APOE_1( (MMGC.APOE_2 ~= 0) & (MMGC.AD==1)))
title('MISMATCH & CASE & GOOD COHORT'); xlim([21 45]);

subplot(2,2,2); histogram(MMBC.APOE_1( (MMBC.APOE_2 ~= 0) & (MMBC.AD==1)))
title('MISMATCH & CASE & BAD COHORT'); xlim([21 45])

subplot(2,2,3); histogram(MMGC.APOE_1( (MMGC.APOE_2 ~= 0) & (MMGC.AD==0)))
title('MISMATCH & CTRL & GOOD COHORT'); xlim([21 45])

subplot(2,2,4); histogram(MMBC.APOE_1( (MMBC.APOE_2 ~= 0) & (MMBC.AD==0)))
title('MISMATCH & CTRL & BAD COHORT'); xlim([21 45])


   

%==========================================================================
%%      MAKE  RECTANGLE  NEURAL NET  VARIANT MATRIX
%==========================================================================
%     APOE 2/2	45411941 ( REF , REF )	45412079 ( ALT , ALT )
%     APOE 2/3	45411941 ( REF , REF )	45412079 ( ALT , REF )
%     APOE 2/4	45411941 ( REF , ALT )	45412079 ( ALT , REF )
%     APOE 3/3	45411941 ( REF , REF )	45412079 ( REF , REF )
%     APOE 3/4	45411941 ( REF , ALT )	45412079 ( REF , REF )
%     APOE 4/4	45411941 ( ALT , ALT )	45412079 ( REF , REF )
% 
%     APOE 2/2	45411941 ( -2 )	45412079 (  2 )
%     APOE 2/3	45411941 ( -2 )	45412079 (  1 )
%     APOE 2/4	45411941 (  1 )	45412079 (  1 )
%     APOE 3/3	45411941 ( -2 )	45412079 ( -2 )
%     APOE 3/4	45411941 (  1 )	45412079 ( -2 )
%     APOE 4/4	45411941 (  2 )	45412079 ( -2 )
clc; clearvars -except ADSP LOCI CASE CTRL PHEN USNP VLOCI VCASE VCTRL VUSNP PVTAB

LOCI = ADSP.LOCI;
CASE = ADSP.CASE;
CTRL = ADSP.CTRL;
USNP = ADSP.USNP;
PHEN = ADSP.PHEN;


APO = [22 23 24 33 34 44];

% 190045411941	85852	19	45411941
% 190045412079	85853	19	45412079

APOE22CASE = PVTAB((PVTAB.APOE == 22) & (PVTAB.AD == 1), :);
APOE22CTRL = PVTAB((PVTAB.APOE == 22) & (PVTAB.AD == 0), :);
APOE23CASE = PVTAB((PVTAB.APOE == 23) & (PVTAB.AD == 1), :);
APOE23CTRL = PVTAB((PVTAB.APOE == 23) & (PVTAB.AD == 0), :);
APOE24CASE = PVTAB((PVTAB.APOE == 24) & (PVTAB.AD == 1), :);
APOE24CTRL = PVTAB((PVTAB.APOE == 24) & (PVTAB.AD == 0), :);
APOE33CASE = PVTAB((PVTAB.APOE == 33) & (PVTAB.AD == 1), :);
APOE33CTRL = PVTAB((PVTAB.APOE == 33) & (PVTAB.AD == 0), :);
APOE34CASE = PVTAB((PVTAB.APOE == 34) & (PVTAB.AD == 1), :);
APOE34CTRL = PVTAB((PVTAB.APOE == 34) & (PVTAB.AD == 0), :);
APOE44CASE = PVTAB((PVTAB.APOE == 44) & (PVTAB.AD == 1), :);
APOE44CTRL = PVTAB((PVTAB.APOE == 44) & (PVTAB.AD == 0), :);

USNP{85852} = [];
USNP{85853} = [];

%% APOE22
%--------------------------------------------------------------------------
%     APOE 2/2	45411941 ( REF , REF )	45412079 ( ALT , ALT )
%     APOE 2/2	{85852} ( ~ )	{85853} (  2 )
% 
% REM [SRR ~] FR CASE{85852} CTRL{85852} 
% ADD [SRR 2] TO CASE{85853} CTRL{85853} 
%----------------------------------
disp('APOE22 BEFORE FIX')
CAi = ismember(CASE{85852}(:,1) , APOE22CASE.SRR); disp(sum(CAi))
COi = ismember(CTRL{85852}(:,1) , APOE22CTRL.SRR); disp(sum(COi))
CAj = ismember(CASE{85853}(:,1) , APOE22CASE.SRR); disp(sum(CAj))
COj = ismember(CTRL{85853}(:,1) , APOE22CTRL.SRR); disp(sum(COj))

CASE{85852}(CAi,:) = [];
CTRL{85852}(COi,:) = [];
CASE{85853}(CAj,:) = [];
CTRL{85853}(COj,:) = [];

CASE{85853}( (end+1):(end+numel(APOE22CASE.SRR)),:) = ...
    [APOE22CASE.SRR zeros(size(APOE22CASE.SRR,1),1)+2   ];

CTRL{85853}( (end+1):(end+numel(APOE22CTRL.SRR)),:) = ...
    [APOE22CTRL.SRR zeros(size(APOE22CTRL.SRR,1),1)+2   ];


disp('APOE22 AFTER FIX')
CAi = ismember(CASE{85852}(:,1) , APOE22CASE.SRR); disp(sum(CAi))
COi = ismember(CTRL{85852}(:,1) , APOE22CTRL.SRR); disp(sum(COi))
CAj = ismember(CASE{85853}(:,1) , APOE22CASE.SRR); disp(sum(CAj))
COj = ismember(CTRL{85853}(:,1) , APOE22CTRL.SRR); disp(sum(COj))
%--------------------------------------------------------------------------


%% APOE23
%--------------------------------------------------------------------------
%     APOE 2/2	45411941 ( REF , REF )	45412079 ( ALT , ALT )
%     APOE 2/2	{85852} ( ~ )	{85853} (  2 )
% 
% REM [SRR ~] FR CASE{85852} CTRL{85852} 
% ADD [SRR 2] TO CASE{85853} CTRL{85853} 
%----------------------------------
disp('APOE23 BEFORE FIX')
CAi = ismember(CASE{85852}(:,1) , APOE23CASE.SRR); disp(sum(CAi))
COi = ismember(CTRL{85852}(:,1) , APOE23CTRL.SRR); disp(sum(COi))
CAj = ismember(CASE{85853}(:,1) , APOE23CASE.SRR); disp(sum(CAj))
COj = ismember(CTRL{85853}(:,1) , APOE23CTRL.SRR); disp(sum(COj))

CASE{85852}(CAi,:) = [];
CTRL{85852}(COi,:) = [];
CASE{85853}(CAj,:) = [];
CTRL{85853}(COj,:) = [];

CASE{85853}( (end+1):(end+numel(APOE23CASE.SRR)),:) = ...
    [APOE23CASE.SRR zeros(size(APOE23CASE.SRR,1),1)+1   ];

CTRL{85853}( (end+1):(end+numel(APOE23CTRL.SRR)),:) = ...
    [APOE23CTRL.SRR zeros(size(APOE23CTRL.SRR,1),1)+1   ];


disp('APOE23 AFTER FIX')
CAi = ismember(CASE{85852}(:,1) , APOE23CASE.SRR); disp(sum(CAi))
COi = ismember(CTRL{85852}(:,1) , APOE23CTRL.SRR); disp(sum(COi))
CAj = ismember(CASE{85853}(:,1) , APOE23CASE.SRR); disp(sum(CAj))
COj = ismember(CTRL{85853}(:,1) , APOE23CTRL.SRR); disp(sum(COj))
%--------------------------------------------------------------------------


%% APOE24
%--------------------------------------------------------------------------
%     APOE 2/2	45411941 ( REF , REF )	45412079 ( ALT , ALT )
%     APOE 2/2	{85852} ( ~ )	{85853} (  2 )
% 
% REM [SRR ~] FR CASE{85852} CTRL{85852} 
% ADD [SRR 2] TO CASE{85853} CTRL{85853} 
%----------------------------------
disp('APOE24 BEFORE FIX')
CAi = ismember(CASE{85852}(:,1) , APOE24CASE.SRR); disp(sum(CAi))
COi = ismember(CTRL{85852}(:,1) , APOE24CTRL.SRR); disp(sum(COi))
CAj = ismember(CASE{85853}(:,1) , APOE24CASE.SRR); disp(sum(CAj))
COj = ismember(CTRL{85853}(:,1) , APOE24CTRL.SRR); disp(sum(COj))

CASE{85852}(CAi,:) = [];
CTRL{85852}(COi,:) = [];
CASE{85853}(CAj,:) = [];
CTRL{85853}(COj,:) = [];

CASE{85852}( (end+1):(end+numel(APOE24CASE.SRR)),:) = ...
    [APOE24CASE.SRR zeros(size(APOE24CASE.SRR,1),1)+1   ];

CTRL{85852}( (end+1):(end+numel(APOE24CTRL.SRR)),:) = ...
    [APOE24CTRL.SRR zeros(size(APOE24CTRL.SRR,1),1)+1   ];

CASE{85853}( (end+1):(end+numel(APOE24CASE.SRR)),:) = ...
    [APOE24CASE.SRR zeros(size(APOE24CASE.SRR,1),1)+1   ];

CTRL{85853}( (end+1):(end+numel(APOE24CTRL.SRR)),:) = ...
    [APOE24CTRL.SRR zeros(size(APOE24CTRL.SRR,1),1)+1   ];


disp('APOE24 AFTER FIX')
CAi = ismember(CASE{85852}(:,1) , APOE24CASE.SRR); disp(sum(CAi))
COi = ismember(CTRL{85852}(:,1) , APOE24CTRL.SRR); disp(sum(COi))
CAj = ismember(CASE{85853}(:,1) , APOE24CASE.SRR); disp(sum(CAj))
COj = ismember(CTRL{85853}(:,1) , APOE24CTRL.SRR); disp(sum(COj))
%--------------------------------------------------------------------------



%% APOE33
%--------------------------------------------------------------------------
%     APOE 2/2	45411941 ( REF , REF )	45412079 ( ALT , ALT )
%     APOE 2/2	{85852} ( ~ )	{85853} (  2 )
% 
% REM [SRR ~] FR CASE{85852} CTRL{85852} 
% ADD [SRR 2] TO CASE{85853} CTRL{85853} 
%----------------------------------
disp('APOE33 BEFORE FIX')
CAi = ismember(CASE{85852}(:,1) , APOE33CASE.SRR); disp(sum(CAi))
COi = ismember(CTRL{85852}(:,1) , APOE33CTRL.SRR); disp(sum(COi))
CAj = ismember(CASE{85853}(:,1) , APOE33CASE.SRR); disp(sum(CAj))
COj = ismember(CTRL{85853}(:,1) , APOE33CTRL.SRR); disp(sum(COj))

CASE{85852}(CAi,:) = [];
CTRL{85852}(COi,:) = [];
CASE{85853}(CAj,:) = [];
CTRL{85853}(COj,:) = [];


disp('APOE33 AFTER FIX')
CAi = ismember(CASE{85852}(:,1) , APOE33CASE.SRR); disp(sum(CAi))
COi = ismember(CTRL{85852}(:,1) , APOE33CTRL.SRR); disp(sum(COi))
CAj = ismember(CASE{85853}(:,1) , APOE33CASE.SRR); disp(sum(CAj))
COj = ismember(CTRL{85853}(:,1) , APOE33CTRL.SRR); disp(sum(COj))
%--------------------------------------------------------------------------


%% APOE34
%--------------------------------------------------------------------------
%     APOE 2/2	45411941 ( REF , REF )	45412079 ( ALT , ALT )
%     APOE 2/2	{85852} ( ~ )	{85853} (  2 )
% 
% REM [SRR ~] FR CASE{85852} CTRL{85852} 
% ADD [SRR 2] TO CASE{85853} CTRL{85853} 
%----------------------------------
disp('APOE34 BEFORE FIX')
CAi = ismember(CASE{85852}(:,1) , APOE34CASE.SRR); disp(sum(CAi))
COi = ismember(CTRL{85852}(:,1) , APOE34CTRL.SRR); disp(sum(COi))
CAj = ismember(CASE{85853}(:,1) , APOE34CASE.SRR); disp(sum(CAj))
COj = ismember(CTRL{85853}(:,1) , APOE34CTRL.SRR); disp(sum(COj))

CASE{85852}(CAi,:) = [];
CTRL{85852}(COi,:) = [];
CASE{85853}(CAj,:) = [];
CTRL{85853}(COj,:) = [];

CASE{85852}( (end+1):(end+numel(APOE34CASE.SRR)),:) = ...
    [APOE34CASE.SRR zeros(size(APOE34CASE.SRR,1),1)+1   ];

CTRL{85852}( (end+1):(end+numel(APOE34CTRL.SRR)),:) = ...
    [APOE34CTRL.SRR zeros(size(APOE34CTRL.SRR,1),1)+1   ];


disp('APOE34 AFTER FIX')
CAi = ismember(CASE{85852}(:,1) , APOE34CASE.SRR); disp(sum(CAi))
COi = ismember(CTRL{85852}(:,1) , APOE34CTRL.SRR); disp(sum(COi))
CAj = ismember(CASE{85853}(:,1) , APOE34CASE.SRR); disp(sum(CAj))
COj = ismember(CTRL{85853}(:,1) , APOE34CTRL.SRR); disp(sum(COj))
%--------------------------------------------------------------------------



%% APOE44
%--------------------------------------------------------------------------
%     APOE 2/2	45411941 ( REF , REF )	45412079 ( ALT , ALT )
%     APOE 2/2	{85852} ( ~ )	{85853} (  2 )
% 
% REM [SRR ~] FR CASE{85852} CTRL{85852} 
% ADD [SRR 2] TO CASE{85853} CTRL{85853} 
%----------------------------------
disp('APOE44 BEFORE FIX')
CAi = ismember(CASE{85852}(:,1) , APOE44CASE.SRR); disp(sum(CAi))
COi = ismember(CTRL{85852}(:,1) , APOE44CTRL.SRR); disp(sum(COi))
CAj = ismember(CASE{85853}(:,1) , APOE44CASE.SRR); disp(sum(CAj))
COj = ismember(CTRL{85853}(:,1) , APOE44CTRL.SRR); disp(sum(COj))

CASE{85852}(CAi,:) = [];
CTRL{85852}(COi,:) = [];
CASE{85853}(CAj,:) = [];
CTRL{85853}(COj,:) = [];

CASE{85852}( (end+1):(end+numel(APOE44CASE.SRR)),:) = ...
    [APOE44CASE.SRR zeros(size(APOE44CASE.SRR,1),1)+2   ];

CTRL{85852}( (end+1):(end+numel(APOE44CTRL.SRR)),:) = ...
    [APOE44CTRL.SRR zeros(size(APOE44CTRL.SRR,1),1)+2   ];


disp('APOE44 AFTER FIX')
CAi = ismember(CASE{85852}(:,1) , APOE44CASE.SRR); disp(sum(CAi))
COi = ismember(CTRL{85852}(:,1) , APOE44CTRL.SRR); disp(sum(COi))
CAj = ismember(CASE{85853}(:,1) , APOE44CASE.SRR); disp(sum(CAj))
COj = ismember(CTRL{85853}(:,1) , APOE44CTRL.SRR); disp(sum(COj))
%--------------------------------------------------------------------------

%% SAVE DATASET

clc; clearvars -except LOCI PHEN CASE CTRL USNP

% save('GENOSDATAFINAL.mat')
save('GENOSDATAFINAL_MINI.mat')




