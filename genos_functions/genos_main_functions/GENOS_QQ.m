% GENETIC VARIANT PROCESSING PIPELINE TO PREP FOR NEURAL NET TRAINING
% DATASET: ADSP WES VCF (LATEST RELEASE)
%
% Overview of goals, issues, and processing steps:
%{
The neural network classifier requires a 2D rectangular matrix for 
training. This PERSON x LOCI matrix, will indicate for each person 
at each locus (chr:pos) whether they are homozygous reference ('0'),
heterozygous ('1'), or homozygous alternate ('2').


The pipeline will prepare to matrices, each composed of a unique group
of people. This will provide one dataset for neural net training, and
another for testing the accuracy of the classifier.


Instead of listing every processing detail here, those will be explained
in each block of code along the way. Though, I will explain a data
stratification issue we have previously identified so it's clear where
things are going out-of-the-gate. 


ADSP has combined sequencing data from 25 different consortium studies
(wherever you see the code mention 'STUDY' or 'COHORT', assume they are 
all synonymously referring to these 25 data-sources). Because of this, 
major stratification issues lurk in this dataset. Further, many of these 
PHEN provide very unequal numbers of Alzheimer's patients (CASE) and 
controls (CTRL).


Here is a summary table of those 25 STUDY-PHEN, and their 
share of participant contributions...

%% COID CONSORT	STUDY_NAME  STUDY   FAM PEOPLE  ENR     CASE	CTRL	TOT    %CA %CO  EQL	CULL	NOTES
01  DGC     Adult Chng  ACT     0	1277	0       323     945     1268	25	75	323	-622	KEEP
02  ADGC	AD Centers  ADC     0	3263	0       2438	817     3255	75	25	817	1621	KEEP
03  CHARGE	Athrosclro  ARIC	0	57      0       39      18      57      68	32	18	21      REMOVE
04  CHARGE	Aus Stroke	ASPS	0	126     0       121     5       126     96	4	5	116     REMOVE
05  ADGC	Chic Aging  CHAP	0	231     0       27      204     231     12	88	27	-177	REMOVE
06  CHARGE	CardioHlth  CHS     0	840     0       250     583     833     30	70	250	-333	KEEP
07  ADGC	Hispanic S  CUHS	361	0       331     160     171     331     48	52	160	-11     HISPANIC
08  CHARGE	Erasmus Er  ERF     30	45      0       45      0       45      100	0	0	45      REMOVE
09  CHARGE	Framingham  FHS     0	581     0       157     424     581     27	73	157	-267	KEEP
10	ADGC	Gene Diffs  GDF     0	207     0       111     96      207     54	46	96	15      KEEP
11	ADGC	NIA - LOAD  LOAD	80  113     363     367     109     476     77	23	109	258     KEEP
12	ADGC	Aging Proj  MAP     0	415     0       138     277     415     33	67	138	-139	KEEP
13	ADGC	Mayo Clini  MAYO	0	349     0       250     99      349     72	28	99	151     KEEP
14	ADGC	Miami Univ  MIA     61	181     19      186     14      200     93	7	14	172     REMOVE
15	ADGC	AD Genetic  MIR     0	284     47      316     15      331     95	5	15	301     REMOVE
16	ADGC	Mayo cl PD  MPD     0	20      0       0       20      20      0	100	0	-20     REMOVE
17	ADGC	NationC AD  NCRD	18	108     52      160     0       160     100	0	0	160     REMOVE
18	ADGC	Wash Unive  RAS     0	0       46      46      0       46      100	0	0	46      REMOVE
19	ADGC	Relig Ordr  ROS     0	351     0       154     197     351     44	56	154	-43     KEEP
20	CHARGE	RotterdamS  RS      0	1089	0       276     813     1089	25	75	276	-537	KEEP
21	ADGC	Texas AD S  TARC	0	144     0       132     12      144     92	8	12	120     REMOVE
22	ADGC	Un Toronto  TOR     0	0       9       9       0       9       100	0	0	9       REMOVE
23	ADGC	Vanderbilt  VAN     6	235     1       210     26      236     89	11	26	184     REMOVE
24	ADGC	WashNY Age  WCAP	0	147     3       34      116     150     23	77	34	-82     REMOVE
25	ADGC 	Univ Penns  UPN     40 	0       3       0       0       0       0	0	0	0       REMOVE



Along the way, the code will take steps to mitigate stratification
issues due to COHORT-specific effects. 

...I'm trying to think if there is anything else you should know...
I can't, sooo, I guess we should just go-ahead and start with step #1.
(aka the first chance for something to go wrong)...

%}

%% STEP-1: LOAD THE DATASET

clc; close all; clear; rng('shuffle')
cd(fileparts(which('GENOS_QQ.m')));


MATDATA = 'ADSP_ALT_UNC.mat';
which(MATDATA)
load(MATDATA)


clc; rng('shuffle');
clearvars -except ADSP




%% CARBON COPY MAIN VARIABLES FROM ADSP.STRUCT

LOCI = ADSP.LOCI;
CASE = ADSP.CASE;
CTRL = ADSP.CTRL;
PHEN = ADSP.PHEN;

clearvars -except ADSP LOCI CASE CTRL PHEN





%% MAKE RAND GROUP BASED ON SHUFFLED DATASET


[RCASE, RCTRL, RPHEN] = makerandgroups(PHEN, CASE, CTRL);


clearvars -except ADSP LOCI CASE CTRL PHEN RCASE RCTRL RPHEN




%###############################################################
%%          COUNT NUMBER OF VARIANTS PER LOCI
%###############################################################

% The varsum() function will go through each known variant loci
% and check whether anyone's SRR ID from your subset of IDs match
% all known SRR IDs for that loci. It will then sum the total
% number of alleles (+1 for hetzy-alt, +2 for homzy-alt) for each
% loci and return the totals.


% TO CONSERVE CODING EFFORTS, TREAT RAND GROUP AS NN TEST GROUP

NCASESRR = PHEN.SRR(PHEN.AD==1);
NCTRLSRR = PHEN.SRR(PHEN.AD==0);

RCASESRR = RPHEN.SRR(RPHEN.AD==1);
RCTRLSRR = RPHEN.SRR(RPHEN.AD==0);



[NCASEN, NCTRLN] = varsum( CASE, NCASESRR,  CTRL, NCTRLSRR);
[RCASEN, RCTRLN] = varsum(RCASE, RCASESRR, RCTRL, RCTRLSRR);




USNP = ADSP.UNSNP;

[TRCASEUN, TRCTRLUN, TECASEUN, TECTRLUN] = uncsum(USNP,...
                NCASESRR, NCTRLSRR, RCASESRR, RCTRLSRR);



% SAVE COUNTS AS NEW TABLE COLUMNS
LOCI.TRCASEREF = (numel(NCASESRR)*2) - (TRCASEUN.*2) - NCASEN;
LOCI.TRCTRLREF = (numel(NCTRLSRR)*2) - (TRCTRLUN.*2) - NCTRLN;
LOCI.TRCASEALT = NCASEN;
LOCI.TRCTRLALT = NCTRLN;


LOCI.TECASEREF = (numel(RCASESRR)*2) - (TECASEUN.*2) - RCASEN;
LOCI.TECTRLREF = (numel(RCTRLSRR)*2) - (TECTRLUN.*2) - RCTRLN;
LOCI.TECASEALT = RCASEN;
LOCI.TECTRLALT = RCTRLN;



%%
close all
fh01 = figure('Units','normalized','OuterPosition',[.01 .25 .48 .62],'Color','w');
ax01 = axes('Position',[.06 .06 .9 .9],'Color','none');
ph01 = histogram(NCASEN(NCASEN>5), 40,'DisplayStyle','stairs',...
    'LineWidth',3,'Normalization','count'); hold on;
ph02 = histogram(NCTRLN(NCTRLN>5), 40,'DisplayStyle','stairs',...
    'LineWidth',3,'Normalization','count'); hold on;
ph03 = histogram(RCASEN(RCASEN>5), 40,'DisplayStyle','stairs',...
    'LineWidth',3,'Normalization','count'); hold on;
ph04 = histogram(RCTRLN(RCTRLN>5), 40,'DisplayStyle','stairs',...
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



clearvars -except ADSP LOCI CASE CTRL PHEN RCASE RCTRL RPHEN

















%###############################################################
%%               COMPUTE FISHER'S P-VALUE
%###############################################################


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




close all
x = -log(LOCI.TRFISHP);
y = -log(LOCI.TEFISHP);
histogram(x(x>6.9)); hold on
histogram(y(y>6.9));


clearvars -except ADSP LOCI CASE CTRL PHEN RCASE RCTRL RPHEN






















%%
%{

%###############################################################
%%          DISPLAY TOTAL VARIANT COUNTS PER CASE V CTRL
%###############################################################



sum(LOCI.TRCASEALT( strcmp(cellstr(LOCI.EFFECT),'SYN') ))
sum(LOCI.TRCTRLALT( strcmp(cellstr(LOCI.EFFECT),'SYN') ))

sum(LOCI.TRCASEALT( strcmp(cellstr(LOCI.EFFECT),'MIS') ))
sum(LOCI.TRCTRLALT( strcmp(cellstr(LOCI.EFFECT),'MIS') ))



sum(LOCI.TRCASEALT( ~strcmp(cellstr(LOCI.EFFECT),'SYN') &...
                    ~strcmp(cellstr(LOCI.EFFECT),'MIS')   ))

sum(LOCI.TRCTRLALT( ~strcmp(cellstr(LOCI.EFFECT),'SYN') &...
                    ~strcmp(cellstr(LOCI.EFFECT),'MIS')   ))







%% SUM OVER ASSYMETRIC VARIANTS

PASS = ((LOCI.TRCASEALT > ((LOCI.TRCTRLALT+1)*3)) |...
        (LOCI.TRCTRLALT > ((LOCI.TRCASEALT+1)*3)));
sum(PASS)


sum(LOCI.TRCASEALT( strcmp(cellstr(LOCI.EFFECT),'SYN') & PASS ))
sum(LOCI.TRCTRLALT( strcmp(cellstr(LOCI.EFFECT),'SYN') & PASS ))

sum(LOCI.TRCASEALT( strcmp(cellstr(LOCI.EFFECT),'MIS') & PASS ))
sum(LOCI.TRCTRLALT( strcmp(cellstr(LOCI.EFFECT),'MIS') & PASS ))



sum(LOCI.TRCASEALT( (~strcmp(cellstr(LOCI.EFFECT),'SYN')  &...
                     ~strcmp(cellstr(LOCI.EFFECT),'MIS')) &...
                    PASS   ))

sum(LOCI.TRCTRLALT( (~strcmp(cellstr(LOCI.EFFECT),'SYN')  &...
                     ~strcmp(cellstr(LOCI.EFFECT),'MIS')) &...
                     PASS   ))






















%% COMPUTE CHI SQUARE VALUE



% COMPUTE CHI SQUARED STATISTICS FOR THE TRAINING GROUP
[CHIP, CHIOR] = chisq(LOCI.TRCASEREF,LOCI.TRCASEALT,...
                      LOCI.TRCTRLREF,LOCI.TRCTRLALT);

LOCI.TRCHIP  = CHIP;
LOCI.TRCHIOR = CHIOR;



% COMPUTE CHI SQUARED STATISTICS FOR THE TRAINING GROUP
[CHIP, CHIOR] = chisq(LOCI.TECASEREF, LOCI.TECASEALT,...
                      LOCI.TECTRLREF, LOCI.TECTRLALT);

LOCI.TECHIP  = CHIP;
LOCI.TECHIOR = CHIOR;




%###############################################################
%###############################################################
%###############################################################
%###############################################################
%% PREPARE DATA FOR NEURAL NET CLASSIFIER SUPERVISED LEARNING
%###############################################################
%###############################################################
%###############################################################
%###############################################################
clc;




%% STORE VARIABLES FOR NEURAL NETWORK TRAINING AS 'AMX'

AMX         = LOCI;
AMXCASE     = CASE;
AMXCTRL     = CTRL;
AMXTRCASE   = CTRCASE;
AMXTRCTRL   = CTRCTRL;
AMXTECASE   = CTECASE;
AMXTECTRL   = CTECTRL;


clearvars -except ADSP GENB LOCI CASE CTRL PHEN...
CTRCASE CTRCTRL CTECASE CTECTRL TRCASEN TRCTRLN TECASEN TECTRLN...
AMX AMXCASE AMXCTRL AMXTRCASE AMXTRCTRL AMXTECASE AMXTECTRL








%% FILTER VARIANTS BASED ALT > REF

PASS = (AMX.TRCASEREF > AMX.TRCASEALT./1.5) | (AMX.TRCTRLREF > AMX.TRCTRLALT./1.5);
sum(~PASS)

AMX      = AMX(PASS,:);
AMXCASE  = AMXCASE(PASS);
AMXCTRL  = AMXCTRL(PASS);
AMX.VID  = (1:size(AMX,1))';




clearvars -except ADSP GENB LOCI CASE CTRL PHEN...
CTRCASE CTRCTRL CTECASE CTECTRL TRCASEN TRCTRLN TECASEN TECTRLN...
AMX AMXCASE AMXCTRL AMXTRCASE AMXTRCTRL AMXTECASE AMXTECTRL






%% TAG ASYMETRIC VARIANTS


AMX.MORECASE = (AMX.TRCASEALT > ((AMX.TRCTRLALT+1).*2.1));
AMX.MORECTRL = (AMX.TRCTRLALT > ((AMX.TRCASEALT+1).*2.1));

AMX.ASYM = AMX.MORECASE | AMX.MORECTRL;



clearvars -except ADSP GENB LOCI CASE CTRL PHEN...
CTRCASE CTRCTRL CTECASE CTECTRL TRCASEN TRCTRLN TECASEN TECTRLN...
AMX AMXCASE AMXCTRL AMXTRCASE AMXTRCTRL AMXTECASE AMXTECTRL









%###############################################################
%%          COUNT NUMBER OF VARIANTS PER GENE
%###############################################################

% The varsum() function will go through each known variant loci
% and check whether anyone's SRR ID from your subset of IDs match
% all known SRR IDs for that loci. It will then sum the total
% number of alleles (+1 for hetzy-alt, +2 for homzy-alt) for each
% loci and return the totals.


TRCASEREF = AMX.TRCASEREF;
TRCTRLREF = AMX.TRCTRLREF;
TRCASEALT = AMX.TRCASEALT;
TRCTRLALT = AMX.TRCTRLALT;

TECASEREF = AMX.TECASEREF;
TECTRLREF = AMX.TECTRLREF;
TECASEALT = AMX.TECASEALT;
TECTRLALT = AMX.TECTRLALT;

GTRCASEREF = AMX.TRCASEREF;
GTRCTRLREF = AMX.TRCTRLREF;
GTRCASEALT = AMX.TRCASEALT;
GTRCTRLALT = AMX.TRCTRLALT;

GTECASEREF = AMX.TECASEREF;
GTECTRLREF = AMX.TECTRLREF;
GTECASEALT = AMX.TECASEALT;
GTECTRLALT = AMX.TECTRLALT;



u = unique(AMX.GENEi);

for nn = 1:numel(u)

    i = AMX.GENEi == u(nn);


    GTRCASEREF(i) = sum(TRCASEREF(i));
    GTRCTRLREF(i) = sum(TRCTRLREF(i));
    GTRCASEALT(i) = sum(TRCASEALT(i));
    GTRCTRLALT(i) = sum(TRCTRLALT(i));

    GTECASEREF(i) = sum(TECASEREF(i));
    GTECTRLREF(i) = sum(TECTRLREF(i));
    GTECASEALT(i) = sum(TECASEALT(i));
    GTECTRLALT(i) = sum(TECTRLALT(i));

    if ~mod(nn,100); disp(nn/numel(u)); end
end


AMX.GTRCASEREF = GTRCASEREF;
AMX.GTRCTRLREF = GTRCTRLREF;
AMX.GTRCASEALT = GTRCASEALT;
AMX.GTRCTRLALT = GTRCTRLALT;

AMX.GTECASEREF = GTECASEREF;
AMX.GTECTRLREF = GTECTRLREF;
AMX.GTECASEALT = GTECASEALT;
AMX.GTECTRLALT = GTECTRLALT;



clearvars -except ADSP GENB LOCI CASE CTRL PHEN...
CTRCASE CTRCTRL CTECASE CTECTRL TRCASEN TRCTRLN TECASEN TECTRLN...
AMX AMXCASE AMXCTRL AMXTRCASE AMXTRCTRL AMXTECASE AMXTECTRL










%% GET CHI SQUARE P-VALUE FOR EACH GENE-LEVEL VARIANT SUM



[TRCHIP,TRCHIOR]=chisq(AMX.GTRCASEREF,AMX.GTRCASEALT,AMX.GTRCTRLREF,AMX.GTRCTRLALT);
[TECHIP,TECHIOR]=chisq(AMX.GTECASEREF,AMX.GTECASEALT,AMX.GTECTRLREF,AMX.GTECTRLALT);

AMX.TRCHIP  = TRCHIP;
AMX.TRCHIOR = TRCHIOR;

AMX.TECHIP  = TECHIP;
AMX.TECHIOR = TECHIOR;




clearvars -except ADSP GENB LOCI CASE CTRL PHEN...
CTRCASE CTRCTRL CTECASE CTECTRL TRCASEN TRCTRLN TECASEN TECTRLN...
AMX AMXCASE AMXCTRL AMXTRCASE AMXTRCTRL AMXTECASE AMXTECTRL


















%}






%% ########################################################################
%
%
%
%                 QQ-PLOTS   &   MANHATTAN PLOTS
%
%
%
%##########################################################################









%% MAKE SURE MAIN AMX VARIABLES ARE SET


VX = LOCI;
RX = LOCI;



VX.FISHP      = LOCI.TRFISHP;
VX.FISHOR     = LOCI.TRFISHOR;
VX.CASEREF    = LOCI.TRCASEREF;
VX.CASEALT    = LOCI.TRCASEALT;
VX.CTRLREF    = LOCI.TRCTRLREF;
VX.CTRLALT    = LOCI.TRCTRLALT;



RX.FISHP      = LOCI.TEFISHP;
RX.FISHOR     = LOCI.TEFISHOR;
RX.CASEREF    = LOCI.TECASEREF;
RX.CASEALT    = LOCI.TECASEALT;
RX.CTRLREF    = LOCI.TECTRLREF;
RX.CTRLALT    = LOCI.TECTRLALT;


VX.nLogP = -log(VX.FISHP);
RX.nLogP = -log(RX.FISHP);




[~,i] = sort(VX.FISHP);
VX = VX(i,:);

[~,i] = sort(RX.FISHP);
RX = RX(i,:);


clearvars -except ADSP LOCI CASE CTRL PHEN RCASE RCTRL RPHEN VX RX



%% CREATE ADMX TABLES FOR VX & RX

MIS = (strcmp(cellstr(VX.EFFECT),'MIS'));
SYN = (strcmp(cellstr(VX.EFFECT),'SYN'));
SPA = (strcmp(cellstr(VX.EFFECT),'SPA'));
SPD = (strcmp(cellstr(VX.EFFECT),'SPD'));
STG = (strcmp(cellstr(VX.EFFECT),'STG'));
STL = (strcmp(cellstr(VX.EFFECT),'STL'));
LOF = sum([SPA SPD STG STL],2)>0;

LCASYN = VX(SYN,:);
LCAMIS = VX(MIS,:);
LCALOF = VX(LOF,:);


MIS = (strcmp(cellstr(RX.EFFECT),'MIS'));
SYN = (strcmp(cellstr(RX.EFFECT),'SYN'));
SPA = (strcmp(cellstr(RX.EFFECT),'SPA'));
SPD = (strcmp(cellstr(RX.EFFECT),'SPD'));
STG = (strcmp(cellstr(RX.EFFECT),'STG'));
STL = (strcmp(cellstr(RX.EFFECT),'STL'));
LOF = sum([SPA SPD STG STL],2)>0;

LRASYN = RX(SYN,:);
LRAMIS = RX(MIS,:);
LRALOF = RX(LOF,:);


% LOCI-LEVEL TABLES
% LCASYN  Loci CASE-CTRL (ALL)  SYN
% LCAMIS  Loci CASE-CTRL (ALL)  MIS
% LCALOF  Loci CASE-CTRL (ALL)  LOF
% LRASYN  Loci RAND-RAND (ALL)  SYN
% LRAMIS  Loci RAND-RAND (ALL)  MIS
% LRALOF  Loci RAND-RAND (ALL)  LOF


clc; clearvars -except ADSP LOCI CASE CTRL PHEN RCASE RCTRL RPHEN VX RX...
LCASYN LCAMIS LCALOF LRASYN LRAMIS LRALOF



% save('ADSPvsRAND.mat')




%#########################################################################
%% REMOVE APOE ROWS OTHERWISE GRAPHS HAVE HEAVY SKEW
%#########################################################################

doRemove = 1;

if doRemove == 1
LCASYN(contains(cellstr(LCASYN.GENE),'APOE'),:)=[];
LCAMIS(contains(cellstr(LCAMIS.GENE),'APOE'),:)=[];
LCALOF(contains(cellstr(LCALOF.GENE),'APOE'),:)=[];
LRASYN(contains(cellstr(LRASYN.GENE),'APOE'),:)=[];
LRAMIS(contains(cellstr(LRAMIS.GENE),'APOE'),:)=[];
LRALOF(contains(cellstr(LRALOF.GENE),'APOE'),:)=[];


LCASYN(contains(cellstr(LCASYN.GENE),'TOMM40'),:)=[];
LCAMIS(contains(cellstr(LCAMIS.GENE),'TOMM40'),:)=[];
LCALOF(contains(cellstr(LCALOF.GENE),'TOMM40'),:)=[];
LRASYN(contains(cellstr(LRASYN.GENE),'TOMM40'),:)=[];
LRAMIS(contains(cellstr(LRAMIS.GENE),'TOMM40'),:)=[];
LRALOF(contains(cellstr(LRALOF.GENE),'TOMM40'),:)=[];


LCASYN(contains(cellstr(LCASYN.GENE),'MUC'),:)=[];
LCAMIS(contains(cellstr(LCAMIS.GENE),'MUC'),:)=[];
LCALOF(contains(cellstr(LCALOF.GENE),'MUC'),:)=[];
LRASYN(contains(cellstr(LRASYN.GENE),'MUC'),:)=[];
LRAMIS(contains(cellstr(LRAMIS.GENE),'MUC'),:)=[];
LRALOF(contains(cellstr(LRALOF.GENE),'MUC'),:)=[];

% LRASYN(contains(cellstr(LRASYN.GENE),'FLG'),:)=[];
% LRAMIS(contains(cellstr(LRAMIS.GENE),'FLG'),:)=[];

end



clc; clearvars -except ADSP LOCI CASE CTRL PHEN RCASE RCTRL RPHEN VX RX...
LCASYN LCAMIS LCALOF LRASYN LRAMIS LRALOF




%% REMOVE APOE ROWS FROM GENE-LEVEL VARIABLES

if doRemove == 1

% REMOVE TOP-2 OUTLIERS FROM RAND-RAND GENE LEVEL

[i,j] = min(LRASYN.FISHP);    LRASYN(j,:) = [];
[i,j] = min(LRASYN.FISHP);    LRASYN(j,:) = [];
[i,j] = min(LRAMIS.FISHP);    LRAMIS(j,:) = [];
[i,j] = min(LRAMIS.FISHP);    LRAMIS(j,:) = [];

clc; clearvars -except ADSP LOCI CASE CTRL PHEN RCASE RCTRL RPHEN VX RX...
LCASYN LCAMIS LCALOF LRASYN LRAMIS LRALOF

end









%% ########################################################################
%
%
%
%                 QQ-PLOTS   &   MANHATTAN PLOTS
%
%
%
%##########################################################################






























%#############################################################
%%   QQ-PLOTS   LOCI-LEVEL   ALL & ASYMMETRIC
%#############################################################

clc; close all;
figure('Units','normalized','OuterPosition',[.05 .05 .85 .91],'Color','w');
hax1 = axes('Position',[.06 .56 .90 .39],'Color','none');
hax2 = axes('Position',[.06 .06 .90 .39],'Color','none');


% [ALL][SYNONYMOUS][LOCI]

axes(hax1)
qq1 = qqplot( -log10(LRASYN.FISHP) , -log10(LCASYN.FISHP) );
xlabel('-log10 Fishers Exact P RAND:RAND')
ylabel('-log10 Fishers Exact P CASE:CTRL')
title('ALL  |  SYNONYMOUS  |  LOCI')
qq1(3).Color = [.98 .98 .98];
qq1(1).Marker = 'o'; grid on;
qq1(1).MarkerSize = 9;
qq1(1).MarkerFaceColor = [.2 .3 .9];
qq1(1).MarkerEdgeColor = 'none';
% % Plot Polyfit line
% Y = quantile(-log10(LCASYN.FISHP),linspace(.1,.9998,20));
% X = quantile(-log10(LRASYN.FISHP),linspace(.1,.9999,20));
% fitpoly2=fit(X',Y','poly2');
% fX = [0 X max(-log10(LRASYN.FISHP))]';
% fY = fitpoly2(fX);
% line(hax1,fX,fY,'LineStyle','--');

hax1.XLim(1) = 0;
hax1.YLim(1) = 0;

axis tight
hax1.XLim(2) = hax1.YLim(2);

line(hax1,[0 hax1.XLim(2)],[0 hax1.YLim(2)],'LineStyle','-');


% refline(hax1,[0 mean(-log10(LCASYN.FISHP))])


% figure
% loglog(qq1(1).XData, qq1(1).YData)



% [ALL][MISSENSE][LOCI]

axes(hax2)
qq2 = qqplot( -log10(LRAMIS.FISHP) , -log10(LCAMIS.FISHP) );
xlabel('-log10 Fishers Exact P RAND:RAND')
ylabel('-log10 Fishers Exact P CASE:CTRL')
title('ALL  |  MISSENSE  |  LOCI')
qq2(3).Color = [.98 .98 .98];
qq2(1).Marker = 'o'; grid on;
qq2(1).MarkerSize = 9;
qq2(1).MarkerFaceColor = [.9 .3 .3];
qq2(1).MarkerEdgeColor = 'none';
% % Plot Polyfit line
% Y = quantile(-log10(LCAMIS.FISHP),linspace(.1,.9998,20));
% X = quantile(-log10(LRAMIS.FISHP),linspace(.1,.9999,20));
% fitpoly2=fit(X',Y','poly2');
% fX = [0 X max(-log10(LRAMIS.FISHP))]';
% fY = fitpoly2(fX);
% line(hax2,fX,fY,'LineStyle','--');

hax2.XLim(1) = 0;
hax2.YLim(1) = 0;

axis tight
hax2.XLim(2) = hax2.YLim(2);

line(hax2,[0 hax2.XLim(2)],[0 hax2.YLim(2)],'LineStyle','-');









% ADD TEST LABLES OF GENE NAMES TO QQ PLOTS
[~,j] = sort(LCASYN.FISHP);  LCASYN_GENE  =  LCASYN.GENE(j,:);
[~,j] = sort(LCAMIS.FISHP);  LCAMIS_GENE  =  LCAMIS.GENE(j,:);

text(hax1, fliplr(qq1(1).XData(end-5:end)) , fliplr(qq1(1).YData(end-5:end)) ,...
     LCASYN_GENE(1:6,:) , 'HorizontalAlignment','right');
text(hax2, fliplr(qq2(1).XData(end-5:end)) , fliplr(qq2(1).YData(end-5:end)) ,...
     LCAMIS_GENE(1:6,:) , 'HorizontalAlignment','right');


clc; clearvars -except ADSP LOCI CASE CTRL PHEN RCASE RCTRL RPHEN VX RX...
LCASYN LCAMIS LCALOF LRASYN LRAMIS LRALOF


% saveas(gcf,'APOE_QQ_LOCI_ALL_ASYM.png')












%#############################################################
%%   QQ-PLOTS   LOCI-LEVEL   ALL & ASYMMETRIC
%#############################################################

clc; close all;
figure('Units','normalized','OuterPosition',[.05 .05 .85 .91],'Color','w');
hax1 = axes('Position',[.06 .06 .90 .90],'Color','none');

% [ALL][SYNONYMOUS][LOCI]

axes(hax1)
qq1 = qqplot( -log10(LRASYN.FISHP) , -log10(LCASYN.FISHP) );
xlabel('-log10 Fishers Exact P RAND:RAND')
ylabel('-log10 Fishers Exact P CASE:CTRL')
title('ALL  |  SYNONYMOUS  |  LOCI')
qq1(3).Color = [.98 .98 .98];
qq1(1).Marker = 'o'; grid on;
qq1(1).MarkerSize = 9;
qq1(1).MarkerFaceColor = [.2 .3 .9];
qq1(1).MarkerEdgeColor = 'none';
% % Plot Polyfit line
% Y = quantile(-log10(LCASYN.FISHP),linspace(.1,.9998,20));
% X = quantile(-log10(LRASYN.FISHP),linspace(.1,.9999,20));
% fitpoly2=fit(X',Y','poly2');
% fX = [0 X max(-log10(LRASYN.FISHP))]';
% fY = fitpoly2(fX);
% line(hax1,fX,fY,'LineStyle','--');

hax1.XLim(1) = 0;
hax1.YLim(1) = 0;

axis tight
hax1.XLim(2) = hax1.YLim(2);

line(hax1,[0 hax1.XLim(2)],[0 hax1.YLim(2)],'LineStyle','-');

% ADD TEST LABLES OF GENE NAMES TO QQ PLOTS

[~,j] = sort(LCASYN.FISHP);  LCASYN_GENE  =  LCASYN.GENE(j,:);

text(hax1, fliplr(qq1(1).XData(end-5:end)) , fliplr(qq1(1).YData(end-5:end)) ,...
     LCASYN_GENE(1:6,:) , 'HorizontalAlignment','right');













% [ALL][MISSENSE][LOCI]
figure('Units','normalized','OuterPosition',[.05 .05 .85 .91],'Color','w');
hax2 = axes('Position',[.06 .06 .90 .90],'Color','none');




axes(hax2)
qq2 = qqplot( -log10(LRAMIS.FISHP) , -log10(LCAMIS.FISHP) );
xlabel('-log10 Fishers Exact P RAND:RAND')
ylabel('-log10 Fishers Exact P CASE:CTRL')
title('ALL  |  MISSENSE  |  LOCI')
qq2(3).Color = [.98 .98 .98];
qq2(1).Marker = 'o'; grid on;
qq2(1).MarkerSize = 9;
qq2(1).MarkerFaceColor = [.9 .3 .3];
qq2(1).MarkerEdgeColor = 'none';
hax2.XLim(1) = 0;
hax2.YLim(1) = 0;
axis tight
hax2.XLim(2) = hax2.YLim(2);

line(hax2,[0 hax2.XLim(2)],[0 hax2.YLim(2)],'LineStyle','-');



% ADD TEST LABLES OF GENE NAMES TO QQ PLOTS
[~,j] = sort(LCAMIS.FISHP);  LCAMIS_GENE  =  LCAMIS.GENE(j,:);
text(hax2, fliplr(qq2(1).XData(end-5:end)) , fliplr(qq2(1).YData(end-5:end)) ,...
     LCAMIS_GENE(1:6,:) , 'HorizontalAlignment','right');


clc; clearvars -except ADSP LOCI CASE CTRL PHEN RCASE RCTRL RPHEN VX RX...
LCASYN LCAMIS LCALOF LRASYN LRAMIS LRALOF


% saveas(gcf,'APOE_QQ_LOCI_ALL_ASYM.png')

