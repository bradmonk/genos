function [] = genosmanhattan(LOCI)


% keyboard

%% STORE VARIABLES FOR NEURAL NETWORK TRAINING AS 'VLOCI'
%clc; clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL

% SET MAIN FISHP TO TRAINING GROUP FISHP
VLOCI = LOCI;
VLOCI.FISHP      = LOCI.TRFISHP;
VLOCI.FISHOR     = LOCI.TRFISHOR;

VLOCI.CASEREF    = LOCI.TRCASEREF;
VLOCI.CASEALT    = LOCI.TRCASEALT;
VLOCI.CTRLREF    = LOCI.TRCTRLREF;
VLOCI.CTRLALT    = LOCI.TRCTRLALT;


% SORT VARIANTS BY CHR:POS
% [~,i]  = sort(VLOCI.TRFISHP);
[~,i]  = sort(VLOCI.CHRPOS);
VLOCI  = VLOCI(i,:);




%----------------
clc; close all;
fh1 = figure('Units','normalized','OuterPosition',[.01 .05 .9 .9],...
             'Color','w','MenuBar','none');
ax1 = axes('Position',[.05 .55 .9 .4],'Color','none',...
           'XColor','none','XLabel',[]); hold on;
ax2 = axes('Position',[.05 .05 .9 .4],'Color','none',...
           'XColor','none','XLabel',[]); hold on;


%----------------
    AH = -log(VLOCI.TRFISHP);
    PAH = AH > 20;
    AH(PAH) = 20-rand*2;


axes(ax1)
ph1 = scatter(1:size(VLOCI,1),  AH  ,300, VLOCI.CHR , '.');

    ylabel('\fontsize{22}-log(P)')
    ax1.FontSize = 16;
    bonfapha = 1e-04;
    line(ax1,[1 size(VLOCI,1)],-log([bonfapha bonfapha]),...
        'Color','red','LineStyle','--','LineWidth',2)
    axis tight
    ax1.YLim = [0 25];
%----------------


%----------------
    MH = -log10(VLOCI.TEFISHP);
    PMH = MH > 20;
    MH(PMH) = 20-rand*2;

axes(ax2)
ph2 = scatter(1:size(VLOCI,1),  MH  ,300, VLOCI.CHR , '.');

    ylabel('\fontsize{22}-log(P)')
    ax2.FontSize = 16;

    bonfapha = 1e-04; %.05/size(VLOCI,1); %
    line(ax2,[1 size(VLOCI,1)],-log([bonfapha bonfapha]),...
        'Color','red','LineStyle','--','LineWidth',2)
    axis tight
    ax2.YLim = [0 25];
%----------------




% clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
% VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL



%---- ADD GENE NAME TEXT LABELS TO GRAPH
nL = 90;
TRnlpCUTOFF = -log(5e-05);
TEnlpCUTOFF = -log(5e-05);

VLOCI.ORD = (1:size(VLOCI,1))';


% SORT TABLE BY P-VALUE
[i,j] = sort(VLOCI.TRFISHP);
GENOGENE = VLOCI.GENE(j,:);
GENOORD  = VLOCI.ORD(j,:);
GENOP    = -log(VLOCI.TRFISHP(j,:));
I = GENOP>25;
K = (exp(exp(log(log(GENOP(I)))))./25-1)./2;
GENOP(I) = 24 + K;


[i,j] = sort(VLOCI.TEFISHP);
GENODGENE = VLOCI.GENE(j,:);
GENODORD  = VLOCI.ORD(j,:);
GENODP    = -log(VLOCI.TEFISHP(j,:));
I = GENODP>25;
K = (exp(exp(log(log(GENODP(I)))))./25-1)./2;
GENODP(I) = 24 + K;


TRi = GENOP  > TRnlpCUTOFF;
TEi = GENODP > TEnlpCUTOFF;


Gi = GENOGENE(TRi,:);
Xi = GENOORD(TRi,:);
Yi = GENOP(TRi,:);

Gj = GENODGENE(TEi,:);
Xj = GENODORD(TEi,:);
Yj = GENODP(TEi,:);

axes(ax1)
text(Xi , Yi , char(Gi) , 'HorizontalAlignment','right');


axes(ax2)
text(Xj , Yj , char(Gj) , 'HorizontalAlignment','right');



% clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
% VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL...
% fh1 ax1 ax2 Gi Xi Yi Gj Xj Yj



% %---------- EXPORT IMAGE TO DESKTOP FOLDER -----------------------------
% dt=char(datetime(datetime,'Format','yyyy-MM-dd-HH-mm-ss'));
% print(['/Users/bradleymonk/Desktop/MANHATTAN ' dt '.png'],'-dpng');
% %-----------------------------------------------------------------------

% clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
% VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL








%% DETERMINE HAPLOTYPES WITH LINKAGE DISEQUILIBRIUM
%{
clc; clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL



% SORT VARIANTS BY FISHERS P-VALUE
% [~,i]  = sort(VLOCI.TRFISHP);
i = VLOCI.TRFISHP < .005;
VLOCI  = VLOCI(i,:);
VCASE  = VCASE(i);
VCTRL  = VCTRL(i);
VUSNP  = VUSNP(i);



% COMPUTE WHETHER P-VALUES UPSTREAM AND DOWNSTREAM ARE DECREASING
PPOS = VLOCI.TRFISHP(2:end-1);
PPUP = VLOCI.TRFISHP(3:end);
PPDN = VLOCI.TRFISHP(1:end-2);


PREVP = PPOS < PPUP;
NEXTP = PPOS < PPDN;
VLOCI.PREVP = [0;PREVP;0];
VLOCI.NEXTP = [0;NEXTP;0];
VLOCI.MINIMA = VLOCI.PREVP & VLOCI.NEXTP;

VLOCI.MINIMA(1) = VLOCI.TRFISHP(1) < VLOCI.TRFISHP(2);
VLOCI.MINIMA(end) = VLOCI.TRFISHP(end) < VLOCI.TRFISHP(end-1);


LOLOC = VLOCI(VLOCI.MINIMA==1,:);
%}




%% COMPUTE PROXIMITY TO UPSTREAM AND DOWNSTREAM VARIANT
%{

% LPOS = LOC.CHRPOS(2:end-1);
% LPUP = LOC.CHRPOS(3:end);
% LPDN = LOC.CHRPOS(1:end-2);


APOE_TOMM40_DIST = 412079 - 395714;

% SORT BY CHRPOS
[~,i]  = sort(VLOCI.CHRPOS);
VLOCI  = VLOCI(i,:);
VCASE  = VCASE(i);
VCTRL  = VCTRL(i);
VUSNP  = VUSNP(i);



POS = VLOCI.POS(2:end-1);
NEXTPOS = VLOCI.POS(3:end);
PREVPOS = VLOCI.POS(1:end-2);

DISTPREV = POS - PREVPOS;
DISTNEXT = NEXTPOS - POS;

VLOCI.DISTPREV       = [0;DISTPREV;0];
VLOCI.DISTNEXT       = [0;DISTNEXT;0];
VLOCI.DISTNEXT(1)    = VLOCI.POS(2)-VLOCI.POS(1);

VLOCI.DISTPREV(end)  = VLOCI.POS(end)-VLOCI.POS(end-1);
VLOCI.DISTNEXT(1)    = VLOCI.POS(2)-VLOCI.POS(1);






OK1 = (VLOCI.DISTPREV < APOE_TOMM40_DIST);
OK2 = (VLOCI.DISTNEXT < APOE_TOMM40_DIST);
BAD = BAD1 | BAD2;

VLOCI(BAD,:) = [];
VCASE(i) = [];
VCTRL(i) = [];
VUSNP(i) = [];

% SORT BY CHRPOS
[~,i]  = sort(VLOCI.CHRPOS);
VLOCI  = VLOCI(i,:);
VCASE  = VCASE(i);
VCTRL  = VCTRL(i);
VUSNP  = VUSNP(i);



BAD = (VLOCI.DISTPREV < APOE_TOMM40_DIST) & (VLOCI.MINIMA == 0);
VLOCI(BAD,:) = [];
VCASE(i) = [];
VCTRL(i) = [];
VUSNP(i) = [];

% SORT BY CHRPOS
[~,i]  = sort(VLOCI.CHRPOS);
VLOCI  = VLOCI(i,:);
VCASE  = VCASE(i);
VCTRL  = VCTRL(i);
VUSNP  = VUSNP(i);






close all
% subplot(2,1,1); histogram(-log(1./double(LOC.DISTNEXT)))
% line([-log(1/600) -log(1/600)],[0 3000])


figure
subplot(2,1,1); histogram(log(double(VLOCI.DISTNEXT)))
% line([log(600) log(600)],[0 3000])
subplot(2,1,2); histogram(double(VLOCI.DISTNEXT(VLOCI.DISTNEXT<1e5)))






figure
subplot(2,1,1); histogram(log(double(LOX.DISTNEXT)))
% line([log(600) log(600)],[0 3000])
subplot(2,1,2); histogram(double(LOX.DISTNEXT(LOX.DISTNEXT<1e5)))




% SORT VARIANTS BY FISHERS P-VALUE
[~,i]  = sort(VLOCI.TRFISHP);
VLOCI  = VLOCI(i,:);
VCASE  = VCASE(i);
VCTRL  = VCTRL(i);
VUSNP  = VUSNP(i);


[~,i]  = sort(LOX.CHRPOS);
LOX    = LOX(i,:);
VCASE  = VCASE(i);
VCTRL  = VCTRL(i);
VUSNP  = VUSNP(i);



%}




%% ONLY KEEP ONE OR TWO VARIANTS PER GENE (THE ONES WITH LOWEST P-VALS)
%{

% SORT VARIANTS BY FISHERS P-VALUE
[~,i]  = sort(VLOCI.TRFISHP);
VLOCI  = VLOCI(i,:);
VCASE  = VCASE(i);
VCTRL  = VCTRL(i);
VUSNP  = VUSNP(i);



% GET FIRST TOP RATED VARIANT FOR EACH GENE
[~,i,~] = unique(string(VLOCI.GENE),'stable');

VLO1  = VLOCI(i,:);
VCA1  = VCASE(i);
VCO1  = VCTRL(i);
VUS1  = VUSNP(i);

VLOCI(i,:) = [];
VCASE(i)   = [];
VCTRL(i)   = [];
VUSNP(i)   = [];


% GET SECOND TOP RATED VARIANT FOR EACH GENE
[~,i,~] = unique(string(VLOCI.GENE),'stable');

VLO2  = VLOCI(i,:);
VCA2  = VCASE(i);
VCO2  = VCTRL(i);
VUS2  = VUSNP(i);

VLOCI = [VLO1 ; VLO2];
VCASE = [VCA1 ; VCA2];
VCTRL = [VCO1 ; VCO2];
VUSNP = [VUS1 ; VUS2];



% SORT VARIANTS BY FISHERS P-VALUE
[~,i]  = sort(VLOCI.TRFISHP);
VLOCI  = VLOCI(i,:);
VCASE  = VCASE(i);
VCTRL  = VCTRL(i);
VUSNP  = VUSNP(i);




clc; disp(VLOCI(1:10,:))
clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
VLOCI VCASE VCTRL VUSNP VTRCASE VTRCTRL VTECASE VTECTRL
%}

end