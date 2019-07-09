%==========================================================================
%% STEP-1: LOAD THE DATASET
%==========================================================================
close all; clear; clc; rng('shuffle');
P.home = fileparts(which('GENOS.m')); cd(P.home);
P.P1 = [P.home filesep 'genos_functions'];
P.P3 = [P.P1 filesep 'genos_main_functions'];
P.P4 = [P.home filesep 'genos_other'];
addpath(join(string(struct2cell(P)),pathsep,1))
cd(P.home); P.f = filesep;


%==========================================================================
%% (GENOX) PERTURB TEST OPTIONS & FOLDERS 
%==========================================================================
clc; clearvars -except P


% SET DEFAULT FOLDERS 
%------------------------------------------------------
%{
% P.genox.DIRroot = [P.home P.f 'genos_data' P.f 'GENOX'];
% P.genox.DIRmat  = [P.genox.DIRroot P.f 'PLO_MAT'];
% P.genox.DIRimg  = [P.genox.DIRroot P.f 'PLO_IMG'];
% P.genox.DIRstatsmat = [P.genox.DIRroot P.f 'PLO_STATS_MAT'];
% P.genox.DIRstatsimg = [P.genox.DIRroot P.f 'PLO_STATS_IMG'];
% 
% 
% P.genox.DIRroot = [P.home P.f 'genos_data' P.f 'GENOX'];
% P.genox.DIRmat  = [P.genox.DIRroot P.f 'ORLO_MAT'];
% P.genox.DIRimg  = [P.genox.DIRroot P.f 'ORLO_IMG'];
% P.genox.DIRstatsmat = [P.genox.DIRroot P.f 'ORLO_STATS_MAT'];
% P.genox.DIRstatsimg = [P.genox.DIRroot P.f 'ORLO_STATS_IMG'];
% 
% 
% P.genox.DIRroot = [P.home P.f 'genos_data' P.f 'GENOX'];
% P.genox.DIRmat  = [P.genox.DIRroot P.f 'ORHI_MAT'];
% P.genox.DIRimg  = [P.genox.DIRroot P.f 'ORHI_IMG'];
% P.genox.DIRstatsmat = [P.genox.DIRroot P.f 'ORHI_STATS_MAT'];
% P.genox.DIRstatsimg = [P.genox.DIRroot P.f 'ORHI_STATS_IMG'];
% 
% 
% P.genox.DIRroot = [P.home P.f 'genos_data' P.f 'GENOX'];
% P.genox.DIRmat  = [P.genox.DIRroot P.f 'SYN_MAT'];
% P.genox.DIRimg  = [P.genox.DIRroot P.f 'SYN_IMG'];
% P.genox.DIRstatsmat = [P.genox.DIRroot P.f 'SYN_STATS_MAT'];
% P.genox.DIRstatsimg = [P.genox.DIRroot P.f 'SYN_STATS_IMG'];
%}


P.genox.DIRroot = 'F:\GENOSDATA\GENOX\PLO';
P.genox.DIRmat  = [P.genox.DIRroot P.f 'PLO_MAT'];
P.genox.DIRimg  = [P.genox.DIRroot P.f 'PLO_MAT'];
P.genox.DIRstatsmat = [P.genox.DIRroot P.f 'PLO_STATS_MAT'];
P.genox.DIRstatsimg = [P.genox.DIRroot P.f 'PLO_STATS_IMG'];




% IMPORT MAT FILES
%------------------------------------------------------
P.PERT.w = what(P.genox.DIRmat);
P.PERT.finfo = dir(P.PERT.w.path);
P.PERT.finames = {P.PERT.finfo.name};
c=~cellfun(@isempty,regexp(P.PERT.finames,'((\S)+(\.mat+))'));
P.PERT.finames = string(P.PERT.finames(c)');
P.PERT.folder = P.PERT.finfo.folder;
P.PERT.fipaths = fullfile(P.PERT.folder,P.PERT.finames);
disp(P.PERT.fipaths); disp(P.PERT.finames);

PERT = load(P.PERT.fipaths{1});

clearvars -except P PERT


%==========================================================================
%% GET PERTURB STATS
%==========================================================================
clc; clearvars -except P PERT


for ij = 1:numel(P.PERT.fipaths)
clc; clearvars -except P ADSP INFO PER ij
disp(ij); disp(P.PERT.finames{ij})


% LOAD SNP FILE
%-----------------------------------
PERT = load(P.PERT.fipaths{ij});


% DETERMINE GENE, CHRPOS, AND NLOOPS
%-----------------------------------
PER(1).GENE   = PERT.INFO.TOPLOCI{1,1}.GENE(1);
PER(1).CHRPOS = PERT.INFO.TOPLOCI{1,1}.CHRPOS(1);
NLoops = size(PERT.LOOPDATA.AREA_HO,3);




% NAT/NAT PERFORMANCE
%-----------------------------------
PER(1).area = mean(PERT.LOOPDATA.AREA_HO(:,:,NLoops,1) ,3);
PER(1).caseMu = midgravity(PER(1).area(:,1), PER(1).area(:,2), PER(1).area(:,4));
PER(1).ctrlMu = midgravity(PER(1).area(:,1), PER(1).area(:,3), PER(1).area(:,5));
PER(1).cacoZd = 0;


% A = mean(PERT.LOOPDATA.AREA_HO(:,:,1:10,1) ,3);
% close all
% area(A(:,1),A(:,2:end))


% REF/REF PERFORMANCE
%-----------------------------------
PER(2).area = mean(PERT.LOOPDATA.AREA_HO(:,:,NLoops,2) ,3);
PER(2).caseMu = midgravity(PER(2).area(:,1), PER(2).area(:,2), PER(2).area(:,4));
PER(2).ctrlMu = midgravity(PER(2).area(:,1), PER(2).area(:,3), PER(2).area(:,5));
PER(2).cacoZd = (abs(PER(2).caseMu - PER(1).caseMu)+...
                 abs(PER(2).ctrlMu - PER(1).ctrlMu))/2;


% REF/ALT PERFORMANCE
%-----------------------------------
PER(3).area = mean(PERT.LOOPDATA.AREA_HO(:,:,NLoops,3) ,3);
PER(3).caseMu = midgravity(PER(3).area(:,1), PER(3).area(:,2), PER(3).area(:,4));
PER(3).ctrlMu = midgravity(PER(3).area(:,1), PER(3).area(:,3), PER(3).area(:,5));
PER(3).cacoZd = (abs(PER(3).caseMu - PER(1).caseMu)+...
                 abs(PER(3).ctrlMu - PER(1).ctrlMu))/2;





% ALT/ALT PERFORMANCE
%-----------------------------------
PER(4).area = mean(PERT.LOOPDATA.AREA_HO(:,:,NLoops,4) ,3);
PER(4).caseMu = midgravity(PER(4).area(:,1), PER(4).area(:,2), PER(4).area(:,4));
PER(4).ctrlMu = midgravity(PER(4).area(:,1), PER(4).area(:,3), PER(4).area(:,5));
PER(4).cacoZd = (abs(PER(4).caseMu - PER(1).caseMu)+...
                 abs(PER(4).ctrlMu - PER(1).ctrlMu))/2;



%------------------------------------------%
P.SAVEPATH = [P.genox.DIRstatsmat P.f ...
              PER(1).GENE{1} '_' num2str(PER(1).CHRPOS) '.mat'];
% save(P.SAVEPATH,'P','PERT','PER');
save(P.SAVEPATH,'PER');
disp('File saved...'); disp(P.SAVEPATH)
%------------------------------------------%

end


clc; clearvars -except P PERT PER
%==========================================================================
%% QUANTIFY SHIFTS AND EXPORT TO EXCEL
%==========================================================================
clc; clearvars -except P ADSP INFO PERT PER

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
STAT(1).GENE   = GENOS.PER.GENE;
STAT(1).CHRPOS = GENOS.PER.CHRPOS;
STAT(1).AREA   = GENOS.PER.area;
STAT(1).CASEMU = GENOS.PER.caseMu;
STAT(1).CTRLMU = GENOS.PER.ctrlMu;
STAT(1).CACOMU = GENOS.PER.cacoZd;


GENOS = load(P.PERT.fipaths{1});
%---
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


STATISTICS = STATS(:,[1 2 4 5 6 7]);
STATISTICS.CASEHIST = zeros(size(STATISTICS,1),50);
STATISTICS.CTRLHIST = zeros(size(STATISTICS,1),50);


for k = 1:size(STATISTICS,1)

    STATISTICS.CASEHIST(k,:) = STATS.area{k,1}(:,4)' + STATS.area{k,1}(:,2)';

    STATISTICS.CTRLHIST(k,:) = STATS.area{k,1}(:,5)' + STATS.area{k,1}(:,3)';

end

writetable(STATISTICS, [P.genox.DIRroot P.f 'SNPS.xlsx'],'Sheet','RAW')


clc; clearvars -except P PERT PER STATS STATISTICS


%==========================================================================
%% PLOT CONFUSION PERFORMANCE HISTOGRAMS
%==========================================================================
clc; clearvars -except P PERT PER STATS STATISTICS

% IMPORT PEROX FILES
%-----------------------------------
P.PERX.w = what(P.genox.DIRstatsmat);
P.PERX.finfo = dir(P.PERX.w.path);
P.PERX.finames = {P.PERX.finfo.name};
c=~cellfun(@isempty,regexp(P.PERX.finames,'((\S)+(\.mat+))'));
P.PERX.finames = string(P.PERX.finames(c)');
P.PERX.folder = P.PERX.finfo.folder;
P.PERX.fipaths = fullfile(P.PERX.folder,P.PERX.finames);
disp(P.PERX.fipaths); disp(P.PERX.finames);


for ij = 1:numel(P.PERX.fipaths)
clc; clearvars -except P PERT PER STATS STATISTICS ij
disp(ij);



load(P.PERX.fipaths{ij});


f = char(P.PERX.finames(ij));
gnam = f(1:end-4);



close all
fh01 = figure('Units','normalized','OuterPosition',[.01 .06 .97 .88],'Color','w');

ax10 = axes('Position',[.05 .59 .31 .37],'Color','none');
ax11 = axes('Position',[.42 .59 .06 .37],'Color','none');

ax20 = axes('Position',[.54 .59 .31 .37],'Color','none');
ax21 = axes('Position',[.91 .59 .06 .37],'Color','none');

ax30 = axes('Position',[.05 .09 .31 .37],'Color','none');
ax31 = axes('Position',[.42 .09 .06 .37],'Color','none');

ax40 = axes('Position',[.54 .09 .31 .37],'Color','none');
ax41 = axes('Position',[.91 .09 .06 .37],'Color','none');




%--------------------------------------------------------------------------
%---   NAT/NAT | PLOT CONFUSION PERFORMANCE HISTOGRAM
%--------------------------------------------------------------------------

%---STACKED HISTOGRAM ---------------------
    axes(ax10)
ph0=area(PER(1).area(:,1),PER(1).area(:,2:5)); hold on;
    ax10.ColorOrderIndex = 1;
ph1=bar( PER(1).area(:,1), PER(1).area(:,2:5) ,'stacked');
    ylabel('Count'); 
    ax10.FontSize=16; 
    lh=legend([ph1],{'Case Miss','Ctrl Miss','Case Hit','Ctrl Hit'},...
    'Location','best');
    title([gnam ' | APOE33 | SNPxx'],'interpreter','none')
    hold on
    %pause(.1); 


%---MID LINES ON HISTOGRAM ---------------------
    Y = ax10.YLim(1):50:ax10.YLim(2);
    Xa = repmat(PER(1).caseMu,1,numel(Y));
    Xb = repmat(PER(1).ctrlMu,1,numel(Y));
ph2=plot(Xa,Y,'.-',Xb,Y,'.-','MarkerSize',30,'LineWidth',6);
    ph2(1).Color = [.5 .1 .1];
    ph2(2).Color = [.1 .2 .6];
    lh.String(5:6) = {'Case Mu','Ctrl Mu'};
    hold on
    %pause(.2); 


%---BAR GRAPH ---------------------
    axes(ax11)
ph3=bar([.6 2 3.4],[ PER(1).caseMu , PER(1).ctrlMu ,  PER(1).cacoZd],...
    .5,'FaceColor',[.31 .31 .31]); 
    %---BAR GRAPH FORMATTING ------
    ylabel('deviation from zero')
    ax11.YGrid = 'on';
    ax11.YLim = [-.5 .5]; 
    ax11.XTickLabels = {'CASE','CTRL','AbsMd'};
    ax11.XTickLabelRotation = 45;
    ax11.FontSize=16;
    % ax11.YTick = [0 25 50 75 100]; 
    %title('Performance Summary')





%--------------------------------------------------------------------------
%---   REF/REF | PLOT CONFUSION PERFORMANCE HISTOGRAM
%--------------------------------------------------------------------------


%---STACKED HISTOGRAM ---------------------
    axes(ax20)
ph1=area( PER(2).area(:,1), PER(2).area(:,2:5) ); hold on;
    ax20.ColorOrderIndex = 1;
ph1=bar( PER(2).area(:,1), PER(2).area(:,2:5) ,'stacked');
    ylabel('Count'); 
    ax20.FontSize=16; 
    lh=legend([ph1],{'Case Miss','Ctrl Miss','Case Hit','Ctrl Hit'},...
    'Location','best');
    title([gnam ' | APOE33 | SNP33'],'interpreter','none')
    hold on
    %pause(.1); 


%---MID LINES ON HISTOGRAM ---------------------
    Y = ax20.YLim(1):50:ax20.YLim(2);
    Xa = repmat(PER(2).caseMu,1,numel(Y));
    Xb = repmat(PER(2).ctrlMu,1,numel(Y));
ph2=plot(Xa,Y,'.-',Xb,Y,'.-','MarkerSize',30,'LineWidth',6);
    ph2(1).Color = [.5 .1 .1];
    ph2(2).Color = [.1 .2 .6];
    lh.String(5:6) = {'Case Mu','Ctrl Mu'};
    hold on
    %pause(.2); 



%---BAR GRAPH ---------------------
    axes(ax21)
ph3=bar([.6 2 3.4], [ PER(2).caseMu , PER(2).ctrlMu ,  PER(2).cacoZd],...
    .5,'FaceColor',[.31 .31 .31]); 
    %---BAR GRAPH FORMATTING ------
    grid on; ylabel('zero deviance weight')
    ax21.YLim = [-.5 .5]; 
    ax21.XTickLabels = {'CASE','CTRL','AbsMd'};
    ax21.XTickLabelRotation = 45;
    ax21.FontSize=16;
    % ax21.YTick = [0 25 50 75 100]; 
    %title('Performance Summary')







%--------------------------------------------------------------------------
%---   REF/ALT | PLOT CONFUSION PERFORMANCE HISTOGRAM
%--------------------------------------------------------------------------

%---STACKED HISTOGRAM ---------------------
    axes(ax30)
ph1=area(PER(3).area(:,1),PER(3).area(:,2:5)); hold on;
    ax30.ColorOrderIndex = 1;
ph1=bar( PER(3).area(:,1), PER(3).area(:,2:5) ,'stacked');
    ylabel('Count'); 
    ax30.FontSize=16; 
    lh=legend([ph1],{'Case Miss','Ctrl Miss','Case Hit','Ctrl Hit'},...
    'Location','best'); 
    title([gnam ' | APOE33 | SNP34'],'interpreter','none')
    hold on
    %pause(.1); 
    

%---MID LINES ON HISTOGRAM ---------------------
    Y = ax30.YLim(1):50:ax30.YLim(2);
    Xa = repmat(PER(3).caseMu,1,numel(Y));
    Xb = repmat(PER(3).ctrlMu,1,numel(Y));
ph2=plot(Xa,Y,'.-',Xb,Y,'.-','MarkerSize',30,'LineWidth',6);
    ph2(1).Color = [.5 .1 .1];
    ph2(2).Color = [.1 .2 .6];
    lh.String(5:6) = {'Case Mu','Ctrl Mu'};
    hold on
    %pause(.2); 


%---BAR GRAPH ---------------------
    axes(ax31)
ph3=bar([.6 2 3.4], [ PER(3).caseMu , PER(3).ctrlMu ,  PER(3).cacoZd],...
    .5,'FaceColor',[.31 .31 .31]); 
    %---BAR GRAPH FORMATTING ------
    grid on; ylabel('zero deviance weight')
    ax31.YLim = [-.5 .5]; 
    ax31.XTickLabels = {'CASE','CTRL','AbsMd'};
    ax31.XTickLabelRotation = 45;
    ax31.FontSize=16;
    % ax31.YTick = [0 25 50 75 100]; 
    %title('Performance Summary')





%--------------------------------------------------------------------------
%---   ALT/ALT | PLOT CONFUSION PERFORMANCE HISTOGRAM
%--------------------------------------------------------------------------

%---STACKED HISTOGRAM ---------------------
    axes(ax40)
ph1=area(PER(4).area(:,1),PER(4).area(:,2:5)); hold on;
    ax40.ColorOrderIndex = 1;
ph1=bar( PER(4).area(:,1), PER(4).area(:,2:5) ,'stacked');
    ylabel('Count'); 
    ax40.FontSize=16; 
    lh=legend([ph1],{'Case Miss','Ctrl Miss','Case Hit','Ctrl Hit'},...
    'Location','best');
    title([gnam ' | APOE33 | SNP44'],'interpreter','none');
    hold on
    %pause(.1); 


%---MID LINES ON HISTOGRAM ---------------------
    Y = ax40.YLim(1):50:ax40.YLim(2);
    Xa = repmat(PER(4).caseMu,1,numel(Y));
    Xb = repmat(PER(4).ctrlMu,1,numel(Y));
ph2=plot(Xa,Y,'.-',Xb,Y,'.-','MarkerSize',30,'LineWidth',6);
    ph2(1).Color = [.5 .1 .1];
    ph2(2).Color = [.1 .2 .6];
    lh.String(5:6) = {'Case Mu','Ctrl Mu'};
    hold on
    %pause(.2); 


%---BAR GRAPH ---------------------
    axes(ax41)
ph3=bar([.6 2 3.4], [ PER(4).caseMu , PER(4).ctrlMu ,  PER(4).cacoZd],...
    .5,'FaceColor',[.31 .31 .31]); 
    %---BAR GRAPH FORMATTING ------
    grid on; ylabel('zero deviance weight')
    ax41.YLim = [-.5 .5]; 
    ax41.XTickLabels = {'CASE','CTRL','AbsMd'};
    ax41.XTickLabelRotation = 45;
    ax41.FontSize=16;
    % ax41.YTick = [0 25 50 75 100]; 
    %title('Performance Summary')
%----------------------------------------------------------------------
pause(1)
set(gcf, 'PaperPositionMode', 'auto');
finam = char(P.PERX.finames(ij));
saveas(gcf, [P.genox.DIRstatsimg P.f finam(1:end-4) '.png']);
pause(1)
%---------------------------------------------------------------------- 


end




%==========================================================================
%% xxx
%==========================================================================
