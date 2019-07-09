function [] = nnperfplots(YL,MX,ACT,NETDAT)
%==========================================================================
%% PLOT DATA
%=====================================================================
% clc; clearvars -except YL MX net VIN ACT

% keyboard

%%


OUT = ACT - 0;
nCACO = numel(YL);
nCASE = sum(YL);
nCTRL = sum(~YL);



% PLOT 1
%-----------------------------------
ECDFcenters = (.01:.01:.99);
[ECDFf,ECDFx] = ecdf(OUT);

KERNELedges = (0:.01:1);
pdCACO  = fitdist(OUT','Kernel');
pdCTRL  = fitdist(OUT(1,YL(1,:)~=1)','Kernel');
pdCASE  = fitdist(OUT(1,YL(1,:)==1)','Kernel');
pdCACOy = pdf( pdCACO, KERNELedges );
pdCTRLy = pdf( pdCTRL, KERNELedges );
pdCASEy = pdf( pdCASE, KERNELedges );





% PLOT 2
%-----------------------------------
P2_edges = (-.5:.02:.5);
xCASE = OUT(1,YL(1,:)==1) - .5;
xCTRL = OUT(1,YL(1,:)~=1) - .5;




% PLOT 3
%-----------------------------------
edges = (-.5:.02:.5);
bCASE = histcounts(xCASE,edges);
bCTRL = histcounts(xCTRL,edges);
BARX = [(-.5:.02:.49); (-.49:.02:.5)]';
BARY = [bCTRL; bCASE]';



% PLOT 4
%-----------------------------------
% [KSf,KSi,KSb] = ksdensity(OUT,'Support','positive','Function','cdf');
xCACO = [xCASE(1:min([nCASE nCTRL])); xCTRL(1:min([nCASE nCTRL]))]';
gridx1 = -.5:.01:.5;
gridx2 = -.5:.01:.5;
[x1,x2] = meshgrid(gridx1, gridx2);
x1 = x1(:);
x2 = x2(:);
xi = [x1 x2];





%==========================================================================
close all;
fh01 = figure('Units','normalized','OuterPosition',[.01 .04 .98 .93],...
              'Color','w','MenuBar','none');

ax07 = axes('Position',[.03 .78 .42 .20],'Color','none');
ax08 = axes('Position',[.50 .78 .48 .20],'Color','none');

ax01 = axes('Position',[.03 .42 .28 .31],'Color','none');
ax02 = axes('Position',[.36 .42 .28 .31],'Color','none');
ax03 = axes('Position',[.68 .42 .28 .31],'Color','none');

ax04 = axes('Position',[.03 .05 .30 .33],'Color','none');
ax05 = axes('Position',[.39 .05 .35 .33],'Color','none');
ax06 = axes('Position',[.82 .05 .13 .33],'Color','none');
%==========================================================================



MXA = MX;
MXB = MX';

axes(ax07);
%-------------------------------
bar( ([sum(MXA==-1);sum(MXA==0);sum(MXA==2);sum(MXA==3)]') ,'stacked');axis tight


axes(ax08);
%-------------------------------
bar( ([sum(MXB==-1);sum(MXB==0);sum(MXB==2);sum(MXB==3)]') ,'stacked');axis tight



axes(ax01);
%-------------------------------
ecdfhist(ECDFf, ECDFx, ECDFcenters)
h = findobj(gca,'Type','patch');
h.FaceColor = [.3 .3 .3];
hold on;
%plot(KERNELedges,pdCACOy,'LineWidth',6,'Color',[.4 .9 .4])
plot(KERNELedges,pdCTRLy,'LineWidth',6,'Color',[.4 .6 .9])
plot(KERNELedges,pdCASEy,'LineWidth',6,'Color',[.9 .4 .4])



axes(ax02);
%-------------------------------
histogram(xCTRL,P2_edges)
hold on
histogram(xCASE,P2_edges)




axes(ax03);
%-------------------------------
bar(BARX(:,1),BARY(:,1),0.45); 
hold on;
bar(BARX(:,2),BARY(:,2),0.45);




axes(ax04);
%-------------------------------
ksdensity(xCACO,xi);
xlabel('CASE')
ylabel('CTRL')
view([70,33])







axes(ax05);
%---PLOT  AREAGRAM ---------------------
area(NETDAT.XAREA,NETDAT.YAREA)
legend({'Case Miss','Ctrl Miss',...
        'Case Hit','Ctrl Hit'},...
        'Location','best');
ylabel('Count')
ax1.FontSize=16;






axes(ax06);
%---BAR GRAPH ---------------------
MU_CASE    = NETDAT.MU(1);
MU_CTRL    = NETDAT.MU(2);
MU_CACO    = NETDAT.MU(3);
MU_CASE_HI = NETDAT.MU(4);
MU_CTRL_HI = NETDAT.MU(5);
MU_CACO_HI = NETDAT.MU(6);

bar(([ MU_CASE_HI , MU_CTRL_HI , MU_CACO_HI ].*100),.5, 'FaceColor',[.31 .31 .31]); 
    hold on; 
bar(([ MU_CASE, MU_CTRL, MU_CACO ].*100),.20,'FaceColor',[.95 .85 .50]);


grid on; ylabel('Pct. Correct')
legend({'Top 25%','All'},'Location','Northwest','NumColumns',2)
ax06.YLim = [0 116]; 
ax06.YTick = [0 25 50 75 100]; 
ax06.XTickLabels = {'CASE','CTRL','ALL'};
ax06.XTickLabelRotation = 33;
ax06.FontSize=16;
pause(.2)
%----------------------------------------------------------------------











end