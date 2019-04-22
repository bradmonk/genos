function [] = pcaplotter(PC,PHECOLORS,cbTITLE,CMAP,cbTICKS,cbTICKLABELS,...
                PClast,AZEL,do2D,doAnimate,doSave,varargin)


% keyboard

if do2D
fh1=figure('Units','normalized','OuterPosition',[.01 .05 .98 .86],'Color','w');
hax1 = axes('Position',[.04 .55 .22 .4],'Color','none');
hax2 = axes('Position',[.28 .55 .22 .4],'Color','none');
hax3 = axes('Position',[.04 .04 .22 .4],'Color','none');
hax4 = axes('Position',[.28 .04 .22 .4],'Color','none');
hax5 = axes('Position',[.60 .15 .35 .7],'Color','none');

if strcmp(cbTITLE,'AGE')
    doCLabs = 0;
else
    doCLabs = 1;
end


    axes(hax1); colormap(hax1, CMAP); hold on;
ph1 = scatter(PC(:,1), PC(:,2),50, PHECOLORS,'filled','MarkerFaceAlpha',.2);
    xlabel('\bf PC1');ylabel('\bf PC2'); title('PC1 & PC2')
    if ~doCLabs; cb=colorbar('Direction','reverse'); else
    cb=colorbar('Ticks',cbTICKS,'TickLabels',cbTICKLABELS,'Direction','reverse'); end
    cb.Label.String = cbTITLE; cb.Label.FontSize = 14; cb.Label.Rotation = -90;
    cb.Label.VerticalAlignment = 'baseline'; axis off;


    axes(hax2); colormap(hax2, CMAP); hold on;
ph2 = scatter(PC(:,1), PC(:,3),50, PHECOLORS,'filled','MarkerFaceAlpha',.2);
    xlabel('\bf PC1');ylabel('\bf PC3'); title('PC1 & PC3')
    if ~doCLabs; cb=colorbar('Direction','reverse'); else
    cb=colorbar('Ticks',cbTICKS,'TickLabels',cbTICKLABELS,'Direction','reverse'); end
    cb.Label.String = cbTITLE; cb.Label.FontSize = 14; cb.Label.Rotation = -90;
    cb.Label.VerticalAlignment = 'baseline'; axis off;

    axes(hax3); colormap(hax3, CMAP); hold on;
ph3 = scatter(PC(:,2), PC(:,3),50, PHECOLORS,'filled','MarkerFaceAlpha',.2);
    xlabel('\bf PC2');ylabel('\bf PC3'); title('PC2 & PC3')
    if ~doCLabs; cb=colorbar('Direction','reverse'); else
    cb=colorbar('Ticks',cbTICKS,'TickLabels',cbTICKLABELS,'Direction','reverse'); end
    cb.Label.String = cbTITLE; cb.Label.FontSize = 14; cb.Label.Rotation = -90;
    cb.Label.VerticalAlignment = 'baseline'; axis off;

    axes(hax4); colormap(hax4, CMAP); hold on;
ph4 = scatter(PC(:,PClast(1)), PC(:,PClast(2)),50, PHECOLORS,'filled','MarkerFaceAlpha',.2);
    xlabel('\bf PCx');ylabel('\bf PCy'); title('PCx & PCy')
    if ~doCLabs; cb=colorbar('Direction','reverse'); else
    cb=colorbar('Ticks',cbTICKS,'TickLabels',cbTICKLABELS,'Direction','reverse'); end
    cb.Label.String = cbTITLE; cb.Label.FontSize = 14; cb.Label.Rotation = -90;
    cb.Label.VerticalAlignment = 'baseline'; axis off;




    axes(hax5); colormap(hax5, CMAP); hold on;
ph5 = scatter3(PC(:,1), PC(:,2), PC(:,3),50, PHECOLORS,'filled','MarkerFaceAlpha',.2);
    xlabel('\bf PC1');ylabel('\bf PC2');zlabel('\bf PC3');
    axis on; grid on; view([AZEL(1),AZEL(2)]);


if doAnimate
    axes(hax5); colormap(hax5, CMAP); hold on;
ph5 = scatter3(PC(:,1), PC(:,2), PC(:,3),5, PHECOLORS);
    xlabel('\bf PCA1');ylabel('\bf PCA2');zlabel('\bf PCA3');
    view([-60,35]); grid off; axis vis3d off; pause(1)
    for i=1:36
        camorbit(10,0,'camera',[0 0 1])
        drawnow
    end
    cla
    axes(hax5); colormap(hax5, CMAP); hold on;
ph5 = scatter3(PC(:,1), PC(:,2), PC(:,3),50, PHECOLORS,'filled','MarkerFaceAlpha',.2);
    xlabel('\bf PC1');ylabel('\bf PC2');zlabel('\bf PC3');
    axis on; grid on; view([AZEL(1),AZEL(2)]);
end



%%
%---------- EXPORT IMAGE TO DESKTOP FOLDER -----------------------------
if doSave
dt=char(datetime(datetime,'Format','yyyy-MM-dd-HH-mm-ss'));
print(['/Users/bradleymonk/Desktop/PCA/PCA_' cbTITLE dt '.png'],'-dpng');
end
%-----------------------------------------------------------------------

end











%% PLOT PCA --- ALL PCs IN 1D -------------------------------------

% keyboard


if ~do2D

fh1=figure('Units','normalized','OuterPosition',[.01 .05 .98 .86],'Color','w');
hax1 = axes('Position',[.08 .08 .40 .80],'Color','none');
hax2 = axes('Position',[.56 .08 .40 .80],'Color','none');


PE = varargin{1};
XXcbTITLE      = cbTITLE;
XXPHECOLORS    = PHECOLORS;
XXcbTICKLABELS = cbTICKLABELS;
XXcbTICKS      = cbTICKS;
XXCMAP         = CMAP;

%-----------------------------
TITL         = XXcbTITLE.A;
PHECOLORS    = XXPHECOLORS.A;
cbTICKLABELS = XXcbTICKLABELS.A;
cbTICKS      = XXcbTICKS.A;
CMAP         = XXCMAP.A;
%-----------------------------

    nPC = zeros(size(PC(:,1)));
    axes(hax1); colormap(hax1, CMAP); hold on;
ph1 = scatter( nPC+1  , PC(:,1), 50, PHECOLORS,'o','MarkerEdgeAlpha',.02); hold on
ph2 = scatter( nPC+2  , PC(:,2), 50, PHECOLORS,'o','MarkerEdgeAlpha',.02); hold on
ph3 = scatter( nPC+3  , PC(:,3), 50, PHECOLORS,'o','MarkerEdgeAlpha',.02); hold on
ph4 = scatter( nPC+4  , PC(:,4), 50, PHECOLORS,'o','MarkerEdgeAlpha',.02); hold on
ph5 = scatter( nPC+5  , PC(:,PClast(2)), 50, PHECOLORS,'o','MarkerEdgeAlpha',.02);

    xlabel('\bf PCA Components 1-5 ');ylabel('\bf Eigenvector Coefficients');
    cb=colorbar('Ticks',cbTICKS,'TickLabels',cbTICKLABELS,'Direction','reverse');
    cb.Label.String = TITL; cb.Label.FontSize = 14; cb.Label.Rotation = -90;
    cb.Label.VerticalAlignment = 'baseline'; hax1.XLim = [0 6];

text(1,max(ph1.YData),sprintf(['explains: \n ' num2str(round(PE(1),2)) '%% \n']),...
'Interpreter','tex','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',14);
text(2,max(ph2.YData),sprintf(['explains: \n ' num2str(round(PE(2),2)) '%% \n']),...
'Interpreter','tex','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',14);
text(3,max(ph3.YData),sprintf(['explains: \n ' num2str(round(PE(3),2)) '%% \n']),...
'Interpreter','tex','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',14);
text(4,max(ph4.YData),sprintf(['explains: \n ' num2str(round(PE(4),2)) '%% \n']),...
'Interpreter','tex','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',14);
text(5,max(ph5.YData),sprintf(['explains: \n ' num2str(round(PE(PClast(2)),2)) '%% \n']),...
'Interpreter','tex','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',14);


%-----------------------------
% TITL         = 'Case-vs-Ctrl';
% PHECOLORS    = PHE.AD;
% cbTICKLABELS = unique(PHE.AD);
% cbTICKS      = 1:numel(cbTICKLABELS);
% CM = jet;
% CMAP        = [CM(17,:); CM(end-16,:)];
%-----------------------------
%-----------------------------
TITL         = XXcbTITLE.B;
PHECOLORS    = XXPHECOLORS.B;
cbTICKLABELS = XXcbTICKLABELS.B;
cbTICKS      = XXcbTICKS.B;
CMAP         = XXCMAP.B;
%-----------------------------


    nPC = zeros(size(PC(:,1)));
    axes(hax2); colormap(hax2, CMAP); hold on;
ph1 = scatter( nPC+1  , PC(:,1), 50, PHECOLORS,'o','MarkerEdgeAlpha',.02);
hold on
ph2 = scatter( nPC+2  , PC(:,2), 50, PHECOLORS,'o','MarkerEdgeAlpha',.02);
hold on
ph3 = scatter( nPC+3  , PC(:,3), 50, PHECOLORS,'o','MarkerEdgeAlpha',.02);
hold on
ph4 = scatter( nPC+4  , PC(:,4), 50, PHECOLORS,'o','MarkerEdgeAlpha',.02);
hold on
ph5 = scatter( nPC+5  , PC(:,PClast(2)), 50, PHECOLORS,'o','MarkerEdgeAlpha',.02);



xlabel('\bf PCA Components 1-5 ');ylabel('\bf Eigenvector Coefficients');
cb=colorbar('Ticks',cbTICKS,'TickLabels',cbTICKLABELS,'Direction','reverse');
cb.Label.String = TITL; cb.Label.FontSize = 18; cb.Label.Rotation = -90;
cb.Label.VerticalAlignment = 'baseline'; hax2.XLim = [0 6];



%---------- EXPORT IMAGE TO DESKTOP FOLDER -----------------------------
if doSave
dt=char(datetime(datetime,'Format','yyyy-MM-dd-HH-mm-ss'));
print(['/Users/bradleymonk/Desktop/PCA/PCA_1D_' dt '.png'],'-dpng');
end
%-----------------------------------------------------------------------
end




end