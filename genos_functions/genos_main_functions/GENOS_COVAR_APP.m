function [] = GENOS_COVAR_APP()
%% COVARAPP: DISPLAYS CO-OCCURRENCE VARIANT TABLES
clc; close all; clear;


APPDATA = 'GENOS_APPDATA';

load(APPDATA,'ASYMX','CACO_COUNT','CACO_NLOGP','CACO_PVAL',...
'CTRLMX','DFCACO','JCASE','JCTRL','NCACO','NLCACO','PCACO','PCASEb',...
'PCASEh','PCASEx','PCTRLb','PCTRLh','PCTRLx','PHE','VnCASE','VnCTRL');

% load('GENOS_APPDATA.mat')


%% SET DESIRED VARIABLE TO DISPLAY IN THE HEATMAP TABLE
% IN ALL TABLES:   BOTTOM-LEFT:CTRL   TOP-RIGHT:CASE


% OPTIONS FOR DISPLAY IN THE HEATMAP TABLE...

% IMHEAT = NCACO;           % CO-OCCURENCE COUNTS
% IMHEAT = PCACO;           % P(NCACO)
% IMHEAT = NLCACO;          % -log(PCACO)
% IMHEAT = DFCACO;            % NLCACO normalized to opposing group


IMHEAT = DFCACO;



%% MAKE MIDDLE DIAGONAL MORE EASILY VISIBLE

NCACO(diag(diag(NCACO,0).*0+1)==1) = 0;

PCACO(diag(diag(PCACO,0).*0+1)==1) = 0;

NLCACO = round(NLCACO);

DFCACO = round(DFCACO);






%% CREATE TABLES AND PLOTS


% CREATE ROW AND COLUMN NAMES USING GENE NAME + #of ALTS
%--------------------------------------------------------
sz=size(ASYMX,1);
tx = cellstr([num2str([1:sz]') repmat('_',sz,1) ASYMX.GENE repmat('_',sz,1) num2str(ASYMX.CASEALT) ...
              repmat(',',size(ASYMX,1),1) num2str(ASYMX.CTRLALT)  ]);
Rx1 = ' '; Tx1 = ''; Rx2 = '-'; Tx2 = '';
tx = regexprep(tx,Rx1,Tx1);
rtxt = regexprep(tx,Rx2,Tx2);
ctxt = regexprep(tx,Rx2,Tx2);
% ctxt = cellstr(ASYMX.GENE);






% CREATE FIGURE AND TAB GROUP PANELS
%-----------------------------------------------------
fh = figure('Units','normalized','OuterPosition',[.01 .05 .95 .85],'Color','w','Menubar','none');
th = uitabgroup(fh,'Position',[0.01 0.01 0.97 0.97]);
th1 = uitab(th,'Title','Nco','BackgroundColor',[1.00 1.00 1.00]);
th2 = uitab(th,'Title','P(Nco)','BackgroundColor',[1.00 1.00 1.00]);
th3 = uitab(th,'Title','-logP','BackgroundColor',[1.00 1.00 1.00]);
th4 = uitab(th,'Title','Heatmap table','BackgroundColor',[1.00 1.00 1.00]);
th5 = uitab(th,'Title','Heatmap plots','BackgroundColor',[1.00 1.00 1.00]);
th6 = uitab(th,'Title','DIFF-logP','BackgroundColor',[1.00 1.00 1.00]);


%-----------------------------------------
%  Nco TABLE   
%-----------------------------------------
t1 = uitable(th1,'Data',NCACO, 'Units','normalized','Position',[.01 .01 .98 .95]);
t1.ColumnName = ctxt; t1.RowName = rtxt; t1.FontSize=10;
t1.ColumnWidth = num2cell(ones(1,size(NCACO,2)).*50);


%-----------------------------------------
%  P(Nco) TABLE   
%-----------------------------------------
t2 = uitable(th2,'Data',PCACO, 'Units','normalized','Position',[.01 .01 .98 .95]);
t2.ColumnName = ctxt; t2.RowName = rtxt;
t2.ColumnWidth = num2cell(ones(1,size(PCACO,2)).*50);


%-----------------------------------------
%  -logP TABLE   
%-----------------------------------------
t3 = uitable(th3,'Data',NLCACO,'Units','normalized','Position',[.01 .01 .98 .95]);
t3.ColumnName = ctxt; t3.RowName = rtxt;
t3.ColumnWidth = num2cell(ones(1,size(NLCACO,2)).*50);




%-----------------------------------------
%  DIFF(-logP) TABLE   
%-----------------------------------------
t6 = uitable(th6,'Data',round(DFCACO),'Units','normalized','Position',[.01 .01 .98 .95]);
t6.ColumnName = ctxt; t6.RowName = rtxt;
t6.ColumnWidth = num2cell(ones(1,size(DFCACO,2)).*50);





%-----------------------------------------
%  HEATMAP PLOTS
%-----------------------------------------
ax5 = axes('Parent',th5,'Position',[.30 .01 .45 .97],'Color','none',...
           'YDir','reverse','XTick',[],'YTick',[],'XTickLabel','','YTickLabel','');
colormap(parula)
axis tight; hold on;


% IMHEAT = PCASEx;
% IMHEAT = DFCACO;



t5 = imagesc(ax5,IMHEAT);


% get the colormap
C  = colormap(ax5);
L  = size(C,1);
colormap(parula)





% COLOR GENERATOR FUNCTION FOR TABLE HEATMAPS
%-----------------------------------------------------
colergen = @(color,text) ['<html><table border=0 width=300 bgcolor=',...
           color,'><TR><TD>',text,'</TD></TR> </table></html>'];




%-----------------------------------------
%  HEATMAP TABLES
%-----------------------------------------

% Scale the matrix to the range of the map.
Gs = round(interp1(linspace(min(IMHEAT(:)),max(IMHEAT(:)),L),1:L,IMHEAT));

% Make RGB image from scaled.
iRGB = reshape(C(Gs,:),[size(Gs) 3]); 

iRGB = round(iRGB.*255);

HEAT = cell(size(IMHEAT));

sz = size(IMHEAT);

for nn = 1:numel(HEAT)

    [I,J] = ind2sub(sz,nn);

    rgb = squeeze(iRGB(I,J,:))';
    if I==J
    rgb=[0 0 0];
    end

    hex(:,2:7) = reshape(sprintf('%02X',rgb.'),6,[]).'; 
    hex(:,1) = '#';

    %st = sprintf('%0.2f', IMHEAT(nn) );
    st = sprintf('%0.0f', IMHEAT(nn) );

    HEAT{nn} = colergen(  hex  ,  st  );
end

% sprintf('%0.2f',IMHEAT(3,5))
% for nn = 1:numel(HEAT)
%     HEAT{nn} = colergen(  '#FF0000'  ,  num2str(NLCACO(nn))  );
% end

t4 = uitable(th4,'Data',HEAT,'Units','normalized','Position',[.01 .01 .98 .95]);
t4.ColumnName = ctxt; t4.RowName = rtxt;
t4.ColumnWidth = num2cell(ones(1,size(HEAT,2)).*50);

t4.FontWeight= 'bold';
t4.ForegroundColor = 'red';







%-----------------------------------------
%  ADD MORE HEATMAP PLOTS
%-----------------------------------------
ax7  = axes('Parent',th5,'Position',[.03 .53 .25 .40],'Color','none',...
           'YDir','reverse','XTick',[],'YTick',[],'XTickLabel','','YTickLabel','');
            axis tight; hold on;
ax8  = axes('Parent',th5,'Position',[.03 .01 .25 .40],'Color','none',...
           'YDir','reverse','XTick',[],'YTick',[],'XTickLabel','','YTickLabel','');
            axis tight; hold on;
ax9  = axes('Parent',th5,'Position',[.77 .54 .20 .40],'Color','none');
            axis tight; hold on;
ax10 = axes('Parent',th5,'Position',[.77 .04 .20 .40],'Color','none');
            axis tight; hold on;



% axes(ax7); colormap(ax7,parula)
% imagesc(  triu( -log(PCASEh) ,5)  )
% axes(ax8); colormap(ax8,parula)
% imagesc(  triu( -log(PCTRLh) ,5)  )


axes(ax7); colormap(ax7,parula);
s1 = surf( triu( -log(PCASEh) ,5) ); 
title('CASE -logP'); view(25,30)
s1.EdgeColor = 'interp';
ax7.Color = 'none';


axes(ax8); colormap(ax8,parula); %pink%parula
s2 = surf( triu( -log(PCTRLh) ,5) );
title('CTRL -logP'); view(-117,25)
s2.EdgeColor = 'interp';
ax8.Color = 'none';
ax8.ZLim = ax7.ZLim;




casePx = [];
for nn = 1:size(PCASEh,2)/2
    casePx = [casePx; diag(PCASEh,nn)];
end

ctrlPx = [];
for nn = 1:size(PCTRLh,2)/2
    ctrlPx = [ctrlPx; diag(PCTRLh,nn)];
end


axes(ax9)
histogram(-log(casePx))
title('CASE -log(P) distribution')

axes(ax10)
histogram(-log(ctrlPx))
title('CTRL -log(P) distribution')

end
