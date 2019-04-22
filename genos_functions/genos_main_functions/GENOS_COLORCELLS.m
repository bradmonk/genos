%% GENOS_COVAR_ANALYSIS.m

% ##################################################################
% CLEAR CONSOLE AND WORKSPACE
clc; close all; clear;
system('sudo purge')
cd(fileparts(which('GENOS_COLORCELLS.m')));



  MakeOrLoad = 'LOAD';
% MakeOrLoad = 'MAKE';




if strcmp(MakeOrLoad , 'MAKE' )


    f='GENOMICSDATA_EQUAL_PRO.mat'; 
    which(f)


    [ADMX,SNPCASE,SNPCTRL,PHE,CASEID,CTRLID,IDVX, ...
    ASYMX,ASYCASE,ASYCTRL,ADNN,caMX,coMX] = COVAR_PREP(f);


else

    
    GEN=load('ADSP.mat');
    load('GENOS_OUTPUT.mat');
    %load('GENOMIXPROCESSED_NNET.mat','CASEMX','CTRLMX')
    %CASEMX = CASEMX(ismember(CASEMX.SRR,CASEID) , :);
    %CTRLMX = CTRLMX(ismember(CTRLMX.SRR,CTRLID) , :);

end








clc; 
ADMX(1:4 , :)
ASYMX(1:4 , :)
caMX(1:7,1:7)
clearvars -except ADMX SNPCASE SNPCTRL PHE CASEID CTRLID IDVX CASEMX CTRLMX...
ASYMX ASYCASE ASYCTRL ADNN caMX coMX






%% GENERATE COLOCALIZATION MATRIX


V =  caMX';
VnCASE = sum(V,2);
JCASE = V*V';

V =  coMX';
VnCTRL = sum(V,2);
JCTRL = V*V';



clearvars -except ADMX SNPCASE SNPCTRL PHE CASEID CTRLID IDVX  CASEMX CTRLMX...
ASYMX ASYCASE ASYCTRL ADNN caMX coMX VnCASE JCASE VnCTRL JCTRL


close all; figure
imagesc( log(triu(JCASE) + tril(JCTRL)) );
title({'\fontsize{15}log(V*V'')     CTRL  \\  CASE '})
colormap(jet)

%% BOOTSTRAP PROBABILITY CASE GROUP
%{

PCASEh = JCASE.*0;

tot = 4000;
nboots = 1000;

for i = 1:size(JCASE,1)
% ai = 1:VnCASE(i);
parfor j = 1:size(JCASE,2)
    if j>=i
        a = zeros(1,tot);
        b = zeros(1,tot);

        a( 1:VnCASE(i) ) = 1;
        b( 1:VnCASE(j) ) = 1;

        cboot = bootstrp(nboots,@(x,y)sum(sum([x,y],2)>1),a',b');

        [f,x] = ecdf(cboot);
        h = JCASE(i,j) > x;
        k = find(h,1,'last');

        if k
        PCASEh(i,j) = 1-f(k);
        else
        PCASEh(i,j) = 1;
        end
    else
        PCASEh(i,j) = 1;
    end
end
disp(i)
end

% tic
% co = zeros(10000,1);
% for nn = 1:10000
%     a = a(randperm(tot));
%     b = b(randperm(tot));
%     co(nn) = sum(sum([a;b]) > 1);
% end
% toc
% close all;
% figure
% histogram(co)
% figure
% histogram(coboot)


%% BOOTSTRAP PROBABILITY CTRL GROUP


PCTRLh = JCTRL.*0;

tot = 4000;
nboots = 1000;

for i = 1:size(JCTRL,1)
% ai = 1:VnCTRL(i);
parfor j = 1:size(JCTRL,2)
    if j>=i
        a = zeros(1,tot);
        b = zeros(1,tot);

        a( 1:VnCTRL(i) ) = 1;
        b( 1:VnCTRL(j) ) = 1;

        cboot = bootstrp(nboots,@(x,y)sum(sum([x,y],2)>1),a',b');

        [f,x] = ecdf(cboot);
        h = JCTRL(i,j) > x;
        k = find(h,1,'last');

        if k
        PCTRLh(i,j) = 1-f(k);
        else
        PCTRLh(i,j) = 1;
        end
    else
        PCTRLh(i,j) = 1;
    end
end
disp(i)
end





%% PLOT TILE GRAPH

clc; close all;
fh1=figure('Units','normalized','OuterPosition',[.01 .04 .97 .85],'Color','k');
hax1 = axes('Position',[.04 .05 .45 .9],'Color','none');
hax2 = axes('Position',[.53 .05 .45 .9],'Color','none');


% axes(hax1)
% imagesc(  triu( 1-PCASEh, 1)  )
% colormap(hot)
% 
% 
% axes(hax2)
% imagesc(  triu( 1-PCTRLh, 1)  )
% colormap(gfp)


axes(hax1)
s1 = surf( triu( 1-PCASEh, 1) ); 
title('CASE'); view(-36,14)

axes(hax2)
s2 = surf( triu( 1-PCTRLh, 1) );
title('CTRL'); view(-36,14)


colormap(pink)

s1.EdgeColor = 'interp';
s2.EdgeColor = 'interp';

hax1.Color = 'none';
hax2.Color = 'none';


%}











%% COMPUTE PROBABILITIES USING HYPERGEOMETRIC DISTRIBUTION FUNCTION

% Y = hygecdf(X,M,K,N) 
% computes hypergeometric cdf at each X value using:
% M population size (CASE OR CTRL GROUP SIZE)
% K number of desired items in the population (SNPa COUNT)
% N number of samples drawn, without replacement (SNPb COUNT)
%
% NOTE: Comparing two variant loci, the locus with the greater variant
% count is referred to as SNPa, and the lesser is SNPb
%
% The parameters in M, K, and N must all be positive integers, 
% with N ? M. The values in X must be less than or equal to all 
% the parameter values.
%
% M = 4000;   % number of people
% K = 100;    % number of variant A
% N = 50;     % number of variant B
% X = 30;     % number of co-occurrances
%
% Y = hygecdf(X,M,K,N);
%
%
% S = 4000;   % number of people
% A = 100;    % number of variant A
% B = 50;     % number of variant B
% C = 30;     % number of co-occurrances
% 
% P = hygecdf(C,S,A,B);





% pCA=PCASEh;
% pCO=PCTRLh;


PCASEh = JCASE.*0;
PCASEb = JCASE.*0;
PCASEx = JCASE.*0;


groupN = size(CASEID,1);

for i = 1:size(JCASE,1)
for j = 1:size(JCASE,2)

    if j>i

        SNPS = [VnCASE(i) VnCASE(j)];

        S = groupN;     % number of people
        A = max(SNPS);  % number of variant A
        B = min(SNPS);  % number of variant B
        C = JCASE(i,j); % number of co-occurrances

        PCASEh(i,j) = 1 - hygecdf(C,S,A,B);
        PCASEb(i,j) = 1 - binocdf(C,B,A/S);
        PCASEx(i,j) = (C/B) * (1 - A/S) ;

    else
        PCASEh(i,j) = 1;
        PCASEb(i,j) = 1;
        PCASEx(i,j) = 1;
    end
end
disp(i)
end





PCTRLh = JCTRL.*0;
PCTRLb = JCTRL.*0;
PCTRLx = JCTRL.*0;

groupN = size(CTRLID,1);

for i = 1:size(JCTRL,1)
for j = 1:size(JCTRL,2)

    if j>i

        SNPS = [VnCTRL(i) VnCTRL(j)];

        S = groupN;     % number of people
        A = max(SNPS);  % number of variant A
        B = min(SNPS);  % number of variant B
        C = JCTRL(i,j); % number of co-occurrences

        PCTRLh(i,j) = 1 - hygecdf(C,S,A,B);
        PCTRLb(i,j) = 1 - binocdf(C,B,A/S);
        PCTRLx(i,j) = (C/B) * (1 - A/S) ;



    else
        PCTRLh(i,j) = 1;
        PCTRLb(i,j) = 1;
        PCTRLx(i,j) = 1;
    end
end
disp(i)
end



% PCTRLh(i,j) = 1 - hygecdf(C,S,A,B);
% PCTRLh(i,j) = 1 - binopdf(C,B,A/S);
% note this should really be...
% PCTRLh(i,j) = hygecdf(C-1,S,A,B,'upper');
% however, sometimes there are zero co-occurrances
% and the 'upper' option uses a more accurate but
% much slower algorithm


clearvars -except ADMX SNPCASE SNPCTRL PHE CASEID CTRLID IDVX CASEMX CTRLMX...
ASYMX ASYCASE ASYCTRL ADNN caMX coMX VnCASE JCASE VnCTRL JCTRL ...
PCASEh PCASEb PCASEx PCTRLh PCTRLb PCTRLx




%% PLOT PROBABILITY SURFACE GRAPHS


clc; close all;
fh1=figure('Units','normalized','OuterPosition',[.01 .04 .85 .93],'Color','w');
hax1 = axes('Position',[.05 .51 .40 .47],'Color','none');
hax2 = axes('Position',[.55 .51 .40 .47],'Color','none');
hax3 = axes('Position',[.05 .04 .40 .40],'Color','none');
hax4 = axes('Position',[.55 .04 .40 .40],'Color','none');

casePx=[]; ctrlPx=[];
for nn = 1:size(PCASEh,2)/2; casePx = [casePx; diag(PCASEh,nn)]; end
for nn = 1:size(PCTRLh,2)/2; ctrlPx = [ctrlPx; diag(PCTRLh,nn)]; end


axes(hax1); colormap(hax1,parula)
imagesc(  triu( -log(PCASEh) ,5)  )


axes(hax2); colormap(hax2,parula)
imagesc(  triu( -log(PCTRLh) ,5)  )


axes(hax3)
histogram(-log(casePx))
title('CASE -log(P) distribution')


axes(hax4)
histogram(-log(ctrlPx))
title('CTRL -log(P) distribution')



clearvars -except ADMX SNPCASE SNPCTRL PHE CASEID CTRLID IDVX CASEMX CTRLMX...
ASYMX ASYCASE ASYCTRL ADNN caMX coMX VnCASE JCASE VnCTRL JCTRL ...
PCASEh PCASEb PCASEx PCTRLh PCTRLb PCTRLx





%% CREATE TARGET TABLE

% find(isinf(-log(px)))
% diag(PCASEh)
% pCA = triu(PCASEh,1);

pCAlogic = PCASEh<1e-10;
pCAS = pCAlogic .* (PCASEh + 1e-15);

nCAS = pCAlogic .* JCASE;

[row, col, Pcov] = find(pCAS);
[wor, loc, Ncov] = find(nCAS);

GENES = cellstr(ASYMX.GENE);

v1GENE = GENES(row);
v2GENE = GENES(col);


v1CHR = ASYMX.CHR(row);
v2CHR = ASYMX.CHR(col);
v1POS = ASYMX.POS(row);
v2POS = ASYMX.POS(col);

v1chrpos = double(v1CHR).*1e10 + double(v1POS);
v2chrpos = double(v2CHR).*1e10 + double(v2POS);
v1v2dist = abs(v1chrpos - v2chrpos);

v1ALT = double(ASYMX.CASEALTS(row));
v2ALT = double(ASYMX.CASEALTS(col));


TCASE = table(Ncov,v1ALT,v2ALT,Pcov,v1GENE,v1CHR,v1POS,v2GENE,v2CHR,v2POS,v1v2dist);

COCASE = TCASE(v1v2dist>500000 , :);



clearvars -except ADMX SNPCASE SNPCTRL PHE CASEID CTRLID IDVX CASEMX CTRLMX...
ASYMX ASYCASE ASYCTRL ADNN caMX coMX VnCASE JCASE VnCTRL JCTRL ...
PCASEh PCASEb PCASEx PCTRLh PCTRLb PCTRLx ...
TCASE COCASE TCTRL COCTRL




pCOlogic = PCTRLh<1e-10;
pCOS = pCOlogic .* (PCTRLh + 1e-15);

nCOS = pCOlogic .* JCTRL;

[row, col, Pcov] = find(pCOS);
[wor, loc, Ncov] = find(nCOS);

GENES = cellstr(ASYMX.GENE);

v1GENE = GENES(row);
v2GENE = GENES(col);


v1CHR = ASYMX.CHR(row);
v2CHR = ASYMX.CHR(col);
v1POS = ASYMX.POS(row);
v2POS = ASYMX.POS(col);

v1chrpos = double(v1CHR).*1e10 + double(v1POS);
v2chrpos = double(v2CHR).*1e10 + double(v2POS);
v1v2dist = abs(v1chrpos - v2chrpos);

v1ALT = double(ASYMX.CASEALTS(row));
v2ALT = double(ASYMX.CASEALTS(col));

TCTRL = table(Ncov,v1ALT,v2ALT,Pcov,v1GENE,v1CHR,v1POS,v2GENE,v2CHR,v2POS,v1v2dist);

COCTRL = TCTRL(v1v2dist>500000 , :);



clearvars -except ADMX SNPCASE SNPCTRL PHE CASEID CTRLID IDVX CASEMX CTRLMX...
ASYMX ASYCASE ASYCTRL ADNN caMX coMX VnCASE JCASE VnCTRL JCTRL ...
PCASEh PCASEb PCASEx PCTRLh PCTRLb PCTRLx ...
TCASE COCASE TCTRL COCTRL








%% COMBINE CASE AND CONTROL MATRIX


% MAKE P(CACO)
pC = flipud(fliplr(PCTRLh)');
PCACO = pC .* PCASEh;


% MAKE NL(CACO)
pC = flipud(fliplr(PCTRLh)');
NLCACO = -log(pC .* PCASEh);
maxNP = max(NLCACO(~isinf(NLCACO)));
NLCACO(isinf(NLCACO)) = maxNP+1;



% MAKE N(CACO)
NCA = triu( JCASE, 0);
NCO = tril( JCTRL, -1);
NCACO = NCA+NCO;





% MAKE DF(CACO) ...NLCACO-NLCOCA
NLCOCA = flipud(fliplr(NLCACO)');
DFCOCA = NLCACO - NLCOCA;
DFCOCA(DFCOCA<0) = 0;
DFCACO = DFCOCA;



clc; close all;
fh1=figure('Units','normalized','OuterPosition',[.01 .04 .85 .93],'Color','w','Menubar','none');
hax1 = axes('Position',[.05 .53 .40 .44],'Color','none');
hax2 = axes('Position',[.55 .53 .40 .44],'Color','none');
hax3 = axes('Position',[.05 .04 .40 .44],'Color','none');
hax4 = axes('Position',[.55 .04 .40 .44],'Color','none');

axes(hax1);  imagesc( log(NCACO) );  title('log(N)')
axes(hax2);  imagesc(  PCACO     );  title('P')
axes(hax3);  imagesc(  NLCACO    );  title('-log(P)')
axes(hax4);  imagesc(  DFCACO    );  title(' P-P'' ')


clearvars -except ADMX SNPCASE SNPCTRL PHE CASEID CTRLID IDVX  CASEMX CTRLMX...
ASYMX ASYCASE ASYCTRL ADNN caMX coMX VnCASE JCASE VnCTRL JCTRL ...
PCASEh PCASEb PCASEx PCTRLh PCTRLb PCTRLx ...
TCASE COCASE TCTRL COCTRL PCACO NLCACO NCACO DFCACO




%% CREATE TABLE VARIABLES FOR CONCURRENCE: COUNT  PVAL  NLOGP


% CREATE ROW AND COLUMN NAMES USING GENE NAME + POS
tx = cellstr([ASYMX.GENE repmat('_',size(ASYMX,1),1) num2str(ASYMX.POS)]);
Rx1 = ' ';
Tx1 = '';
Rx2 = '-';
Tx2 = '';
tx = regexprep(tx,Rx1,Tx1);
txt = regexprep(tx,Rx2,Tx2);



% CREATE TABLE FOR COUNTS
CACO_COUNT = array2table(NCACO);
CACO_COUNT.Properties.VariableNames = txt;
CACO_COUNT.Properties.RowNames = txt;

% CREATE TABLE FOR PVAL
CACO_PVAL = array2table(PCACO);
CACO_PVAL.Properties.VariableNames = txt;
CACO_PVAL.Properties.RowNames = txt;


% CREATE TABLE FOR NLOGP
CACO_NLOGP = array2table(NLCACO);
CACO_NLOGP.Properties.VariableNames = txt;
CACO_NLOGP.Properties.RowNames = txt;




clearvars -except ADMX SNPCASE SNPCTRL PHE CASEID CTRLID IDVX  CASEMX CTRLMX...
ASYMX ASYCASE ASYCTRL ADNN caMX coMX VnCASE JCASE VnCTRL JCTRL ...
PCASEh PCASEb PCASEx PCTRLh PCTRLb PCTRLx ...
TCASE COCASE TCTRL COCTRL PCACO NLCACO NCACO DFCACO...
CACO_COUNT CACO_PVAL CACO_NLOGP txt DFCACO





%% CREATE HEATMAP CHARTS AND GIVEMESOMETHINGGOOOOD
clc; close all




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
fh = figure('Units','normalized','OuterPosition',[.01 .05 .95 .85],'Color','w');
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
t6 = uitable(th6,'Data',DFCACO,'Units','normalized','Position',[.01 .01 .98 .95]);
t6.ColumnName = ctxt; t6.RowName = rtxt;
t6.ColumnWidth = num2cell(ones(1,size(DFCACO,2)).*50);





%-----------------------------------------
%  HEATMAP PLOTS
%-----------------------------------------
ax5 = axes('Parent',th5,'Position',[.30 .01 .45 .97],'Color','none',...
           'YDir','reverse','XTick',[],'YTick',[],'XTickLabel','','YTickLabel','');
axis tight; hold on;


% IMHEAT = PCASEx;
IMHEAT = DFCACO;

t5 = imagesc(ax5,IMHEAT);


% get the colormap
C  = colormap;  
L  = size(C,1);



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

    st = sprintf('%0.2f', IMHEAT(nn) );

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





% clearvars -except ADMX SNPCASE SNPCTRL PHE CASEID CTRLID IDVX  CASEMX CTRLMX...
% ASYMX ASYCASE ASYCTRL ADNN caMX coMX VnCASE JCASE VnCTRL JCTRL ...
% PCASEh PCASEb PCASEx PCTRLh PCTRLb PCTRLx ...
% TCASE COCASE TCTRL COCTRL PCACO NLCACO NCACO DFCACO...
% CACO_COUNT CACO_PVAL CACO_NLOGP txt DFCACO











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


axes(ax8); colormap(ax8,parula); %pink%gfp
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



clearvars -except ADMX SNPCASE SNPCTRL PHE CASEID CTRLID IDVX  CASEMX CTRLMX...
ASYMX ASYCASE ASYCTRL ADNN caMX coMX VnCASE JCASE VnCTRL JCTRL ...
PCASEh PCASEb PCASEx PCTRLh PCTRLb PCTRLx ...
TCASE COCASE TCTRL COCTRL PCACO NLCACO NCACO DFCACO...
CACO_COUNT CACO_PVAL CACO_NLOGP txt DFCACO









%% SAVE DATASET FOR GENOS_COVAR_APP

save('GENOS_APPDATA.mat','ASYMX','CACO_COUNT','CACO_NLOGP','CACO_PVAL',...
'CTRLMX','DFCACO','JCASE','JCTRL','NCACO','NLCACO','PCACO','PCASEb',...
'PCASEh','PCASEx','PCTRLb','PCTRLh','PCTRLx','PHE','VnCASE','VnCTRL');







return
%% CREATE UIFIGURE (APP-STYLE GUI FEATURE)
close all; clc;

% CREATE ROW AND COLUMN NAMES USING GENE NAME + POS
%-----------------------------------------------------
tx = cellstr([ASYMX.GENE ...
             repmat('_',size(ASYMX,1),1) ...
             num2str(ASYMX.CASEALT) ...
             repmat(',',size(ASYMX,1),1) ...
             num2str(ASYMX.CTRLALT)  ]);
Rx1 = ' ';
Tx1 = '';
Rx2 = '-';
Tx2 = '';
tx = regexprep(tx,Rx1,Tx1);
rtxt = regexprep(tx,Rx2,Tx2);
ctxt = cellstr(ASYMX.GENE);



% DETERMINE POSITIONS OF FIGURE ELEMENTS
%-----------------------------------------------------
props = groot;
ss = props.ScreenSize;
figpos = [20 50 ss(3)*.95 ss(4)*.85];
tabpos = [1 1 figpos(3)-2 figpos(4)-2];
datpos = [1 1 tabpos(3)-2 tabpos(4)-25];



% CREATE FIGURE AND TAB GROUP PANELS
%-----------------------------------------------------
fi  = uifigure('Position',figpos);
th  = uitabgroup(fi,'Position',tabpos);
tb1 = uitab(th,'Title','Nco');
tb2 = uitab(th,'Title','P(Nco)');
tb3 = uitab(th,'Title','-logP(Nco)');
tb4 = uitab(th,'Title','Heatmap');



%-----------------------------------------
%  Nco TABLE   
%-----------------------------------------
t1 = uitable('Parent', tb1,'ColumnEditable', false,'Position',datpos,...   
             'BackgroundColor',[1 1 1; .9 .9 .9],'ColumnFormat',{'numeric'},...
             'Data', NCACO);

t1.ColumnName = ctxt; 
t1.RowName = rtxt;
t1.ColumnWidth = num2cell(ones(1,size(NCACO,2)).*50);
% t1.BackgroundColor=parula(size(NCACO,1));
% BackgroundColor',parula(size(NCACO,1)),...


%-----------------------------------------
%  P(Nco) TABLE
%-----------------------------------------
t2 = uitable('Parent', tb2,'ColumnEditable', false,'Position',datpos,...   
             'BackgroundColor',[1 1 1; .9 .9 .9],'ColumnFormat',{'numeric'},...
             'Data', PCACO);

t2.ColumnName = ctxt; 
t2.RowName = rtxt;
t2.ColumnWidth = num2cell(ones(1,size(PCACO,2)).*50);



%-----------------------------------------
% -logP(Nco) TABLE
%-----------------------------------------
t3 = uitable('Parent', tb3,'ColumnEditable', false,'Position',datpos,...   
             'BackgroundColor',[1 1 1; .9 .9 .9],'ColumnFormat',{'numeric'},...
             'Data', NLCACO);

t3.ColumnName = ctxt; 
t3.RowName = rtxt;
t3.ColumnWidth = num2cell(ones(1,size(NLCACO,2)).*50);



%-----------------------------------------
%  HEATMAP TABLE
%-----------------------------------------
ax1 = uiaxes(tb4);
% t4 = imagesc(NLCACO , 'Parent', tb4);

clearvars -except ADMX SNPCASE SNPCTRL PHE CASEID CTRLID IDVX  CASEMX CTRLMX...
ASYMX ASYCASE ASYCTRL ADNN caMX coMX VnCASE JCASE VnCTRL JCTRL ...
PCASEh PCASEb PCASEx PCTRLh PCTRLb PCTRLx ...
TCASE COCASE TCTRL COCTRL PCACO NLCACO NCACO DFCACO...
CACO_COUNT CACO_PVAL CACO_NLOGP txt DFCACO




















return
% ######################################################################
%%                PRINCIPAL COMPONENTS ANALYSIS
% ######################################################################


% SPLIT SUBJECT ADSP-ID TO ISOLATE COHORT
CASE_ADSPID = split(CASEMX.ADSPID,"-");
CTRL_ADSPID = split(CTRLMX.ADSPID,"-");


% CREATE TABLE OF ALL SUBJECTS WITH SRR COHORT SEX AGE APOE AD
CASE_ADSPMX = table(CASEMX.SRR,CASE_ADSPID,'VariableNames',{'SRR','ADSPID'});
CTRL_ADSPMX = table(CTRLMX.SRR,CTRL_ADSPID,'VariableNames',{'SRR','ADSPID'});
CASE_ADSPMX.COHORT = CASE_ADSPMX.ADSPID(:,2);
CTRL_ADSPMX.COHORT = CTRL_ADSPMX.ADSPID(:,2);

CASE_ADSPMX.SEX  = CASEMX.SEX;
CASE_ADSPMX.AGE  = CASEMX.AGE;
CASE_ADSPMX.APOE = CASEMX.APOE;
CASE_ADSPMX.AD   = CASEMX.AD;

CTRL_ADSPMX.SEX  = CTRLMX.SEX;
CTRL_ADSPMX.AGE  = CTRLMX.AGE;
CTRL_ADSPMX.APOE = CTRLMX.APOE;
CTRL_ADSPMX.AD   = CTRLMX.AD;




% FIND SUBSET OF SUBJECTS INCLUDED IN PCA MATRIX
[i,j] = ismember(CASE_ADSPMX.SRR,CASEID);
CASE_ADSPMX.USED = i;
CASE_ADSPMX.ORDER = j;

[i,j] = ismember(CTRL_ADSPMX.SRR,CTRLID);
CTRL_ADSPMX.USED = i;
CTRL_ADSPMX.ORDER = j;



% FIND SUBSET OF SUBJECTS INCLUDED IN PCA MATRIX
[i,j] = ismember(PHE.SRR,[CASEID;CTRLID]);
PHE.USED = i;
PHE.ORDER = j;



% CONCATINATE CASE AND CONTROL MATRICES INTO 'ADSPMX
ADSPMX = [CASE_ADSPMX; CTRL_ADSPMX];



% CREATE NUMERIC LABELS FOR COHORTS
[i,j,k] = unique(ADSPMX.COHORT);
ADSPMX.COHORTID = k;


% GIVE APOE ID CODES
ADSPMX.APOE(ADSPMX.APOE == 22) = 1;
ADSPMX.APOE(ADSPMX.APOE == 23) = 2;
ADSPMX.APOE(ADSPMX.APOE == 24) = 3;
ADSPMX.APOE(ADSPMX.APOE == 33) = 4;
ADSPMX.APOE(ADSPMX.APOE == 34) = 5;
ADSPMX.APOE(ADSPMX.APOE == 44) = 6;


% PHENAD = [CASE_ADSPMX; CTRL_ADSPMX];

% TRANSPOSE ADNN SO SUBJECTS ARE COLUMNS AND CONVERT TO UINT8
% PCAMX = double(ADNN(2:end,4:end))';


nCACO = (triu(JCASE) + tril(JCTRL));
LOGnCACO = log(triu(JCASE) + tril(JCTRL));



PCAMX = nCACO;



clearvars -except ADMX SNPCASE SNPCTRL PHE CASEID CTRLID IDVX  CASEMX CTRLMX...
ASYMX ASYCASE ASYCTRL ADNN caMX coMX VnCASE JCASE VnCTRL JCTRL ...
PCASEh PCASEb PCASEx PCTRLh PCTRLb PCTRLx ...
TCASE COCASE TCTRL COCTRL PCACO NLCACO NCACO DFCACO...
CACO_COUNT CACO_PVAL CACO_NLOGP txt DFCACO...
CASE_ADSPID CTRL_ADSPID CASE_ADSPMX CTRL_ADSPMX ADSPMX PCAMX




%% RUN PRINCIPLE COMPONENTS ANAYSIS
tic

ss = statset('pca')
ss.Display = 'iter';
ss.MaxIter = 100;
ss.TolFun = 1e4;
ss.TolX = 1e4;
ss.UseParallel = true;


[PCAcoeff,PCAscore,PCAlatent,PCAtsquared,PCAexplained] = pca(  PCAMX ,...
'Options',ss); % 'NumComponents',5,  'Algorithm','eig', Elapsed time is 70.459438 seconds.

toc



clearvars -except ADMX SNPCASE SNPCTRL PHE CASEID CTRLID IDVX  CASEMX CTRLMX...
ASYMX ASYCASE ASYCTRL ADNN caMX coMX VnCASE JCASE VnCTRL JCTRL ...
PCASEh PCASEb PCASEx PCTRLh PCTRLb PCTRLx ...
TCASE COCASE TCTRL COCTRL PCACO NLCACO NCACO DFCACO...
CACO_COUNT CACO_PVAL CACO_NLOGP txt DFCACO...
CASE_ADSPID CTRL_ADSPID CASE_ADSPMX CTRL_ADSPMX ADSPMX PCAMX...
PCAcoeff PCAscore PCAlatent PCAtsquared PCAexplained


% ######################################################################
%%            PLOT PCA RESULTS AND COLOR BY PHENOTYPE
% ######################################################################


%% PLOT ALL PCs IN 1-D
clc; close all;
fh1=figure('Units','normalized','OuterPosition',[.02 .06 .9 .8],'Color','w');
hax1 = axes('Position',[.08 .08 .40 .80],'Color','none');
hax2 = axes('Position',[.56 .08 .40 .80],'Color','none');


Tx = ASYMX.CHR;
Tn = unique(Tx);
Ts = 'CHROMOSOME';

axes(hax1); colormap(hax1, (prism(max(Tn))) );
ph1 = scatter( (1:numel(PCAcoeff(:,1))).*0+1  , PCAcoeff(:,1),50, Tx,'o');
hold on
ph2 = scatter( (1:numel(PCAcoeff(:,2))).*0+2  , PCAcoeff(:,2),50, Tx,'o');
hold on
ph3 = scatter( (1:numel(PCAcoeff(:,3))).*0+3  , PCAcoeff(:,3),900, Tx,'.');
hold on
ph4 = scatter( (1:numel(PCAcoeff(:,4))).*0+4  , PCAcoeff(:,4),900, Tx,'.');
hold on
ph5 = scatter( (1:numel(PCAcoeff(:,5))).*0+5  , PCAcoeff(:,5),900, Tx,'.');


hax1.XLim = [0 6];


text(1,max(ph1.YData),sprintf(['explains: \n ' num2str(round(PCAexplained(1),2)) '%% \n']),...
'Interpreter','tex','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',14);

text(2,max(ph2.YData),sprintf(['explains: \n ' num2str(round(PCAexplained(2),2)) '%% \n']),...
'Interpreter','tex','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',14);

text(3,max(ph3.YData),sprintf(['explains: \n ' num2str(round(PCAexplained(3),2)) '%% \n']),...
'Interpreter','tex','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',14);

text(4,max(ph4.YData),sprintf(['explains: \n ' num2str(round(PCAexplained(4),2)) '%% \n']),...
'Interpreter','tex','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',14);

text(5,max(ph5.YData),sprintf(['explains: \n ' num2str(round(PCAexplained(5),2)) '%% \n']),...
'Interpreter','tex','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',14);




xlabel('\bf PCA1 - PCA5');ylabel('\bf Eigenvector Coefficients');
cb=colorbar('Ticks',unique(Tx),'TickLabels',(unique(Tx)),'Direction','reverse');
cb.Label.String = Ts; cb.Label.FontSize = 18; cb.Label.Rotation = -90;
cb.Label.VerticalAlignment = 'baseline';





axes(hax2); colormap(hax2, (prism(max(Tn))) );

ph1 = scatter( (1:numel(PCAcoeff(:,1))).*0+1  , PCAcoeff(:,1),200, Tx,'o','MarkerEdgeAlpha',.02);
hold on
ph2 = scatter( (1:numel(PCAcoeff(:,2))).*0+2  , PCAcoeff(:,2),200, Tx,'o','MarkerEdgeAlpha',.02);
hold on
ph3 = scatter( (1:numel(PCAcoeff(:,3))).*0+3  , PCAcoeff(:,3),200, Tx,'o','MarkerEdgeAlpha',.02);
hold on
axis tight
hax2.XLim = [.5 3.5];
% hax2.XLim = [0 6];



clearvars -except ADMX SNPCASE SNPCTRL PHE CASEID CTRLID IDVX  CASEMX CTRLMX...
ASYMX ASYCASE ASYCTRL ADNN caMX coMX VnCASE JCASE VnCTRL JCTRL ...
PCASEh PCASEb PCASEx PCTRLh PCTRLb PCTRLx ...
TCASE COCASE TCTRL COCTRL PCACO NLCACO NCACO DFCACO...
CACO_COUNT CACO_PVAL CACO_NLOGP txt DFCACO...
CASE_ADSPID CTRL_ADSPID CASE_ADSPMX CTRL_ADSPMX ADSPMX PCAMX...
PCAcoeff PCAscore PCAlatent PCAtsquared PCAexplained




%% --- CHROMOSOME
clc; close all;
fh1=figure('Units','normalized','OuterPosition',[.02 .06 .97 .8],'Color','w');
hax1 = axes('Position',[.05 .07 .40 .85],'Color','none');
hax2 = axes('Position',[.52 .07 .40 .85],'Color','none');


Tx = ASYMX.CHR;
Tn = unique(Tx);
Ts = 'CHROMOSOME';


axes(hax1); colormap(hax1, (prism(max(Tn))) );
ph1 = scatter(PCAcoeff(:,1), PCAcoeff(:,2),500, Tx,'.');
xlabel('\bf PCA1');ylabel('\bf PCA2');

cb=colorbar('Ticks',Tn,'TickLabels',Tn,'Direction','reverse');
cb.Label.String = Ts; cb.Label.FontSize = 14; cb.Label.Rotation = -90;
cb.Label.VerticalAlignment = 'baseline';


axes(hax2); colormap(hax2, (prism(max(Tn))) );
ph2 = scatter3(PCAcoeff(:,1), PCAcoeff(:,2), PCAcoeff(:,3),500, Tx,'.');
view(-40,30)
xlabel('\bf PCA1');ylabel('\bf PCA2');zlabel('\bf PCA3');



clearvars -except ADMX SNPCASE SNPCTRL PHE CASEID CTRLID IDVX  CASEMX CTRLMX...
ASYMX ASYCASE ASYCTRL ADNN caMX coMX VnCASE JCASE VnCTRL JCTRL ...
PCASEh PCASEb PCASEx PCTRLh PCTRLb PCTRLx ...
TCASE COCASE TCTRL COCTRL PCACO NLCACO NCACO DFCACO...
CACO_COUNT CACO_PVAL CACO_NLOGP txt DFCACO...
CASE_ADSPID CTRL_ADSPID CASE_ADSPMX CTRL_ADSPMX ADSPMX PCAMX...
PCAcoeff PCAscore PCAlatent PCAtsquared PCAexplained


%%





