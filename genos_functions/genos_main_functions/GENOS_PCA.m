% ######################################################################
%%                PRINCIPAL COMPONENTS ANALYSIS
% ######################################################################



%% START HERE IF YOU NEED TO GENERATE P-VALUES AND A PCA MATRIX
% IF YOU HAVE THOSE ALREADY SKIP DOWN

clc; close all; clear; rng('shuffle')
cd(fileparts(which('GENOS_PCA.m')));


MATDATA = 'ADSP.mat';
which(MATDATA)
load(MATDATA)

clearvars -except ADSP



%% CARBON COPY MAIN VARIABLES FROM ADSP.STRUCT

LOCI = ADSP.LOCI(:,1:17);
CASE = ADSP.CASE;
CTRL = ADSP.CTRL;
PHEN = ADSP.PHEN;

clearvars -except ADSP LOCI CASE CTRL PHEN





%###############################################################
%%       DETERMINE WHICH PARTICIPANTS TO KEEP
%###############################################################


% TAG GOOD PHEN
GOODCOH = [1 2 6 9 10 11 12 13 19 20]';
MEHHCOH = [5 7 23 24]'; % 7 IS HISPANIC STUDY


[ai,bi] = ismember(PHEN.COHORTNUM,GOODCOH);
% [ai,bi] = ismember(PHEN.COHORTNUM,[GOODCOH; MEHHCOH]);

PHEN.GOODCOH = ai;


% PUT GOOD PHEN INTO SEPARATE SET
PHE = PHEN(PHEN.GOODCOH==1,:);     % 8824 TOTAL PEOPLE




% ONLY KEEP PEOPLE WITH A REASONABLE NUMBER OF VARIANTS
PHE = PHE(PHE.TOTvars>14000,:);


PHECASE = PHE(PHE.AD==1,:);
PHECTRL = PHE(PHE.AD==0,:);


clearvars -except ADSP LOCI CASE CTRL PHEN PHE PHECASE PHECTRL



%###############################################################
%%          COUNT NUMBER OF VARIANTS PER LOCI
%###############################################################

% The varsum() function will go through each known variant loci
% and check whether anyone's SRR ID from your subset of IDs match
% all known SRR IDs for that loci. It will then sum the total
% number of alleles (+1 for hetzy-alt, +2 for homzy-alt) for each
% loci and return the totals.

UNSNP = ADSP.UNSNP;


[CASEN, CTRLN] = varsum(CASE, PHECASE.SRR, CTRL, PHECTRL.SRR);

[CASEUN, CTRLUN] = uncalledsum(UNSNP, PHECASE.SRR, PHECTRL.SRR);



% SAVE COUNTS AS NEW TABLE COLUMNS
LOCI.CASEREFS = (numel(PHECASE.SRR)*2) - (CASEUN.*2) - CASEN;
LOCI.CTRLREFS = (numel(PHECTRL.SRR)*2) - (CTRLUN.*2) - CTRLN;
LOCI.CASEALTS = CASEN;
LOCI.CTRLALTS = CTRLN;


clearvars -except ADSP LOCI CASE CTRL PHEN PHE PHECASE PHECTRL







%###############################################################
%%               COMPUTE FISHER'S P-VALUE
%###############################################################


% COMPUTE FISHERS STATISTICS FOR THE TRAINING GROUP
[FISHP, FISHOR] = fishp_mex(LOCI.CASEREFS,LOCI.CASEALTS,...
                            LOCI.CTRLREFS,LOCI.CTRLALTS);

LOCI.FISHPS  = FISHP;
LOCI.FISHORS = FISHOR;


clearvars -except ADSP LOCI CASE CTRL PHEN PHE PHECASE PHECTRL





%% MAKE LATEST COUNTS THE MAIN TABLE STATS

LOCI.CASEREF = LOCI.CASEREFS;
LOCI.CTRLREF = LOCI.CTRLREFS;
LOCI.CASEALT = LOCI.CASEALTS;
LOCI.CTRLALT = LOCI.CTRLALTS;
LOCI.FISHP   = LOCI.FISHPS;
LOCI.FISHOR  = LOCI.FISHORS;






%% CLEAN-UP VARIANT TABLE


LOCI.GENE = string(LOCI.GENE);
LOCI.GENE = replace(LOCI.GENE," ","");

LOCI.AAN = string(LOCI.AAN);
LOCI.AAN = replace(LOCI.AAN," ","");
LOCI.AAN = replace(LOCI.AAN,"NA","NaN;NaN");
LOCI.AAN = replace(LOCI.AAN,";",",");
AAN = LOCI.AAN;
AAN = split(AAN,",");
AAN = str2double(AAN);
LOCI.AAN = AAN;

clc; clearvars -except ADSP LOCI CASE CTRL USNP PHEN PHE PHECASE PHECTRL
disp(LOCI(1:9,:))




%% SORT CONTAINERS BY FISHP (SMALLEST TO LARGEST)

USNP = ADSP.UNSNP;

[X,i] = sort(LOCI.FISHP);

LOCI  = LOCI(i,:);
CASE  = CASE(i);
CTRL  = CTRL(i);
USNP  = USNP(i);
LOCI.VID = (1:size(LOCI,1))';


clc; clearvars -except ADSP LOCI CASE CTRL USNP PHEN PHE PHECASE PHECTRL
disp(LOCI(1:9,:))








%% STORE VARIABLES FOR PCA/TSNE AS 'AMX'

AMX         = LOCI;
AMXCASE     = CASE;
AMXCTRL     = CTRL;
AMXUSNP     = USNP;


clearvars -except ADSP LOCI CASE CTRL USNP PHEN PHE PHECASE PHECTRL...
AMX AMXCASE AMXCTRL AMXUSNP





%% FILTER VARIANTS BASED ALT > REF

% PASS = (AMX.CASEREF > AMX.CASEALT./1.5) | (AMX.CTRLREF > AMX.CTRLALT./1.5);
% sum(~PASS)
% 
% AMX      = AMX(PASS,:);
% AMXCASE  = AMXCASE(PASS);
% AMXCTRL  = AMXCTRL(PASS);
% AMXUSNP  = AMXUSNP(PASS);
% AMX.VID  = (1:size(AMX,1))';
% 
% 
% 
% 
% clearvars -except ADSP LOCI CASE CTRL USNP PHEN PHE PHECASE PHECTRL...
% AMX AMXCASE AMXCTRL AMXUSNP





%% TAKE THE TOP N NUMBER OF VARIANTS


N = 8000;


AMX      = AMX(1:N,:);
AMXCASE  = AMXCASE(1:N);
AMXCTRL  = AMXCTRL(1:N);
AMXUSNP  = AMXUSNP(1:N);
AMX.VID  = (1:size(AMX,1))';

fprintf('\n %.0f final loci count \n\n',size(AMX,1))

clearvars -except ADSP LOCI CASE CTRL USNP PHEN PHE PHECASE PHECTRL...
AMX AMXCASE AMXCTRL AMXUSNP








%% MAKE  RECTANGLE  NN  VARIANT MATRIX


[ADNN, caMX, coMX] = unnetpca(AMX,AMXCASE,AMXCTRL,AMXUSNP,PHE);


clearvars -except ADSP LOCI CASE CTRL USNP PHEN PHE PHECASE PHECTRL...
AMX AMXCASE AMXCTRL AMXUSNP ADNN







%% RANDOMIZE ADNN AND REORDER PHE TO MATCH ADNN

ADL = ADNN(1,:);
ADN = ADNN(2:end,:);
i = randperm(size(ADN,1));
ADN = ADN(i,:);
ADNN = [ADL;ADN];


[i,j] = ismember(PHE.SRR, ADN(:,1) );
PHE.USED = i;
PHE.ORDER = j;
PHE = PHE(PHE.USED,:);
PHE = sortrows(PHE,'ORDER');



PCAMX = ADNN(2:end,4:end);


clc; clearvars -except ADSP LOCI CASE CTRL USNP PHEN PHE PHECASE PHECTRL...
AMX AMXCASE AMXCTRL AMXUSNP ADNN PCAMX  



%######################################################################
%%       PCA : PRINCIPAL COMPONENTS ANALYSIS
%######################################################################


ss = statset('pca')
% ss.Display = 'iter';
ss.Display = 'none';
ss.MaxIter = 100;
ss.TolFun = 1e4;
ss.TolX = 1e4;
ss.UseParallel = true;

% 'NumComponents',5,  'Algorithm','eig', Elapsed time is 70.459438 seconds.
% [PCAC,PCAS,PCAlatent,PCAtsquared,PCAE] = pca(  PCAMX' ,'Options',ss); 


% PCAMX(PCAMX>0) = 1;

[PCAC,PCAS,~,~,PCAE] = pca(  PCAMX' , 'Options',ss);


clc; clearvars -except ADSP LOCI CASE CTRL USNP PHEN PHE PHECASE PHECTRL...
AMX AMXCASE AMXCTRL AMXUSNP ADNN PCAMX PCAC PCAS PCAE  







%% RECONSTRUCT ORIGINAL DATA (OPTIONAL TEST)



% MEAN DEVIATE (CENTER) DATA
mu = mean(PCAMX);
Xi = bsxfun(@minus,PCAMX,mu);


% Reconstruct the centered data
Xj = PCAS*PCAC';

% Reconstruct original data
PCAMX_REDUX = round(Xj' + mu);


clc; clearvars -except ADSP LOCI CASE CTRL USNP PHEN PHE PHECASE PHECTRL...
AMX AMXCASE AMXCTRL AMXUSNP ADNN PCAMX PCAC PCAS PCAE PCAMX_REDUX








%######################################################################
%%       PLOT PCA RESULTS AND COLOR BY PHENOTYPE
%######################################################################


%% PLOT PCA --- ALL PCs IN 1D -------------------------------------
clc; close all;
fh1=figure('Units','normalized','OuterPosition',[.02 .06 .9 .8],'Color','w');
hax1 = axes('Position',[.08 .08 .40 .80],'Color','none');
hax2 = axes('Position',[.56 .08 .40 .80],'Color','none');


Tx = PHE.COHORTNUM;
Tn = unique(Tx);
Ts = 'Cohort';

% axes(hax1); colormap(hax1, (prism(max(Tn))) );
axes(hax1); colormap(hax1, flipud(jet(max(Tn))) );

ph1 = scatter( zeros(size(PCAC(:,1)))+1  , PCAC(:,1),50, Tx,'o');
hold on
ph2 = scatter( zeros(size(PCAC(:,2)))+2  , PCAC(:,2),50, Tx,'o');
hold on
ph3 = scatter( zeros(size(PCAC(:,3)))+3  , PCAC(:,3),50, Tx,'o');
hold on
ph4 = scatter( zeros(size(PCAC(:,4)))+4  , PCAC(:,4),900, Tx,'.');
hold on
ph5 = scatter( zeros(size(PCAC(:,5)))+5  , PCAC(:,5),900, Tx,'.');


hax1.XLim = [0 6];


text(1,max(ph1.YData),sprintf(['explains: \n ' num2str(round(PCAE(1),2)) '%% \n']),...
'Interpreter','tex','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',14);

text(2,max(ph2.YData),sprintf(['explains: \n ' num2str(round(PCAE(2),2)) '%% \n']),...
'Interpreter','tex','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',14);

text(3,max(ph3.YData),sprintf(['explains: \n ' num2str(round(PCAE(3),2)) '%% \n']),...
'Interpreter','tex','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',14);

text(4,max(ph4.YData),sprintf(['explains: \n ' num2str(round(PCAE(4),2)) '%% \n']),...
'Interpreter','tex','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',14);

text(5,max(ph5.YData),sprintf(['explains: \n ' num2str(round(PCAE(5),2)) '%% \n']),...
'Interpreter','tex','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',14);




xlabel('\bf PCA Components 1-5 ');ylabel('\bf Eigenvector Coefficients');
cb=colorbar('Ticks',unique(Tx),'TickLabels',(unique(Tx)),'Direction','reverse');
cb.Label.String = Ts; cb.Label.FontSize = 18; cb.Label.Rotation = -90;
cb.Label.VerticalAlignment = 'baseline';







axes(hax2); colormap(hax2, (prism(max(Tn))) );
% axes(hax2); colormap(hax2, flipud(jet(max(Tn))) );

ph1 = scatter( zeros(size(PCAC(:,1)))+1  , PCAC(:,1),50, Tx,'o','MarkerEdgeAlpha',.02);
hold on
ph2 = scatter( zeros(size(PCAC(:,2)))+2  , PCAC(:,2),50, Tx,'o','MarkerEdgeAlpha',.02);
hold on
ph3 = scatter( zeros(size(PCAC(:,3)))+3  , PCAC(:,3),50, Tx,'o','MarkerEdgeAlpha',.02);
hold on
ph4 = scatter( zeros(size(PCAC(:,4)))+4  , PCAC(:,4),50, Tx,'o','MarkerEdgeAlpha',.02);
hold on
ph5 = scatter( zeros(size(PCAC(:,5)))+5  , PCAC(:,5),50, Tx,'o','MarkerEdgeAlpha',.02);

hax2.XLim = [0 6];




clc; clearvars -except ADSP LOCI CASE CTRL USNP PHEN PHE PHECASE PHECTRL...
AMX AMXCASE AMXCTRL AMXUSNP ADNN PCAMX PCAC PCAS PCAE PCAMX_REDUX






%% PLOT PCA --- COHORT ------------------------------------------
clc; close all;
fh1=figure('Units','normalized','OuterPosition',[.02 .06 .97 .8],'Color','w');
hax1 = axes('Position',[.05 .07 .40 .85],'Color','none');
hax2 = axes('Position',[.52 .07 .40 .85],'Color','none');

Ts = 'Cohort';

Tx = PHE.COHORTNUM;

%---Make Tx range from 1:numel(Tn)
Tu = unique(Tx);
[Tx,~] = find((Tx==Tu')');
Tn = unique(Tx);
%-----------------------------

axes(hax1); colormap(hax1, flipud(jet(max(Tn))) );
ph1 = scatter(PCAC(:,1), PCAC(:,2),500, Tx,'.');
xlabel('\bf PCA1');ylabel('\bf PCA2');

cb=colorbar('Ticks',Tn,'TickLabels',Tu,'Direction','reverse');
cb.Label.String = Ts; cb.Label.FontSize = 14; cb.Label.Rotation = -90;
cb.Label.VerticalAlignment = 'baseline';


axes(hax2); colormap(hax2, flipud(jet(max(Tn))) );
ph2 = scatter3(PCAC(:,1), PCAC(:,2), PCAC(:,3),500, Tx,'.');
xlabel('\bf PCA1');ylabel('\bf PCA2');zlabel('\bf PCA3');
view([65,20])


clc; clearvars -except ADSP LOCI CASE CTRL USNP PHEN PHE PHECASE PHECTRL...
AMX AMXCASE AMXCTRL AMXUSNP ADNN PCAMX PCAC PCAS PCAE PCAMX_REDUX







%% PLOT PCA --- CASE-vs-CTRL ------------------------------------------
clc; close all;
fh1=figure('Units','normalized','OuterPosition',[.02 .06 .97 .8],'Color','w');
hax1 = axes('Position',[.05 .07 .40 .85],'Color','none');
hax2 = axes('Position',[.52 .07 .40 .85],'Color','none');


Tx = PHE.AD;
Tn = unique(Tx);
Ts = 'CASE-vs-CTRL';


axes(hax1); colormap(hax1, (lines(max(Tn)+1)) );
ph1 = scatter(PCAC(:,1), PCAC(:,2),500, Tx,'.');
xlabel('\bf PCA1');ylabel('\bf PCA2');

cb=colorbar('Ticks',Tn,'TickLabels',Tn,'Direction','reverse');
cb.Label.String = Ts; cb.Label.FontSize = 14; cb.Label.Rotation = -90;
cb.Label.VerticalAlignment = 'baseline';


axes(hax2); colormap(hax2, (lines(max(Tn)+1)) );
ph2 = scatter3(PCAC(:,1), PCAC(:,2), PCAC(:,3),500, Tx,'.');
xlabel('\bf PCA1');ylabel('\bf PCA2');zlabel('\bf PCA3');
view([65,20])



clc; clearvars -except ADSP LOCI CASE CTRL USNP PHEN PHE PHECASE PHECTRL...
AMX AMXCASE AMXCTRL AMXUSNP ADNN PCAMX PCAC PCAS PCAE PCAMX_REDUX









%% PLOT PCA --- APOE-STATUS ------------------------------------------
clc; close all;
fh1=figure('Units','normalized','OuterPosition',[.02 .06 .97 .8],'Color','w');
hax1 = axes('Position',[.05 .07 .40 .85],'Color','none');
hax2 = axes('Position',[.52 .07 .40 .85],'Color','none');

Ts = 'APOE-STATUS';


Tx = PHE.APOE;

%---Make Tx range from 1:numel(Tn)
Tu = unique(Tx);
[Tx,~] = find((Tx==Tu')');
Tn = unique(Tx);
%-----------------------------



axes(hax1); colormap(hax1, (lines(numel(Tn))) );
% axes(hax1); colormap(hax1, flipud(jet(max(Tn))) );
ph1 = scatter(PCAC(:,1), PCAC(:,2),500, Tx,'.');
xlabel('\bf PCA1');ylabel('\bf PCA2');

cb=colorbar('Ticks',Tn,'TickLabels',Tu,'Direction','reverse');
cb.Label.String = Ts; cb.Label.FontSize = 14; cb.Label.Rotation = -90;
cb.Label.VerticalAlignment = 'baseline';


axes(hax2); colormap(hax2, (lines(numel(Tn))) );
% axes(hax2); colormap(hax2, flipud(jet(max(Tn))) );
ph2 = scatter3(PCAC(:,1), PCAC(:,2), PCAC(:,3),500, Tx,'.');
xlabel('\bf PCA1');ylabel('\bf PCA2');zlabel('\bf PCA3');
view([65,20])


clc; clearvars -except ADSP LOCI CASE CTRL USNP PHEN PHE PHECASE PHECTRL...
AMX AMXCASE AMXCTRL AMXUSNP ADNN PCAMX PCAC PCAS PCAE PCAMX_REDUX







%% PLOT PCA --- AGE ------------------------------------------
clc; close all;
fh1=figure('Units','normalized','OuterPosition',[.02 .06 .97 .8],'Color','w');
hax1 = axes('Position',[.05 .07 .40 .85],'Color','none');
hax2 = axes('Position',[.52 .07 .40 .85],'Color','none');

Ts = 'AGE  (bin)';

AGE = round(PHE.AGE);
[Y,E] = discretize(AGE,8);
for nn = 1:(numel(E)-1)
AGE(Y==nn) = round((E(nn)+E(nn+1))/2);
end

% ONLY INCLUDE PEOPLE OVER 60 YEARS OLD
Tx = AGE(AGE>60);
Tn = unique(Tx);
PCAage = PCAC(AGE>60,:);

% INCLUDE PEOPLE OF ANY AGE
% Tx = AGE;
% Tn = unique(Tx);
% PCAage = PCAC;


axes(hax1); colormap(hax1, (jet(numel(Tn))) );
% axes(hax1); colormap(hax1, flipud(jet(max(Tn))) );
ph1 = scatter(PCAage(:,1), PCAage(:,2),500, Tx,'.');
xlabel('\bf PCA1');ylabel('\bf PCA2');

cb=colorbar('Ticks',Tn,'TickLabels',round(Tn),'Direction','reverse');
cb.Label.String = Ts; cb.Label.FontSize = 14; cb.Label.Rotation = -90;
cb.Label.VerticalAlignment = 'baseline';


axes(hax2); colormap(hax2, (jet(numel(Tn))) );
% axes(hax2); colormap(hax2, flipud(jet(max(Tn))) );
ph2 = scatter3(PCAage(:,1), PCAage(:,2), PCAage(:,3),500, Tx,'.');
xlabel('\bf PCA1');ylabel('\bf PCA2');zlabel('\bf PCA3');
view([65,20])


clc; clearvars -except ADSP LOCI CASE CTRL USNP PHEN PHE PHECASE PHECTRL...
AMX AMXCASE AMXCTRL AMXUSNP ADNN PCAMX PCAC PCAS PCAE PCAMX_REDUX





%% PLOT PCA --- RACE ------------------------------------------
clc; close all;
fh1=figure('Units','normalized','OuterPosition',[.02 .06 .97 .8],'Color','w');
hax1 = axes('Position',[.05 .07 .40 .85],'Color','none');
hax2 = axes('Position',[.52 .07 .40 .85],'Color','none');

Ts = 'RACE';



Tx = PHE.RACE;

Tx(isnan(Tx)) = [];

%---Make Tx range from 1:numel(Tn)
Tu = unique(Tx);
[Tx,~] = find((Tx==Tu')');
Tn = unique(Tx);
%-----------------------------


axes(hax1); colormap(hax1, (lines(max(Tn))) );
ph1 = scatter(PCAC(:,1), PCAC(:,2),500, Tx,'.');
xlabel('\bf PCA1');ylabel('\bf PCA2');

cb=colorbar('Ticks',Tn,'TickLabels',Tn,'Direction','reverse');
cb.Label.String = Ts; cb.Label.FontSize = 14; cb.Label.Rotation = -90;
cb.Label.VerticalAlignment = 'baseline';


axes(hax2); colormap(hax2, (lines(max(Tn))) );
ph2 = scatter3(PCAC(:,1), PCAC(:,2), PCAC(:,3),500, Tx,'.');
xlabel('\bf PCA1');ylabel('\bf PCA2');zlabel('\bf PCA3');
view([65,20])



clc; clearvars -except ADSP LOCI CASE CTRL USNP PHEN PHE PHECASE PHECTRL...
AMX AMXCASE AMXCTRL AMXUSNP ADNN PCAMX PCAC PCAS PCAE PCAMX_REDUX



%% PLOT PCA --- MISC ------------------------------------------
clc; close all;
fh1=figure('Units','normalized','OuterPosition',[.02 .06 .97 .8],'Color','w');
hax1 = axes('Position',[.05 .07 .40 .85],'Color','none');
hax2 = axes('Position',[.52 .07 .40 .85],'Color','none');


Ts = 'ETHNICITY';
Tx = PHE.ETHNIC;
Tu = unique(Tx);
Tu(isnan(Tu)) = [];
[Tx,~] = find((Tx==Tu')');
Tn = unique(Tx);
%-----------------------------



% Ts = 'CONSORTIUM';
% Tx = PHE.CONSORTIUM;
% Txx = ones(size(Tx,1),1);
% Txx(strcmp(Tx,'C')) = 2;
% Tx = Txx;
% Tn = unique(Tx);
% Tn(isnan(Tn)) = [];
%-----------------------------


% Ts = 'BRAAK';
% Tx = PHE.BRAAK;
% Tn = unique(Tx);
% Tn(isnan(Tn)) = [];
%-----------------------------


axes(hax1); colormap(hax1, (lines(max(Tn)+1)) );
ph1 = scatter(PCAC(:,1), PCAC(:,2),500, Tx,'.');
xlabel('\bf PCA1');ylabel('\bf PCA2');

cb=colorbar('Ticks',Tn,'TickLabels',Tn,'Direction','reverse');
cb.Label.String = Ts; cb.Label.FontSize = 14; cb.Label.Rotation = -90;
cb.Label.VerticalAlignment = 'baseline';


axes(hax2); colormap(hax2, (lines(max(Tn)+1)) );
ph2 = scatter3(PCAC(:,1), PCAC(:,2), PCAC(:,3),500, Tx,'.');
xlabel('\bf PCA1');ylabel('\bf PCA2');zlabel('\bf PCA3');
view([65,20])



clc; clearvars -except ADSP LOCI CASE CTRL USNP PHEN PHE PHECASE PHECTRL...
AMX AMXCASE AMXCTRL AMXUSNP ADNN PCAMX PCAC PCAS PCAE PCAMX_REDUX

