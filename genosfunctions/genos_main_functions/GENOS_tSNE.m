% ######################################################################
%%       tSNE : t-Distributed Stochastic Neighbor Embedding
% ######################################################################
clc; close all; clear; rng('shuffle')
cd(fileparts(which('GENOS_tSNE.m')));


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


% [ai,bi] = ismember(PHEN.COHORTNUM,GOODCOH);
[ai,bi] = ismember(PHEN.COHORTNUM,[GOODCOH; MEHHCOH]);

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


AAN = string(LOCI.AAN);
AAN = replace(AAN,"NA","NaN;NaN");
AAN = split(AAN,";");
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


N = 200;


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





%% (OPTIONAL) PRE-PERFORM PCA BEFORE TSNE 

% ss = statset('pca');
% ss.Display = 'iter';
% ss.MaxIter = 100;
% ss.TolFun = 1e4;
% ss.TolX = 1e4;
% ss.UseParallel = true;
% 
% [PCAC,PCAS,~,~,~] = pca(  PCAMX' , 'Options',ss);
% clc; close all; scatter(PCAC(:,1),PCAC(:,2))
%
% % ...,'NumPCAComponents',0,...  means don't use PCA
% tSN = tsne(PCAC(:,1:10),'NumDimensions',2,'Theta',.6,'NumPCAComponents',0);
% 
% clearvars -except ADSP GENB LOCI CASE CTRL USNP PHEN AMX AMXCASE AMXCTRL...
% PHE ADNN PCAMX tSN PCAC PCAS






%######################################################################
%%       tSNE : t-Distributed Stochastic Neighbor Embedding
%######################################################################



tSN = tsne(PCAMX ,'NumDimensions',2,'Theta',.6,'NumPCAComponents',8);



clc; clearvars -except ADSP GENB LOCI CASE CTRL USNP PHEN PHE...
AMX AMXCASE AMXCTRL ADNN PCAMX tSN

disp('done')



%######################################################################
%%       PLOT tSNE RESULTS
%######################################################################




%% PLOT TSNE --- ALZHEIMERS STATUS (CASE/CTRL) --------------------------
close all; 
fh1=figure('Units','normalized','Position',[.05 .05 .70 .84],'Color','w');
ax1=axes('Position',[.05 .02 .9 .9],'Color','none');

ph1 = gscatter(tSN(:,1),tSN(:,2),  PHE.AD, [],'.',15);

title({'\fontsize{16} t-SNE : CASE vs CTRL',' '})
legend(ph1,{'CTRL','CASE'},'FontSize',12,'Box','off','Location','NorthWest');
axis off




clc; clearvars -except ADSP GENB LOCI CASE CTRL USNP PHEN PHE...
AMX AMXCASE AMXCTRL ADNN PCAMX tSN




%% PLOT TSNE --- CONSORTIUM STUDY COHORT (1:24) -------------------------
close all; 
fh1=figure('Units','normalized','Position',[.05 .05 .70 .84],'Color','w');
ax1=axes('Position',[.05 .02 .9 .9],'Color','none');

ph1 = gscatter(tSN(:,1),tSN(:,2),  PHE.COHORT, [],'.',15);


title({'\fontsize{16} t-SNE : STUDY COHORT',' '})
axis off





clc; clearvars -except ADSP GENB LOCI CASE CTRL USNP PHEN PHE...
AMX AMXCASE AMXCTRL ADNN PCAMX tSN





%% PLOT TSNE --- SEX (M/F) ----------------------------------------------
close all; 
fh1=figure('Units','normalized','Position',[.05 .05 .70 .84],'Color','w');
ax1=axes('Position',[.05 .02 .9 .9],'Color','none');

ph1 = gscatter(tSN(:,1),tSN(:,2),  PHE.SEX, [],'.',15);


title({'\fontsize{16} t-SNE : SEX',' '})
legend(ph1,{'Male','Female'},'FontSize',12,'Box','off','Location','NorthWest');
axis off





clc; clearvars -except ADSP GENB LOCI CASE CTRL USNP PHEN PHE...
AMX AMXCASE AMXCTRL ADNN PCAMX tSN





%% PLOT TSNE --- AGE (BINNED AGE) ---------------------------------------
close all; 
fh1=figure('Units','normalized','Position',[.05 .05 .70 .84],'Color','w');
ax1=axes('Position',[.05 .02 .9 .9],'Color','none');

AGE = round(PHE.AGE);
ofAGE = AGE>60;
A = AGE(ofAGE);
[Y,E] = discretize(A,[60 80 90 91]);
% [Y,E] = discretize(A,[60 75 85 90 91]);
for nn = 1:numel(E)
A(Y==nn) = E(nn);
end



ph1 = gscatter(tSN(ofAGE,1),tSN(ofAGE,2),  A, [],'.',15);

title({'\fontsize{16} t-SNE : AGE',' '})
axis off






clc; clearvars -except ADSP GENB LOCI CASE CTRL USNP PHEN PHE...
AMX AMXCASE AMXCTRL ADNN PCAMX tSN






%% PLOT TSNE --- APOE STATUS (22,23,24,33,34,44) ------------------------
close all; 
fh1=figure('Units','normalized','Position',[.05 .05 .70 .84],'Color','w');
ax1=axes('Position',[.05 .02 .9 .9],'Color','none');


ph1 = gscatter(tSN(:,1),tSN(:,2),  PHE.APOE, [],'.',15);


ph1(1).MarkerSize = 35;
ph1(2).MarkerSize = 25;
ph1(2).Color = [.20 .20 .99];
ph1(3).MarkerSize = 35;
ph1(4).Color = [.99 .50 .10];
ph1(5).Color = [.30 .70 .80];
ph1(6).MarkerSize = 25;

title({'\fontsize{16} t-SNE : APOE',' '})
% legend(ph1,{'CTRL','CASE'},'FontSize',12,'Box','off','Location','NorthWest');
axis off






clc; clearvars -except ADSP GENB LOCI CASE CTRL USNP PHEN PHE...
AMX AMXCASE AMXCTRL ADNN PCAMX tSN





%% PLOT TSNE --- RACE * ETHNIC GROUP ------------------------------------------
close all; 
fh1=figure('Units','normalized','Position',[.05 .05 .70 .84],'Color','w');
ax1=axes('Position',[.05 .02 .9 .9],'Color','none');


ph1 = gscatter(tSN(:,1),tSN(:,2),  PHE.RACE .* (PHE.ETHNIC.*2+1), [],'.',15);


title({'\fontsize{16} t-SNE : RACE*ETHINC GROUP',' '})
axis off



clc; clearvars -except ADSP GENB LOCI CASE CTRL USNP PHEN PHE...
AMX AMXCASE AMXCTRL ADNN PCAMX tSN






%% PLOT TSNE --- CONSENT GROUP ------------------------------------------
close all; 
fh1=figure('Units','normalized','Position',[.05 .05 .70 .84],'Color','w');
ax1=axes('Position',[.05 .02 .9 .9],'Color','none');


ph1 = gscatter(tSN(:,1),tSN(:,2),  PHE.Consent, [],'.',15);


title({'\fontsize{16} t-SNE : CONSENT GROUP',' '})
axis off



clc; clearvars -except ADSP GENB LOCI CASE CTRL USNP PHEN PHE...
AMX AMXCASE AMXCTRL ADNN PCAMX tSN


























%######################################################################
%%       tSNE : t-Distributed Stochastic Neighbor Embedding
%######################################################################

% RECOMPUTE TSNE BUT FROM A VARIANT PERSPECTIVE 
% (RATHER THAN A PERSON-PERSPECTIVE; IE COLLAPSE ACROSS PEOPLE)



tSN = tsne(PCAMX' ,'NumDimensions',2,'Theta',.6,'NumPCAComponents',8);


clc; clearvars -except ADSP GENB LOCI CASE CTRL USNP PHEN PHE...
AMX AMXCASE AMXCTRL ADNN PCAMX tSN
disp('done')









%% PLOT TSNE --- CHROMOSOME ------------------------------------------
close all; 
fh1=figure('Units','normalized','Position',[.05 .05 .70 .84],'Color','w');
ax1=axes('Position',[.05 .02 .9 .9],'Color','none');



ph1 = gscatter( tSN(:,1) , tSN(:,2) , AMX.CHR , [],'.',35);

title({'\fontsize{16} t-SNE : CHROMOSOME',' '})
axis off


clc; clearvars -except ADSP GENB LOCI CASE CTRL USNP PHEN PHE...
AMX AMXCASE AMXCTRL ADNN PCAMX tSN



%% PLOT TSNE --- EFFECT ------------------------------------------
close all; 
fh1=figure('Units','normalized','Position',[.05 .05 .70 .84],'Color','w');
ax1=axes('Position',[.05 .02 .9 .9],'Color','none');



ph1 = gscatter( tSN(:,1) , tSN(:,2) , AMX.EFFECT , [],'.',35);

title({'\fontsize{16} t-SNE : CODING EFFECT',' '})
axis off


clc; clearvars -except ADSP GENB LOCI CASE CTRL USNP PHEN PHE...
AMX AMXCASE AMXCTRL ADNN PCAMX tSN




%% PLOT TSNE --- ODDS RATIO > 1 ------------------------------------------
close all; 
fh1=figure('Units','normalized','Position',[.05 .05 .70 .84],'Color','w');
ax1=axes('Position',[.05 .02 .9 .9],'Color','none');



ph1 = gscatter( tSN(:,1) , tSN(:,2) , AMX.FISHOR>1 , [],'.',35);

title({'\fontsize{16} t-SNE : ODDS RATIO > 1',' '})
axis off


clc; clearvars -except ADSP GENB LOCI CASE CTRL USNP PHEN PHE...
AMX AMXCASE AMXCTRL ADNN PCAMX tSN


























%% CLEAN-UP WORKSPACE

clc; clearvars -except ADSP GENB LOCI CASE CTRL USNP PHEN PHE...
AMX AMXCASE AMXCTRL ADNN PCAMX tSN



%% NOTES ON tSNE ADVANCED OPTIONS
%{


%% Use various distance measures to try to obtain a better separation

load fisheriris

rng default % for reproducibility
Y = tsne(meas,'Algorithm','exact','Distance','mahalanobis');
subplot(2,2,1)
gscatter(Y(:,1),Y(:,2),species)
title('Mahalanobis')

rng default % for fair comparison
Y = tsne(meas,'Algorithm','exact','Distance','cosine');
subplot(2,2,2)
gscatter(Y(:,1),Y(:,2),species)
title('Cosine')

rng default % for fair comparison
Y = tsne(meas,'Algorithm','exact','Distance','chebychev');
subplot(2,2,3)
gscatter(Y(:,1),Y(:,2),species)
title('Chebychev')

rng default % for fair comparison
Y = tsne(meas,'Algorithm','exact','Distance','euclidean');
subplot(2,2,4)
gscatter(Y(:,1),Y(:,2),species)
title('Euclidean')



%% Find 2-D and 3-D embeddings of dataset, and compare loss for each embedding

load fisheriris
rng default % for reproducibility
[Y,loss] = tsne(meas,'Algorithm','exact');
rng default % for fair comparison
[Y2,loss2] = tsne(meas,'Algorithm','exact','NumDimensions',3);
fprintf('2-D embedding has loss %g, and 3-D embedding has loss %g.\n',loss,loss2)


gscatter(Y(:,1),Y(:,2),species,eye(3))
title('2-D Embedding')
figure
v = double(categorical(species));
c = full(sparse(1:numel(v),v,ones(size(v)),numel(v),3));
scatter3(Y2(:,1),Y2(:,2),Y2(:,3),15,c,'filled')
title('3-D Embedding')
view(-50,8)




%% tSNE advanced options


%--- ALGORITHM

Y = tsne(x,'Algorithm','barneshut'); % default
Y = tsne(x,'Algorithm','exact');

% The 'exact' algorithm optimizes the Kullback-Leibler divergence of 
% distributions between the original space and the embedded space. The 
% 'barneshut' algorithm performs an approximate optimization that is 
% faster and uses less memory when the number of data rows is large.



%--- DISTANCE

Y = tsne(x,'Distance','euclidean'); % default
Y = tsne(x,'Distance','seuclidean');
Y = tsne(x,'Distance','cityblock');
Y = tsne(x,'Distance','etc');

% Distance metric, specified by one of the following. For definitions 
% of the distance metrics, see pdist.
% 
% 'euclidean'   Euclidean distance.
% 
% 'seuclidean'  Standardized Euclidean distance. Each coordinate 
%               difference between rows in X and the query matrix is 
%               scaled by dividing by the corresponding element of the 
%               standard deviation computed from S = nanstd(X).
% 
% 'cityblock'   City block distance.
% 
% 'chebychev'   Chebychev distance, which is the maximum coordinate 
%               difference.
% 
% 'minkowski'   Minkowski distance with exponent 2. This is the same as 
%               Euclidean distance.
% 
% 'mahalanobis' Mahalanobis distance, computed using the positive 
%               definite covariance matrix nancov(X).
% 
% 'cosine'      One minus the cosine of the included angle between 
%               observations (treated as vectors).
% 
% 'correlation' One minus the sample linear correlation between 
%               observations (treated as sequences of values).
% 
% 'spearman'    One minus the sample Spearman's rank correlation between 
%               observations (treated as sequences of values).
% 
% 'hamming'     Hamming distance, which is the percentage of coordinates 
%               that differ.
% 
% 'jaccard'     One minus the Jaccard coefficient, which is the percentage 
%               of nonzero coordinates that differ.
% 
% 'custom'      A distance function specified using @ (e.g. @distfun).
% 
% In all cases, tsne uses squared pairwise distances to calculate the 
% Gaussian kernel in the joint distribution of X.






%--- EXAGGERATION

Y = tsne(x,'Exaggeration',4); % default
Y = tsne(x,'Exaggeration',10);


% Size of natural clusters in data, specified as a scalar value 1 
% or greater.
% 
% A large exaggeration makes tsne learn larger joint probabilities of 
% Y and creates relatively more space between clusters in Y. tsne uses 
% exaggeration in the first 99 optimization iterations.
% 
% If the value of Kullback-Leibler divergence increases in the early 
% stage of the optimization, try reducing the exaggeration.




%--- NUM DIMENSIONS

Y = tsne(x,'NumDimensions',2); % default
Y = tsne(x,'NumDimensions',3);

% Dimension of the output Y, specified as a positive integer. Generally, 
% set NumDimensions to 2 or 3.





%--- NUM PCA COMPONENTS

Y = tsne(x,'NumPCAComponents',0); % default
Y = tsne(x,'NumPCAComponents',40);


% PCA dimension reduction, specified as a nonnegative integer. Before 
% tsne embeds the high-dimensional data, it first reduces the data's
% dimensionality to NumPCAComponents using the pca function. 
% When NumPCAComponents is 0, tsne does not use PCA.





%--- PERPLEXITY

Y = tsne(x,'Perplexity',30); % default
Y = tsne(x,'Perplexity',40);


% Effective number of local neighbors of each point, specified as a 
% positive scalar. See t-SNE Algorithm.
% 
% Larger perplexity causes tsne to use more points as nearest neighbors. 
% Use a larger value of Perplexity for a large dataset. Typical Perplexity 
% values are from 5 to 50. In the Barnes-Hut algorithm, tsne uses 
% min(3*Perplexity,N-1) as the number of nearest neighbors.





%--- LEARNING RATE

Y = tsne(x,'LearnRate',500); % default
Y = tsne(x,'LearnRate',2000);



% Learning rate for optimization process, specified as a positive scalar. 
% Typically, set values from 100 through 1000.
% 
% When LearnRate is too small, tsne can converge to a poor local minimum. 
% When LearnRate is too large, the optimization can initially have the 
% Kullback-Leibler divergence increase rather than decrease.





%--- OUTPUT OPTIONS

[Y, loss] = tsne(x);

% Y: Embedded points, returned as an n-by-NumDimensions matrix. Each row 
% represents one embedded point. n is the number of rows of data X that 
% do not contain any NaN entries. See Plot Results with NaN Input Data.
%
% loss: Kullback-Leibler divergence between modeled input and output 
% distributions, returned as a nonnegative scalar. 


%}
