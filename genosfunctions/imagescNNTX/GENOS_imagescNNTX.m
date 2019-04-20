%################################################################
%%  ASSESS VARIANT PROFILES OF HIGH-CONFIDENCE POPULATION
%################################################################

close all; clear; clc; rng('shuffle');
try
genosdir = fileparts(which('GENOS_imagescNNTX.m'));
cd(genosdir); addpath(genosdir)
which('NNTX.mat'); load('NNTX.mat')
catch EX
disp('either this script or the dataset is not in your active PATH')
end




% All of these items will load into your workspace:
clearvars -except LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
AMX AMXCASE AMXCTRL AMXUSNP AMXTRCASE AMXTRCTRL AMXTECASE AMXTECTRL...
TRNN TENN TRAINMX TRAINLAB TESTMX TESTLAB TRAINX TESTX...
NN net NNPerf HC


% Note that 'NNTX' is not one of the variables that loads into 
% your workspace. Instead, we use the variable TRAINMX to make
% it clear that we are using the NN training data to understand
% what variants NN considers useful for making high-confidence
% predictions. Also note that the script below can simply be
% appended to the bottom of the main GENOS_NEURALNET.m script
% if for some reason you want to start with freshly generated
% training and testing groups.




%% TRAIN THE NEURAL NET (optional)
% This is optional because the mat data you loaded has a well-trained
% neural net included. Also, note that netstats() is a custom function
% so make sure you have the latest version from the genos google drive
% folder: genos/code/genosfunctions
% link: https://goo.gl/Ht8s3Q


NN = patternnet([500 40 30]);

NN = configure(NN,TRAINMX,TRAINLAB);

net = train(NN,TRAINMX,TRAINLAB);


HC = [.80 .20]; % high confidence bounds

[NNPerf] = netstats(net,TRAINMX,TRAINLAB,TESTMX,TESTLAB,.05,HC(1),HC(2),1);









%% PLOT IMAGESC OF A HIGH-CONFIDENCE-ONLY-NN-MATRIX

ACTIVATION = net(TRAINMX);
GUESSES    = round(ACTIVATION);
CLASSES    = vec2ind(ACTIVATION);


CORRECT = sum(TRAINLAB == GUESSES)>0;

HICONF  = sum(ACTIVATION>HC(1) | ACTIVATION<HC(2))>0;

HICONFCASE  = ACTIVATION(1,:)>HC(1);
HICONFCTRL  = ACTIVATION(2,:)>HC(1);

LOCONF  = sum((ACTIVATION<.6)&(ACTIVATION>.4))>0;


HICASE = TRAINMX(:,HICONFCASE);
HICTRL = TRAINMX(:,HICONFCTRL);
HICASE = HICASE(:,1:200);
HICTRL = HICTRL(:,1:200);
% HICASE = mean(HICASE,2);
% HICTRL = mean(HICTRL,2);


%%---------

HIC = HICASE';
close all
fh1=figure('Units','normalized','OuterPosition',[.01 .05 .95 .6],'Color','w');
t = uitable(fh1,'Data',[cellstr(AMX.GENE'); cellstr(string(HIC))],...
'Units','normalized','Position',[.01 .01 .98 .95]);
t.ColumnEditable = true;



%%---------

close all
fh1 = figure('Units','normalized','OuterPosition',[.01 .05 .95 .92],'Color','w','MenuBar','none');
ax1 = axes('Position',[.06 .05 .43 .90],'Color','none');
ax2 = axes('Position',[.55 .05 .43 .90],'Color','none');


axes(ax1)
imagesc(HICASE)
title('High Confidence Cases')
ax1.YTick = 1:size(HICASE,1);
YLabs = cellstr([repmat('\fontsize{9} ',length(AMX.GENE),1), char(AMX.GENE)])
ax1.YTickLabel = cellstr(YLabs);

axes(ax2)
imagesc(HICTRL)
title('High Confidence Controls')
ax2.YTick = 1:size(HICTRL,1);
ax2.YTickLabel = cellstr(YLabs);


map =  [.10  .10  .60;
        .20  .20  .80;
        .95  .40  .20;
        .99  .99  .60];
colormap(map)



%% PLOT IMAGESC OF MEAN(HIGH-CONFIDENCE-ONLY-NN-MATRIX)

HICASE = TRAINMX(:,HICONFCASE);
HICTRL = TRAINMX(:,HICONFCTRL);

HICASE = mean(HICASE,2);
HICTRL = mean(HICTRL,2);



close all
fh1 = figure('Units','normalized','OuterPosition',[.01 .05 .95 .92],'Color','w','MenuBar','none');
ax1 = axes('Position',[.06 .05 .25 .90],'Color','none');
ax2 = axes('Position',[.33 .05 .27 .90],'Color','none');
ax3 = axes('Position',[.66 .05 .25 .90],'Color','none');


axes(ax1)
imagesc(HICASE)
title('Average High Confidence Case')
ax1.YTick = 1:size(HICASE,1);
YLabs = cellstr([repmat('\fontsize{9} ',length(AMX.GENE),1), char(AMX.GENE)])
ax1.YTickLabel = cellstr(YLabs);
ax1.XTickLabel = [];

axes(ax2)
imagesc(HICTRL)
title('Average High Confidence Control')
ax2.YTick = 1:size(HICTRL,1);
ax2.YTickLabel = [];
ax2.XTickLabel = [];

C = jet();
colormap(C(1:end-6,:))
colorbar


axes(ax3)
imagesc( abs(HICASE-HICTRL))
title('abs |CASE mean - CTRL mean|')
ax3.YTick = 1:size(HICTRL,1);
ax3.YTickLabel = cellstr(YLabs);
ax3.XTickLabel = [];


C = jet();
colormap(ax3, C(12:end-4,:))
colorbar





%% CLEAN UP WORKSPACE

clearvars -except LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
AMX AMXCASE AMXCTRL AMXUSNP AMXTRCASE AMXTRCTRL AMXTECASE AMXTECTRL...
TRNN TENN TRAINMX TRAINLAB TESTMX TESTLAB TRAINX TESTX...
NN net NNPerf HC
