function [] = netplot(PCTCORR,NNP)

% keyboard

%%

[ui,uj,uk] = unique(NNP.BackpropFcn);


NET1 = PCTCORR(:,:,uk==1);
NET2 = PCTCORR(:,:,uk==4);
NET3 = PCTCORR(:,:,uk==3);
NET4 = PCTCORR(:,:,uk==2);


N = sum(uk==1);


%% PLOT MODEL PERFORMANCE
AXMAX = 100;
AXMIN = 50;

%------------------------------------
% PLOT HISTOGRAMS AND ROC CURVE

close all;
fh1=figure('Units','normalized','OuterPosition',[.10 .05 .85 .93],'Color','w','MenuBar','none');
h1 = axes('Position',[.05 .54 .4 .4],'Color','none'); hold on;
h2 = axes('Position',[.53 .54 .4 .4],'Color','none'); hold on;
h3 = axes('Position',[.05 .05 .4 .4],'Color','none'); hold on;
h4 = axes('Position',[.53 .05 .4 .4],'Color','none'); hold on;


%---------- GRADIENT DESCENT ----------------
Nmu = nanmean(NET1.* 100,3);
Nsd =  nanstd(NET1.* 100,[],3) ./ sqrt(N);
if min(Nmu(:,3)) < 50; AXMIN = min(Nmu(:,3)); end

axes(h1)
ph1 = scatter( Nmu(:,2) , Nmu(:,3), 2000 , Nmu(:,1), '.');
ph1.MarkerEdgeColor = [0 0 0]; 
hold on;
eh1 = errorbar(Nmu(:,2) , Nmu(:,3) , Nsd(:,2), Nsd(:,2), Nsd(:,3), Nsd(:,3));
eh1.LineStyle = 'none';
eh1.LineWidth = 2;
hold on
ph1 = scatter( Nmu(:,2) , Nmu(:,3), 1200 , Nmu(:,1), '.');
h1.YLim = [AXMIN AXMAX];
h1.XLim = [0 100]; 
h1.XLabel.String = 'Pct. Population';
h1.YLabel.String = 'Pct. Correct';
cb = colorbar; cb.Label.String = 'Classifier Confidence Threshold';
cb.Label.Rotation = 270;
cb.Label.VerticalAlignment = 'bottom';
cb.Label.FontSize = 14;
colormap(h1,parula)
title('Gradient descent')





%---------- SCALED CONJUGATE GRADIENT ----------------
Nmu = nanmean(NET2.* 100,3);
Nsd =  nanstd(NET2.* 100,[],3) ./ sqrt(N);
% Nsd =  std(NET2.* 100,[],3);

if min(Nmu(:,3)) < 50; AXMIN = min(Nmu(:,3)); end

axes(h2)
ph1 = scatter( Nmu(:,2) , Nmu(:,3), 2000 , Nmu(:,1), '.');
ph1.MarkerEdgeColor = [0 0 0]; 
hold on;
eh1 = errorbar(Nmu(:,2) , Nmu(:,3) , Nsd(:,2), Nsd(:,2), Nsd(:,3), Nsd(:,3));
eh1.LineStyle = 'none';
eh1.LineWidth = 2;
hold on
ph1 = scatter( Nmu(:,2) , Nmu(:,3), 1200 , Nmu(:,1), '.');
h2.YLim = [AXMIN AXMAX]; 
h2.XLim = [0 100]; 
h2.XLabel.String = 'Pct. Population';
h2.YLabel.String = 'Pct. Correct';
cb = colorbar; cb.Label.String = 'Classifier Confidence Threshold';
cb.Label.Rotation = 270;
cb.Label.VerticalAlignment = 'bottom';
cb.Label.FontSize = 14;
colormap(h2,parula)
title('Scaled conjugate gradient')




%---------- RESILIENT GRADIENT  ----------------
Nmu = nanmean(NET3.* 100,3);
Nsd =  nanstd(NET3.* 100,[],3) ./ sqrt(N);
% Nsd =  std(NET3.* 100,[],3);

if min(Nmu(:,3)) < 50; AXMIN = min(Nmu(:,3)); end

axes(h3)
ph1 = scatter( Nmu(:,2) , Nmu(:,3), 2000 , Nmu(:,1), '.');
ph1.MarkerEdgeColor = [0 0 0]; 
hold on;
eh1 = errorbar(Nmu(:,2) , Nmu(:,3) , Nsd(:,2), Nsd(:,2), Nsd(:,3), Nsd(:,3));
eh1.LineStyle = 'none';
eh1.LineWidth = 2;
hold on
ph1 = scatter( Nmu(:,2) , Nmu(:,3), 1200 , Nmu(:,1), '.');
h3.YLim = [AXMIN AXMAX]; 
h3.XLim = [0 100]; 
h3.XLabel.String = 'Pct. Population';
h3.YLabel.String = 'Pct. Correct';
cb = colorbar; cb.Label.String = 'Classifier Confidence Threshold';
cb.Label.Rotation = 270;
cb.Label.VerticalAlignment = 'bottom';
cb.Label.FontSize = 14;
colormap(h3,parula)
title('Resilient')






%---------- ONE-STEP SECANT  ----------------
Nmu = nanmean(NET4.* 100,3);
Nsd =  nanstd(NET4.* 100,[],3) ./ sqrt(N);
% Nsd =  std(NET4.* 100,[],3);

if min(Nmu(:,3)) < 50; AXMIN = min(Nmu(:,3)); end

axes(h4)
ph1 = scatter( Nmu(:,2) , Nmu(:,3), 2000 , Nmu(:,1), '.');
ph1.MarkerEdgeColor = [0 0 0]; 
hold on;
eh1 = errorbar(Nmu(:,2) , Nmu(:,3) , Nsd(:,2), Nsd(:,2), Nsd(:,3), Nsd(:,3));
eh1.LineStyle = 'none';
eh1.LineWidth = 2;
hold on
ph1 = scatter( Nmu(:,2) , Nmu(:,3), 1200 , Nmu(:,1), '.');
h4.YLim = [AXMIN AXMAX]; 
h4.XLim = [0 100]; 
h4.XLabel.String = 'Pct. Population';
h4.YLabel.String = 'Pct. Correct';
cb = colorbar; cb.Label.String = 'Classifier Confidence Threshold';
cb.Label.Rotation = 270;
cb.Label.VerticalAlignment = 'bottom';
cb.Label.FontSize = 14;
colormap(h4,parula)
title('One-step secant')














%%
end


%% DEEP LEARNING: CONVOLUTIONAL NEURAL NETWORK
%{
clc
clearvars -except ADMX SNPCASE SNPCTRL PHE CASEID CTRLID CASEMX CTRLMX...
COHORTS AMX AMXCASE AMXCTRL ADNN caMX coMX Pfilter...
TRAINX TRAINL TESTX TESTL TRAINMX TRAINLAB TESTMX TESTLAB net

load lettersTrainSet

TRAINMX  = TRAINX';
TRAINLAB = TRAINL';

TESTMX  = TESTX';
TESTLAB = TESTL';


[TTACTUAL,~,~] = find(TRAINLAB);
zTRAINLAB = categorical(TTACTUAL);


zTRAINMX = TRAINMX;
sz = size(zTRAINMX,1);
re = sz-1800;
rp = randperm(sz,re);
zTRAINMX(rp,:) = [];
zTRAINMX = reshape(zTRAINMX,50,36,1,size(zTRAINMX,2));



clc
size(zTRAINMX)
size(zTRAINLAB)
size(XTrain)
size(TTrain)


layers = [imageInputLayer([50 36 1]);
          convolution2dLayer(3,3);
          reluLayer();
          maxPooling2dLayer(2,'Stride',2);
          fullyConnectedLayer(size(TRAINLAB,1));
          softmaxLayer();
          classificationLayer()];

% convolution2dLayer(5,16);
% options = trainingOptions('sgdm');

options = trainingOptions('sgdm',...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropFactor',0.2,...
    'LearnRateDropPeriod',5,...
    'MaxEpochs',20,...
    'MiniBatchSize',64,...
    'Verbose',false,...
    'Plots','training-progress',...
    'OutputFcn',@(info)stopIfAccuracyNotImp(info,3));


net = trainNetwork(zTRAINMX,zTRAINLAB,layers,options);



[TTACTUAL,~,~] = find(TESTLAB);
zTESTLAB = categorical(TTACTUAL);


zTESTMX = TESTMX;
sz = size(zTESTMX,1);
re = sz-1800;
rp = randperm(sz,re);
zTESTMX(rp,:) = [];
zTESTMX = reshape(zTESTMX,50,36,1,size(zTESTMX,2));

YTest = classify(net,zTESTMX);

accuracy = sum(YTest == zTESTLAB)/numel(zTESTLAB);

disp(accuracy)



%% DEEP LEARNING: LSTM NETWORK
%}