function [varargout] = multinetplot(gMX,NV,g,NETZ)
%%

keyboard





%%


%------------------------------------------
% q = fliplr(1:numel(NV));
m=0;
for j = 1:numel(NV)
m=m+1;

    MX = gMX{g};

    vn = NV(j);

    TRX = MX.PVTRX(:,1:vn+9);
    TEX = MX.PVTEX(:,1:vn+9);
    TRL = dummyvar(categorical(TRX(:,2)~=1)); 
    TEL = dummyvar(categorical(TEX(:,2)~=1));


    TX = TRX(:,10:end)';
    TL = TRL';
    HX = TEX(:,10:end)';
    HL = TEL';


    % NNN = (numel(NETZ)/2);
    NNN = numel(NETZ);


    TRLABEL = vec2ind(TL);
    TELABEL = vec2ind(HL);

    TRmu = zeros(NNN,101);
    TEmu = zeros(NNN,101);
    TRpop = TRmu;
    TEpop = TEmu;



    net = NETZ{j};  
    net.inputs{1}.size


    TRvec = net(TX);
    TEvec = net(HX);

    TRCONF = rescale(max(TRvec),0,1);
    TECONF = rescale(max(TEvec),0,1);

    TRGUESS = vec2ind(TRvec);
    TEGUESS = vec2ind(TEvec);


        n=0;
        for k = 0:.01:1
        n=n+1;

        tr = (TRCONF>=k);
        te = (TECONF>=k);


        TRpop(m,n) = sum(tr)./numel(tr);
        TEpop(m,n) = sum(te)./numel(te);

        TRmu(m,n) = mean(TRLABEL(tr) == TRGUESS(tr));
        TEmu(m,n) = mean(TELABEL(te) == TEGUESS(te));

        end

end




TRMEANALL = nanmean(TRmu);
TEMEANALL = nanmean(TEmu);

TRSEMALL =  nanstd(TRmu) ./ sqrt(NNN);
TESEMALL =  nanstd(TEmu) ./ sqrt(NNN);

TRPOPALL = nanmean(TRpop);
TEPOPALL = nanmean(TEpop);






%--- GET TOP VALUES FOR GRAPH

TRmu(isnan(TRmu)) = 0;
TEmu(isnan(TEmu)) = 0;


TOPN = 10;
[TRTOP,i] = topkrows(TRmu,TOPN);
[TETOP,j] = topkrows(TEmu,TOPN);

TRTOP(TRTOP<=.5) = NaN;
TETOP(TETOP<=.5) = NaN;

% TRPOPTOP = TRpop(i,:);
% TEPOPTOP = TEpop(j,:);
[TRPOPTOP,i] = topkrows(TRpop,TOPN);
[TEPOPTOP,j] = topkrows(TEpop,TOPN);


TRMEANBEST = nanmean(TRTOP);
TEMEANBEST = nanmean(TETOP);


TRNTOP = TOPN - sum(isnan(TRTOP));
TENTOP = TOPN - sum(isnan(TETOP));

TRSEMBEST =  nanstd(TRTOP) ./ sqrt(TRNTOP);
TESEMBEST =  nanstd(TETOP) ./ sqrt(TENTOP);

TRPOPBEST = nanmean(TRPOPTOP);
TEPOPBEST = nanmean(TEPOPTOP);




%%
% TESEM = TESEM + .01;

%-----------------------------------------------------------------
fh1=figure('Units','normalized','OuterPosition',[.05 .05 .80 .90],...
'Color','w','MenuBar','none');
ah11 = axes('Position',[.06 .54 .4 .4],'Color','none'); hold on;
ah21 = axes('Position',[.57 .54 .4 .4],'Color','none'); hold on;
ah31 = axes('Position',[.06 .08 .4 .4],'Color','none'); hold on;
ah41 = axes('Position',[.57 .08 .4 .4],'Color','none'); hold on;

%----- TRAINING ALL CLASSIFIERS
axes(ah11)
eh11 = errorbar( (0:.01:1) , TRMEANALL , TRSEMALL);
eh11.LineStyle = 'none'; 
eh11.LineWidth = 1.5; 
eh11.Color = [.2 .2 .8]; hold on;
ph11 = scatter( (0:.01:1) , TRMEANALL , 1500 , '.'); hold on
ph11.MarkerEdgeColor = [0 0 0]; hold on;
ph11 = scatter( (0:.01:1) , TRMEANALL , 1000 , TRPOPALL , '.'); hold on
line([0 1],[.5 .5],'Color','k','LineStyle',':','LineWidth',1)
ah11.YLim = [.3 1];
ah11.XLim = [0 1];
ah11.XLabel.String = 'Classifier Confidence';
ah11.YLabel.String = 'Proportion Correct';
cb11 = colorbar; 
cb11.Label.String = 'Proportion of Population';
cb11.Label.Rotation = 270;
cb11.Label.VerticalAlignment = 'bottom';
cb11.Label.FontSize = 14;
colormap(ah11,parula)
title('Gradient descent')


%----- TESTING ALL CLASSIFIERS
axes(ah21)
eh21 = errorbar( (0:.01:1) , TEMEANALL , TESEMALL);
eh21.LineStyle = 'none'; 
eh21.LineWidth = 1.5; 
eh21.Color = [.2 .2 .8]; hold on;
ph21 = scatter( (0:.01:1) , TEMEANALL , 1500 , '.'); hold on
ph21.MarkerEdgeColor = [0 0 0]; hold on;
ph21 = scatter( (0:.01:1) , TEMEANALL , 1000 , TEPOPALL , '.'); hold on
line([0 1],[.5 .5],'Color','k','LineStyle',':','LineWidth',1)
ah21.YLim = [.3 1];
ah21.XLim = [0 1];
ah21.XLabel.String = 'Classifier Confidence';
ah21.YLabel.String = 'Proportion Correct';
cb21 = colorbar; 
cb21.Label.String = 'Proportion of Population';
cb21.Label.Rotation = 270;
cb21.Label.VerticalAlignment = 'bottom';
cb21.Label.FontSize = 14;
colormap(ah11,parula)
title('Gradient descent')


%----- TRAINING BEST CLASSIFIERS
axes(ah31)
eh31 = errorbar( (0:.01:1) , TRMEANBEST , TRSEMBEST);
eh31.LineStyle = 'none'; 
eh31.LineWidth = 1.5; 
eh31.Color = [.2 .2 .8]; hold on;
ph31 = scatter( (0:.01:1) , TRMEANBEST , 1500 , '.'); hold on
ph31.MarkerEdgeColor = [0 0 0]; hold on;
ph31 = scatter( (0:.01:1) , TRMEANBEST , 1000 , TRPOPBEST , '.'); hold on
line([0 1],[.5 .5],'Color','k','LineStyle',':','LineWidth',1)
ah31.YLim = [.3 1];
ah31.XLim = [0 1];
ah31.XLabel.String = 'Classifier Confidence';
ah31.YLabel.String = 'Proportion Correct';
cb31 = colorbar; 
cb31.Label.String = 'Proportion of Population';
cb31.Label.Rotation = 270;
cb31.Label.VerticalAlignment = 'bottom';
cb31.Label.FontSize = 14;
colormap(ah31,parula)
title('Gradient descent')


%----- TESTING BEST CLASSIFIERS
axes(ah41)
eh41 = errorbar( (0:.01:1) , TEMEANBEST , TESEMBEST);
eh41.LineStyle = 'none'; 
eh41.LineWidth = 1.5; 
eh41.Color = [.2 .2 .8]; hold on;
ph41 = scatter( (0:.01:1) , TEMEANBEST , 1500 , '.'); hold on
ph41.MarkerEdgeColor = [0 0 0]; hold on;
ph41 = scatter( (0:.01:1) , TEMEANBEST , 1000 , TEPOPBEST , '.'); hold on
line([0 1],[.5 .5],'Color','k','LineStyle',':','LineWidth',1)
ah41.YLim = [.3 1];
ah41.XLim = [0 1];
ah41.XLabel.String = 'Classifier Confidence';
ah41.YLabel.String = 'Proportion Correct';
cb41 = colorbar; 
cb41.Label.String = 'Proportion of Population';
cb41.Label.Rotation = 270;
cb41.Label.VerticalAlignment = 'bottom';
cb41.Label.FontSize = 14;
colormap(ah11,parula)
title('Gradient descent')
%-----------------------------------------------------------------



%%

NNRESULTS.TRmu       = TRmu;
NNRESULTS.TEmu       = TEmu;
NNRESULTS.TRpop      = TRpop;
NNRESULTS.TEpop      = TEpop;
NNRESULTS.TRMEAN     = TRMEANALL;
NNRESULTS.TEMEAN     = TEMEANALL;
NNRESULTS.TRSEM      = TRSEMALL;
NNRESULTS.TESEM      = TESEMALL;
NNRESULTS.TRPOP      = TRPOPALL;
NNRESULTS.TEPOP      = TEPOPALL;

NNRESULTS.TRTOP      = TRTOP;
NNRESULTS.TETOP      = TETOP;
NNRESULTS.TRMEANBEST = TRMEANBEST;
NNRESULTS.TEMEANBEST = TEMEANBEST;
NNRESULTS.TRSEMBEST  = TRSEMBEST;
NNRESULTS.TESEMBEST  = TESEMBEST;
NNRESULTS.TRPOPBEST  = TRPOPBEST;
NNRESULTS.TEPOPBEST  = TEPOPBEST;


varargout = {NNRESULTS};

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