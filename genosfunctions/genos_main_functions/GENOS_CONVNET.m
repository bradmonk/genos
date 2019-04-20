% #######################################################################
%%         BUILD NEURAL NET TRAINING AND TESTING MATRICES
% #######################################################################
clc; clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
AMX AMXCASE AMXCTRL AMXUSNP AMXTRCASE AMXTRCTRL AMXTECASE AMXTECTRL


AMXTRCACO = [AMXTRCASE; AMXTRCTRL];
AMXTECACO = [AMXTECASE; AMXTECTRL];


doQ = 0;
[TRNN, ~, ~] = nnmake(AMX,AMXCASE,AMXCTRL,AMXUSNP,AMXTRCACO,doQ);
disp(TRNN(1:20,1:10))

[TENN, ~, ~] = nnmake(AMX,AMXCASE,AMXCTRL,AMXUSNP,AMXTECACO,doQ);
disp(TENN(1:20,1:10))


% SCRAMBLE TRAINING MATRIX
TRAIN  = TRNN(2:end,1:end);
N      = size(TRAIN,1);            % Total number of people
xi     = randperm(N)';             % Get N random ints in range 1:N
TRAIN  = TRAIN(xi,:);              % Scramble NN matrix

SRR.train  = TRAIN(:,1);
TRLabs     = TRAIN(:,2);
SRR.TRLabs = TRLabs;

TRAINL = zeros(size(TRLabs,1),2);  % Make 2-D class label matrix
TRAINL(:,1) = TRLabs==1;           % In col-1 set CASEs index 1's
TRAINL(:,2) = TRLabs==0;           % In col-2 set CTRLs index 1's

TRAINX = TRAIN(:,4:end);           % Strip case/ctrl info

TRAINMX  = TRAINX';
TRAINLAB = TRAINL';




% SCRAMBLE TESTING MATRIX
TEST = TENN(2:end,1:end);
N      = size(TEST,1);             % Total number of people
xi     = randperm(N)';             % Get N random ints in range 1:N
TEST = TEST(xi,:);                 % Scramble NN matrix

SRR.test   = TEST(:,1);
TELabs     = TEST(:,2);    
SRR.TELabs = TELabs;

TESTL = zeros(size(TELabs,1),2);   % Make 2-D class label matrix
TESTL(:,1) = TELabs==1;            % In col-1 set CASEs index 1's
TESTL(:,2) = TELabs==0;            % In col-2 set CTRLs index 1's

TESTX  = TEST(:,4:end);            % Strip case/ctrl info

TESTMX  = TESTX';
TESTLAB = TESTL';




fprintf('  %-4.0f  CASE labels in training dataset \n',sum(TRAINLAB(1,:)))
fprintf('  %-4.0f  CTRL labels in training dataset \n',sum(TRAINLAB(2,:)))
fprintf('  %-4.0f  CASE labels in testing  dataset \n',sum(TESTLAB(1,:)))
fprintf('  %-4.0f  CTRL labels in testing  dataset \n',sum(TESTLAB(2,:)))





clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
AMX AMXCASE AMXCTRL AMXUSNP AMXTRCASE AMXTRCTRL AMXTECASE AMXTECTRL...
SRR TRNN TRAIN TEST TRAINMX TRAINLAB TESTMX TESTLAB


% #######################################################################
%%         MACHINE LEARNING USING ARTIFICIAL NEURAL NETWORKS
% #######################################################################
% TRAIN THE NEURAL NETWORK
clc; clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
AMX AMXCASE AMXCTRL AMXUSNP AMXTRCASE AMXTRCTRL AMXTECASE AMXTECTRL...
SRR TRNN TRAIN TEST TRAINMX TRAINLAB TESTMX TESTLAB


net = patternnet([200 100 50]);

% net.divideParam.trainRatio = 70/100;
% net.divideParam.valRatio = 15/100;
% net.divideParam.testRatio = 15/100;


net = configure(net,TRAINMX,TRAINLAB);

[net,tr] = train(net,TRAINMX,TRAINLAB);

HC = [.80 .20];
[NNPerf] = netstats(net,TRAINMX,TRAINLAB,TESTMX,TESTLAB,.05,HC(1),HC(2),1);


















% #######################################################################
%%                  TRAIN CONVOLUTIONAL NEURAL NET
% #######################################################################
clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
AMX AMXCASE AMXCTRL AMXUSNP AMXTRCASE AMXTRCTRL AMXTECASE AMXTECTRL...
SRR TRNN TRAIN TEST TRAINMX TRAINLAB TESTMX TESTLAB
clc; close all; pause(.5)




ROWPIX = 10;
COLPIX = size(TRAINMX,1) / ROWPIX;

% PREP TRAINING DATA FOR CONVOLUTIONAL NEURAL NET
% Before reshaping, make sure input matrix is VARIANTS-BY-PEOPLE!!!

TRAINMXS  = reshape(TRAINMX , ROWPIX, COLPIX, 1,[]);

TRAINLABS = categorical(TRAINLAB(1,:)');






% (EITHER...) USE HOLD-OUT TEST SET AS VALIDATION SET
TESTMXS   = reshape(TESTMX  , ROWPIX, COLPIX,1,[]);
TESTLABS  = categorical(TESTLAB(1,:)');
VALMATRIX = TESTMXS;
VALLABELS = TESTLABS;


% (OR...) USE RANDOM TRAINING DATA SUBSET AS VALIDATION SET
% idx = randperm(size(TRAINMXS,4),500);
% VALMATRIX = TRAINMXS(:,:,:,idx);
% TRAINMXS(:,:,:,idx) = [];
% VALLABELS = TRAINLABS(idx);
% TRAINLABS(idx) = [];


% COMBINE VALIDATION AND HOLD-OUT TEST GROUPS
% size(TRAINMXS); size(TESTMXS); size(VALMATRIX)
% size(TRAINLABS); size(TESTLABS); size(VALLABELS)
% VALMATRIX = cat(4, VALMATRIX, TESTMXS);
% VALLABELS = cat(1, VALLABELS, TESTLABS);




% NOTES ON BUILDING NN LAYERS AND LAYER PARAMETER
%-------------------------------------------------
%{

%############    convolution2dLayer    ###############
% 
% layer = convolution2dLayer(filterSize,numFilters)
% layer = convolution2dLayer(11,96,'Stride',4,'Padding',1) % example
%
% A 2-D convolutional layer applies sliding filters to the input. 
% The layer convolves the input by moving the filters along the input 
% vertically and horizontally and computing the dot product of the 
% weights and the input, and then adding a bias term.
%
% The example above creates a 2-D convolutional layer that...
% - convolves 96 filters of size [11 11]
% - sliding the filters along rows/cols with a 'stride' of [4 4]
% - adding 1 pixel of 'zero-padding' along all edges of the input image/matrix



%############    EXAMPLE    ###############

layers = [
    imageInputLayer([20 20 1])
    convolution2dLayer(2,10,'Padding',1)
    batchNormalizationLayer
    reluLayer   
    maxPooling2dLayer(2,'Stride',2)
    convolution2dLayer(2,32,'Padding',1)
    batchNormalizationLayer
    reluLayer   
    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer
];


options = trainingOptions('sgdm',...
'MaxEpochs',20, ...
'ValidationData',{TESTMXS,TESTLABS},...
'ValidationFrequency',20,...
'Verbose',false,...
'Plots','training-progress');




options = trainingOptions('sgdm',...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropFactor',0.2,...
    'LearnRateDropPeriod',5,...
    'MaxEpochs',20,...
    'MiniBatchSize',64,...
    'ValidationData',{TESTMXS,TESTLABS},...
    'ValidationFrequency',20,...
    'Verbose',false,...
    'Plots','training-progress')


options = 

  TrainingOptionsSGDM with properties:

                     Momentum: 0.9000
             InitialLearnRate: 0.0100
    LearnRateScheduleSettings: [1x1 struct]
             L2Regularization: 1.0000e-04
                    MaxEpochs: 20
                MiniBatchSize: 64
                      Verbose: 1
             VerboseFrequency: 50
               ValidationData: []
          ValidationFrequency: 50
           ValidationPatience: 5
                      Shuffle: 'once'
               CheckpointPath: ''
         ExecutionEnvironment: 'auto'
                   WorkerLoad: []
                    OutputFcn: []
                        Plots: 'training-progress'
               SequenceLength: 'longest'
         SequencePaddingValue: 0


The 'ValidationPatience' value is the number of times that the loss on the
validation set can be larger than or equal to the previously smallest loss
before network training stops. To turn off automatic validation stopping,
specify Inf as the 'ValidationPatience' value.



net = trainNetwork(TRAINMXS,TRAINLABS,layers,options);
%}
%-------------------------------------------------



% PREPARE CONVNET LAYERS
%-------------------------
layers = [
    imageInputLayer([ROWPIX COLPIX 1])

    convolution2dLayer(1,10,'Padding',1)
    batchNormalizationLayer
    reluLayer

    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,16,'Padding',1)
    batchNormalizationLayer
    reluLayer

    averagePooling2dLayer(2,'Stride',2)

    convolution2dLayer(3,16,'Padding',1)
    batchNormalizationLayer
    reluLayer

    dropoutLayer(.1)

    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer
];




% SET CONVNET OPTIONS
%-------------------------
numValidationsPerEpoch = 2;
miniBatchSize = 100;
validationFrequency = floor(size(TRAINMXS,4)/miniBatchSize/numValidationsPerEpoch);

options = trainingOptions(... 
    'sgdm', ...
    'ValidationData',       {VALMATRIX,VALLABELS}, ...
    'InitialLearnRate',     0.01,...
    'MaxEpochs',            100, ...
    'ValidationPatience',   5, ...
    'ValidationFrequency',  validationFrequency, ...
    'MiniBatchSize',        miniBatchSize, ...
    'L2Regularization',     0.0005, ...
    'Momentum',             0.2, ...
    'SequencePaddingValue', 0, ...
    ...
    'Verbose',true, ...
    'Shuffle', 'every-epoch', ...
    'Plots','training-progress' ...
);
%'OutputFcn',@(info)stopNetTraining(info,3) ...




% RUN CONVNET TRAINING
%-------------------------
net = trainNetwork(TRAINMXS,TRAINLABS,layers,options);

disp(TRAINMXS(1:10,1:10,1,1))
disp(TRAINLABS(1:10,1))




%% PLOT CONVNET OUTPUT

figure
[cmat,classNames] = confusionmat(testLabels,predictedLabels);
h = heatmap(classNames,classNames,cmat);
xlabel('Predicted Class');
ylabel('True Class');
title('Confusion Matrix');



figure
idx = randperm(size(testImages,4),9);
for i = 1:numel(idx)
    subplot(3,3,i)
    imshow(testImages(:,:,:,idx(i)));
    prob = num2str(100*max(probs(idx(i),:)),3);
    predClass = char(predictedLabels(idx(i)));
    label = [predClass,', ',prob,'%'];
    title(label)
end
