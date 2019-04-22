function ObjFcn = makeObjFcn(trainImages,trainLabels,valImages,valLabels)
ObjFcn = @valErrorFun;
    function [valError,cons,fileName] = valErrorFun(optVars)
        
        imageSize = [size(trainImages,1) size(trainImages,2) size(trainImages,3)];
        numClasses = numel(unique(trainLabels));
        initialNumFilters = round(size(trainImages,1)/sqrt(optVars.NetworkDepth));
        layers = [
            imageInputLayer(imageSize)

            % The spatial input and output sizes of these convolutional
            % layers are 32-by-32, and the following max pooling layer
            % reduces this to 16-by-16.
            convBlock(3,initialNumFilters,optVars.NetworkDepth)
            maxPooling2dLayer(2,'Stride',2)

            % The spatial input and output sizes of these convolutional
            % layers are 16-by-16, and the following max pooling layer
            % reduces this to 8-by-8.
            convBlock(3,2*initialNumFilters,optVars.NetworkDepth)
            maxPooling2dLayer(2,'Stride',2)

            % The spatial input and output sizes of these convolutional
            % layers are 8-by-8.
            convBlock(3,4*initialNumFilters,optVars.NetworkDepth)

            % Add the fully connected layer and the final softmax and
            % classification layers.
            fullyConnectedLayer(numClasses)
            softmaxLayer
            classificationLayer];


        miniBatchSize = 100;
        numValidationsPerEpoch = 3;
        validationFrequency = floor(size(trainImages,4)/miniBatchSize/numValidationsPerEpoch);
        options = trainingOptions('sgdm',...
            'InitialLearnRate',optVars.InitialLearnRate,...
            'Momentum',optVars.Momentum,...
            'MaxEpochs',100, ...
            'MiniBatchSize',miniBatchSize,...
            'L2Regularization',optVars.L2Regularization,...
            'Shuffle','every-epoch',...
            'Verbose',true,... % 'OutputFcn',@plotTrainingProgress,...
            'ValidationData',{valImages,valLabels},...
            'ValidationPatience',4,...
            'ValidationFrequency',validationFrequency);


        pixelRange = [-4 4];
        imageAugmenter = imageDataAugmenter(...
            'RandXReflection',true,...
            'RandXTranslation',pixelRange,...
            'RandYTranslation',pixelRange);
        datasource = augmentedImageSource(imageSize,trainImages,trainLabels,...
            'DataAugmentation',imageAugmenter,...
            'OutputSizeMode','randcrop');


        trainedNet = trainNetwork(datasource,layers,options);
        close



        options = trainingOptions('sgdm',...
            'InitialLearnRate',optVars.InitialLearnRate/10,...
            'Momentum',optVars.Momentum,...
            'MaxEpochs',100, ...
            'MiniBatchSize',miniBatchSize,...
            'L2Regularization',optVars.L2Regularization,...
            'Shuffle','every-epoch',...
            'Verbose',false,...
            'OutputFcn',@plotTrainingProgress,...
            'ValidationData',{valImages,valLabels},...
            'ValidationPatience',4,...
            'ValidationFrequency',validationFrequency);
        trainedNet = trainNetwork(datasource,trainedNet.Layers,options);
        close


        predictedLabels = classify(trainedNet,valImages);
        valAccuracy = mean(predictedLabels == valLabels);
        valError = 1 - valAccuracy;



        fileName = num2str(valError,10) + ".mat";
        save(fileName,'trainedNet','valError','options')
        cons = [];


    end
end



function layers = convBlock(filterSize,numFilters,numConvLayers)
layers = [
    convolution2dLayer(filterSize,numFilters,'Padding',(filterSize-1)/2)
    batchNormalizationLayer
    reluLayer];
layers = repmat(layers,numConvLayers,1);
end