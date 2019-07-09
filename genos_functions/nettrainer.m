function [NN] = nettrainer(Y,X,PARAMS)
%======================================================================

% keyboard

% DEAL INPUT PARAMETERS
%======================================================================
NEURONS = PARAMS.NEURONS;
LOOPS = PARAMS.LOOPS;



% INIT PATTERNNET
%======================================================================
net = patternnet(NEURONS,'trainscg','crossentropy');



% PATTERNNET TRAINING PARAMETERS
%======================================================================
net.trainParam.epochs = 1000;             %Maximum number of epochs to train
net.trainParam.show = 25;                 %Epochs between displays (NaN=no displays)
net.trainParam.showCommandLine = false;   %Generate command-line output
net.trainParam.showWindow = false;        %Show training GUI
net.trainParam.goal = 0;                  %Performance goal
net.trainParam.time = inf;                %Maximum time to train in seconds
net.trainParam.min_grad = 1e-6;           %Minimum performance gradient
net.trainParam.max_fail = 50;             %Maximum validation failures
net.trainParam.sigma = 5.0e-5;            %Change in weight for 2nd derivative
net.trainParam.lambda = 5.0e-7;           %Regulator of indefiniteness of the Hessian



% PATTERNNET PERFORMANCE PARAMETERS
%======================================================================
net.performParam.regularization = 0.1;
net.performParam.normalization = 'none';




%======================================================================
% LOOP PARAMETERS
%======================================================================
net_BEST = [];
pctCorrect_BEST = 0;





% RETRAIN PATTERNNET N TIMES AND KEEP BEST
%======================================================================
for nn = 1:LOOPS


    net = train(net,X',Y');
    [ERR,~,~,~] = confusion(Y',net(X'));
    pctCorrect = 1-ERR;

    
    % DETERMINE IF NEW NET > BEST NET
    %-----------------------------------------------
    if (pctCorrect > pctCorrect_BEST)

        net_BEST = net;
        pctCorrect_BEST = pctCorrect;
        net_BEST.userdata.note = sprintf('PCT-CORRECT: %.1f',...
            pctCorrect_BEST*100);
    
    end
    %-----------------------------------------------


    fprintf('\n NN-LOOP: %.0f  \n PCT-CORRECT: %.1f \n\n', nn, pctCorrect_BEST*100);
end

fprintf('\n BEST PERFORMANCE: %0.1f%%\n',  pctCorrect_BEST*100)
%======================================================================


NN = net_BEST;


%======================================================================
end