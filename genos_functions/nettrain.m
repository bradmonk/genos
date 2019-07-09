function [net] = nettrain(TRAINVARS,TESTVARS,PARAMS)
%======================================================================

% keyboard

%======================================================================
%%
%======================================================================

% PVXt, VXt, DVXt, DXt, Yt, YDt
[PVXt, VXt, DVXt, DXt, Yt, YDt] = TRAINVARS{:};
[PVXh, VXh, DVXh, DXh, Yh, YDh] = TESTVARS{:};


doDX = PARAMS.doDX;
LOOPS = PARAMS.LOOPS;
NNEURONS =  PARAMS.NNEURONS;

if doDX==1
    MXt = DXt;
    MXh = DXh;
else
    MXt = VXt;
    MXh = VXh;
end


%======================================================================
% SET NEURAL NET PARAMETERS
%======================================================================
NN = patternnet(NNEURONS,'trainscg','crossentropy');
NN.trainParam.max_fail = 50;
NN.trainParam.showWindow = 0;
NN.performParam.regularization = 0.1;
NN.performParam.normalization = 'none';


net_BEST = [];
net_pctCorrect_BEST = 0;


%======================================================================
% RETRAIN NETS 3X AND KEEP BEST PERFORMER
%-----------------------------------------------
for nn = 1:LOOPS

    %----- TRAIN NET-A
    net = train(NN,MXt',Yt');  % TRAIN NETQ
    [ERR,~,~,~] = confusion(Yh',net(MXh'));
    net_pctCorrect = 1-ERR;

    
    % DETERMINE WHETHER NET-A OR NET-B IS BETTER
    %-----------------------------------------------
    if (net_pctCorrect > net_pctCorrect_BEST)

        net_BEST = net;
        net_pctCorrect_BEST = net_pctCorrect;
        net_BEST.userdata.note = sprintf('PCT-CORRECT: %.1f',...
            net_pctCorrect_BEST*100);
    
    end
    %-----------------------------------------------


    fprintf('\n NN-LOOP: %.0f  \n PCT-CORRECT: %.1f \n\n', nn, net_pctCorrect_BEST*100);
end

fprintf('\n BEST PERFORMANCE: %0.1f%%\n',  net_pctCorrect_BEST*100)
%======================================================================


net = net_BEST;


%======================================================================
end

%======================================================================
% RETRAIN NETS 3X AND KEEP BEST PERFORMER
%-----------------------------------------------
%{
for nn = 1:LOOPS

    %----- TRAIN NET-A
    netA = train(NN,DXt',Yt');  % TRAIN NETQ
    [ERR,~,~,~] = confusion(Yh',netA(DXh'));
    netA_pctCorrect = 1-ERR;


    %----- TRAIN NET-B
    netB = train(NNB,DXt',Yt');  % TRAIN NETQ
    [ERR,~,~,~] = confusion(Yh',netB(DXh'));
    netB_pctCorrect = 1-ERR;

    
    % DETERMINE WHETHER NET-A OR NET-B IS BETTER
    %-----------------------------------------------
    if ((netA_pctCorrect > netB_pctCorrect) &...
       (netA_pctCorrect > netBEST_pctCorrect))

        netBEST = netA;
        netBEST_pctCorrect = netA_pctCorrect;
        netBEST.userdata.note = 'netA is best';
    
    elseif ((netB_pctCorrect > netA_pctCorrect) &...
           (netB_pctCorrect > netBEST_pctCorrect))

        netBEST = netB;
        netBEST_pctCorrect = netB_pctCorrect;
        netBEST.userdata.note = 'netB is best';
    end
    %-----------------------------------------------


    fprintf('\n NN-LOOP: %.0f  \n PCT-CORRECT: %.1f \n\n', nn, netBEST_pctCorrect*100);
end

disp(netBEST.userdata.note)
fprintf('\n It correctly classified: %0.1f%%\n',  netBEST_pctCorrect*100)
%}
%======================================================================

