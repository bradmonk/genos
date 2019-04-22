% #######################################################################
%%         BUILD NEURAL NET TRAINING AND TESTING MATRICES
% #######################################################################
clc; clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
AMX AMXCASE AMXCTRL AMXUSNP AMXTRCASE AMXTRCTRL AMXTECASE AMXTECTRL


AMXTRCACO = [AMXTRCASE; AMXTRCTRL];
AMXTECACO = [AMXTECASE; AMXTECTRL];


doQ = 1;
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


net = patternnet([300 50 50]);

% net.divideParam.trainRatio = 70/100;
% net.divideParam.valRatio = 15/100;
% net.divideParam.testRatio = 15/100;


net = configure(net,TRAINMX,TRAINLAB);

[net,tr] = train(net,TRAINMX,TRAINLAB);

HC = [.80 .20];
[NNPerf] = netstats(net,TRAINMX,TRAINLAB,TESTMX,TESTLAB,.05,HC(1),HC(2),1);


% Test the Network
outputs = net(TESTMX);
errors = gsubtract(TESTLAB,outputs);
performance = perform(net,TESTLAB,outputs)



% close all; figure, plotperform(tr)
% close all; figure, plottrainstate(tr)
figure, plotconfusion(TESTLAB,outputs)
figure, ploterrhist(errors)

clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
AMX AMXCASE AMXCTRL AMXUSNP AMXTRCASE AMXTRCTRL AMXTECASE AMXTECTRL...
SRR TRAIN TEST TRAINMX TRAINLAB TESTMX TESTLAB NN net NNPerf 




%% GET AGES OF PEOPLE IN TRAINING & TESTING MATRICES


% GET PHENOTYPE TABLE OF ONLY TRAINING GROUP, 
% WITH A PARTICIPANT ORDER MATCHING THE TRAINING MATRIX

AMXPHENTRAIN = [AMXTRCASE; AMXTRCTRL];

[C,i,j] = intersect(SRR.train,AMXPHENTRAIN.SRR);
disp(numel(C)); disp(numel(i)); disp(numel(j))

[Ai,Bi] = ismember(SRR.train,AMXPHENTRAIN.SRR);
Bj = Bi(Bi>0);

PHETRAIN = AMXPHENTRAIN(Bj,:);




% GET PHENOTYPE TABLE OF ONLY TEST GROUP, 
% WITH A PARTICIPANT ORDER MATCHING THE TEST MATRIX
AMXPHENTEST  = [AMXTECASE; AMXTECTRL];
[C,i,j] = intersect(SRR.test,AMXPHENTEST.SRR);
disp(numel(C)); disp(numel(i)); disp(numel(j))
[Ai,Bi] = ismember(SRR.test,AMXPHENTEST.SRR);
Bj = Bi(Bi>0);
PHETEST = AMXPHENTEST(Bj,:);






ACTIVATION = net(TRAINMX);
GUESSES    = round(ACTIVATION);
CLASSES    = vec2ind(ACTIVATION);
CORRECT    = sum(TRAINLAB == GUESSES)>0;



clearvars -except ADSP LOCI CASE CTRL PHEN USNP TRCASE TRCTRL TECASE TECTRL...
AMX AMXCASE AMXCTRL AMXUSNP AMXTRCASE AMXTRCTRL AMXTECASE AMXTECTRL...
SRR TRAIN TEST TRAINMX TRAINLAB TESTMX TESTLAB NN net NNPerf...
PHETRAIN PHETEST ACTIVATION GUESSES CLASSES CORRECT







%% PLOT FIT OF AGE BY ACTIVATION

AGE = PHETRAIN.AGE;

AGEjitter = (rand(numel(PHETRAIN.AGE),1)-.5).*1.6;




X = AGE + AGEjitter;

Y = (ACTIVATION(1,:))';

C = CORRECT';


XCASE = X(PHETRAIN.AD==1);

YCASE = Y(PHETRAIN.AD==1);

CCASE = C(PHETRAIN.AD==1);



XCTRL = X(PHETRAIN.AD==0);

YCTRL = Y(PHETRAIN.AD==0);

CCTRL = C(PHETRAIN.AD==0);



XCACO = [XCASE ; XCTRL];

YCACO = [YCASE ; YCTRL];

CCACO = [ (CCASE+1)  ;  (CCTRL.*-1)]+2;









%% PLOT HISTOGRAM DISTRIBUTION OF AGES

X = AGE + AGEjitter;
Y = (ACTIVATION(1,:))';
C = CORRECT';


% AGES 60 AND UNDER
XCASE = X(PHETRAIN.AD==1 & PHETRAIN.AGE<=60);
YCASE = Y(PHETRAIN.AD==1 & PHETRAIN.AGE<=60);
CCASE = C(PHETRAIN.AD==1 & PHETRAIN.AGE<=60);
PCA60 = sum(CCASE) / numel(CCASE);

XCTRL = X(PHETRAIN.AD==0 & PHETRAIN.AGE<=70);
YCTRL = Y(PHETRAIN.AD==0 & PHETRAIN.AGE<=70);
CCTRL = C(PHETRAIN.AD==0 & PHETRAIN.AGE<=70);
PCO60 = sum(CCTRL) / numel(CCTRL);


% AGES 61 - 70
XCASE = X(PHETRAIN.AD==1 & PHETRAIN.AGE>60 & PHETRAIN.AGE<=70);
YCASE = Y(PHETRAIN.AD==1 & PHETRAIN.AGE>60 & PHETRAIN.AGE<=70);
CCASE = C(PHETRAIN.AD==1 & PHETRAIN.AGE>60 & PHETRAIN.AGE<=70);
PCA70 = sum(CCASE) / numel(CCASE);

XCTRL = X(PHETRAIN.AD==0 & PHETRAIN.AGE>60 & PHETRAIN.AGE<=70);
YCTRL = Y(PHETRAIN.AD==0 & PHETRAIN.AGE>60 & PHETRAIN.AGE<=70);
CCTRL = C(PHETRAIN.AD==0 & PHETRAIN.AGE>60 & PHETRAIN.AGE<=70);
PCO70 = sum(CCTRL) / numel(CCTRL);




% AGES 71 - 80
XCASE = X(PHETRAIN.AD==1 & PHETRAIN.AGE>70 & PHETRAIN.AGE<=80);
YCASE = Y(PHETRAIN.AD==1 & PHETRAIN.AGE>70 & PHETRAIN.AGE<=80);
CCASE = C(PHETRAIN.AD==1 & PHETRAIN.AGE>70 & PHETRAIN.AGE<=80);
PCA80 = sum(CCASE) / numel(CCASE);

XCTRL = X(PHETRAIN.AD==0 & PHETRAIN.AGE>68 & PHETRAIN.AGE<=78);
YCTRL = Y(PHETRAIN.AD==0 & PHETRAIN.AGE>68 & PHETRAIN.AGE<=78);
CCTRL = C(PHETRAIN.AD==0 & PHETRAIN.AGE>68 & PHETRAIN.AGE<=78);
PCO80 = sum(CCTRL) / numel(CCTRL);



% AGES 81 - 89
XCASE = X(PHETRAIN.AD==1 & PHETRAIN.AGE>75 & PHETRAIN.AGE<89);
YCASE = Y(PHETRAIN.AD==1 & PHETRAIN.AGE>75 & PHETRAIN.AGE<89);
CCASE = C(PHETRAIN.AD==1 & PHETRAIN.AGE>75 & PHETRAIN.AGE<89);
PCA89 = sum(CCASE) / numel(CCASE);

XCTRL = X(PHETRAIN.AD==0 & PHETRAIN.AGE>75 & PHETRAIN.AGE<80);
YCTRL = Y(PHETRAIN.AD==0 & PHETRAIN.AGE>75 & PHETRAIN.AGE<80);
CCTRL = C(PHETRAIN.AD==0 & PHETRAIN.AGE>75 & PHETRAIN.AGE<80);
PCO89 = sum(CCTRL) / numel(CCTRL);



% AGES 90 AND ABOVE
XCASE = X(PHETRAIN.AD==1 & PHETRAIN.AGE>=90);
YCASE = Y(PHETRAIN.AD==1 & PHETRAIN.AGE>=90);
CCASE = C(PHETRAIN.AD==1 & PHETRAIN.AGE>=90);
PCA90 = sum(CCASE) / numel(CCASE);

XCTRL = X(PHETRAIN.AD==0 & PHETRAIN.AGE>=80);
YCTRL = Y(PHETRAIN.AD==0 & PHETRAIN.AGE>=80);
CCTRL = C(PHETRAIN.AD==0 & PHETRAIN.AGE>=80);
PCO90 = sum(CCTRL) / numel(CCTRL);




CA = [PCA60 PCA70 PCA80 PCA89 PCA90];
CO = [PCO60 PCO70 PCO80 PCO89 PCO90];




Yi = [CA ; CO];

map =  [
        .0  .7  .1 ;
        .9  .5  .1 ;
        .9  .1  .6 ;
        .1  .5  .6 ;
        ];

Ci = {map(4,:), map(2,:)};




close all
fh1=figure('Units','normalized','OuterPosition',[.05 .05 .8 .7],'Color','w','MenuBar','none');
hax1 = axes('Position',[.05 .05 .9 .9],'Color','none');

h = superbar(Yi', 'BarFaceColor', Ci, 'BarEdgeColor', [.5  .5  .5]);
for i = 1:numel(hax1.Children)
% hax1.Children(i).EdgeColor = [.5  .5  .5];
hax1.Children(i).LineWidth = 3;
end

title('Percent correct for CASE and CTRL across AGE');


% XCACO = [XCASE ; XCTRL];
% YCACO = [YCASE ; YCTRL];
% CCACO = [ (CCASE+1)  ;  (CCTRL.*-1)]+2;
% 
% close all; clc;
% fh1 = figure('Units','normalized','OuterPosition',[.1 .04 .7 .8],...
% 'Color','w');
% ax1 = axes('Position',[.09 .1 .84 .86],'Color','none');
% 
% 
% h1 = histogram(round(PHETRAIN.AGE(PHETRAIN.AD==1)));
% 
% hold on
% h2 = histogram(round(PHETRAIN.AGE(PHETRAIN.AD==0)));
% h2.FaceAlpha = .4;
% 
% xlabel('\fontsize{22} Age')
% ylabel('\fontsize{22} People')
% ax1.FontSize = 16;
% 
% 
% legend('CASE', 'CTRL', 'Location', 'NorthWest' );
















%% PLOT HISTOGRAM DISTRIBUTION OF AGES

close all; clc;
fh1 = figure('Units','normalized','OuterPosition',[.1 .04 .7 .8],...
'Color','w');
ax1 = axes('Position',[.09 .1 .84 .86],'Color','none');


h1 = histogram(round(PHETRAIN.AGE(PHETRAIN.AD==1)));

hold on
h2 = histogram(round(PHETRAIN.AGE(PHETRAIN.AD==0)));
h2.FaceAlpha = .4;

xlabel('\fontsize{22} Age')
ylabel('\fontsize{22} People')
ax1.FontSize = 16;


legend('CASE', 'CTRL', 'Location', 'NorthWest' );



%% PLOT FIGURE WITHOUT LINEAR FIT LINE

close all; clc;
fh1 = figure('Units','normalized','OuterPosition',[.05 .04 .9 .92],...
'Color','w');
ax1 = axes('Position',[.09 .08 .84 .88],'Color','none');


ph1 = scatter(XCACO, YCACO,  220, CCACO, '.');
map =  [
        .0  .7  .1 ;
        .9  .5  .1 ;
        .9  .1  .6 ;
        .1  .5  .6 ;
        ];
colormap(map)


cb1 = colorbar;
cbt = linspace(1,4,9);
cb1.Ticks = cbt(2:2:end);
cb1.TickLabels = {'TRUE CTRL','FALSE CTRL','FALSE CASE','TRUE CASE'};
cb1.FontSize = 14;

xlabel('\fontsize{22} Participant Age');
ylabel('\fontsize{22} Neural Net Activation');


ax1.XLim = [58 92];
ax1.FontSize = 16;




%% PLOT DATA WITH LINEAR FIT LINE
clc;


x = XCACO;

y = YCACO;

c = CCACO;

close all
fh2 = figure('Units','normalized','OuterPosition',[.05 .04 .9 .92],...
'Color','w');
ax2 = axes('Position',[.09 .08 .84 .88],'Color','none','FontSize',16);


[xData, yData] = prepareCurveData( x, y );
ft = fittype( 'poly1' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Robust = 'Bisquare';
[fitresult, gof] = fit( xData, yData, ft, opts );


h = plot( fitresult, xData, yData, 'predobs', 0.9 );
% legend( h, 'Y vs. X', 'linear fit', 'CI95 upper',...
%  'CI95 lower', 'Location', 'NorthEast' );
legend off
grid off

ax2.XLim = [59 91];
ax2.YLim = [0 1];
ax2.FontSize = 16;

h(1).MarkerSize = 10;

h(2).LineWidth = 4;
h(3).LineWidth = 2;
h(4).LineWidth = 2;


hold on
ph1 = scatter(x, y,  220, c, '.');
map =  [
        .0  .7  .1 ;
        .9  .5  .1 ;
        .9  .1  .6 ;
        .1  .5  .6 ;
        ];
colormap(map)


cb1 = colorbar;
cbt = linspace(1,4,9);
cb1.Ticks = cbt(2:2:end);
cb1.TickLabels = {'TRUE CTRL','FALSE CTRL','FALSE CASE','TRUE CASE'};
cb1.FontSize = 14;

xlabel('\fontsize{22} Participant Age');
ylabel('\fontsize{22} Neural Net Activation');


ax2.XLim = [58 92];
ax2.FontSize = 16;


