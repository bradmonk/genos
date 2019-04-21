function [varargout] = NNstudystats(net,TRAINLAB,TRAINMX,TESTMX,TESTLAB,Pfilter,HICI,LOCI)

yguess = net(TRAINMX);

[CONFIDENCE,GUESS] = max(yguess);
MEANCONF = mean(yguess);


[ACTUAL,~,~] = find(TRAINLAB);
AGC = [ACTUAL GUESS' CONFIDENCE'];
HITS = AGC(:,1)==AGC(:,2);
mu = mean(HITS);
fprintf('%.2f  Percent correct on all training data \n',mu*100)


hi = AGC(:,3)>.85 | AGC(:,3)<.15;
CONFHITS = HITS(hi);
mu = mean(CONFHITS);
pn = sum(hi) / numel(hi);
disp(' ')
fprintf('%.2f  Percent correct on high-confidence training data \n',mu*100)
fprintf('%.2f  Percent of training data registered high-confidence \n',pn*100)



yguess = net(TESTMX);
[CONFIDENCE,GUESS] = max(yguess);
[ACTUAL,~,~] = find(TESTLAB);
AGC = [ACTUAL GUESS' CONFIDENCE'];
HITS = AGC(:,1)==AGC(:,2);
mu = mean(HITS);
PERF.all = sprintf('\n\n%.2f  Percent correct on all hold-out test data',mu*100);
disp(PERF.all)

hi = AGC(:,3)>.85 | AGC(:,3)<.15;
CONFHITS = HITS(hi);
mu = mean(CONFHITS);
pn = sum(hi) / numel(hi);
disp(' ')
PERF.hic = sprintf('%.2f  Percent correct on high-confidence hold-out test data',mu*100);
PERF.pct = sprintf('%.2f  Percent of hold-out test data registered high-confidence \n',pn*100);
disp(PERF.hic)
disp(PERF.pct)

disp(' ')
disp('**high-confidence definition: .15 > nnet_guess >.85')


end