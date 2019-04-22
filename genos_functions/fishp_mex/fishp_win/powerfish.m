function [pval_left, pval_right] = powerfish(a, b, c, d)
% 
% a = round(abs(randn .* 20 + 300));
% b = round(abs(randn .* 20 + 250));
% c = round(abs(randn .* 20 + 30));
% d = round(abs(randn .* 20 + 40));
% 
% T = table([a;c],[b;d],'VariableNames',{'CASE','CTRL'},'RowNames',{'REFS','ALTS'});
% 
% 
% [~,fB,~] = fishertest(T);
% [~,fL,~] = fishertest(T,'Tail','left');
% [~,fR,~] = fishertest(T,'Tail','right');
% 
% 
% [pL,pR] = FastFisherExactTest(T.CASE(1), T.CTRL(1), T.CASE(2), T.CTRL(2));
% 
% disp([fB; fL; pL; fR; pR])
% 



PosC1 = a+b;
PosC2 = a+c;
Total = a+b+c+d;

if a>=PosC1
	pval_left = 1;
else
	x_left = 0:a;
	lst = -log_f(x_left)-log_f(PosC2-x_left)-log_f(PosC1-x_left)-log_f(Total+x_left-PosC2-PosC1);
	lst_max = max(lst);
	pval_left = exp(log_f(PosC2)+log_f(Total-PosC2)+log_f(PosC1)+log_f(Total-PosC1)-log_f(Total)+lst_max+log(sum(exp(lst-lst_max))));
end

if nargout>1
	if a*Total*PosC1*PosC2 == 0
		pval_right = 1;
	else
		minCT = min(PosC1,PosC2);
		n_lst = minCT-a+1;
		x_right = a:n_lst+a-1;
		lst = -log_f(x_right)-log_f(PosC2-x_right)-log_f(PosC1-x_right)-log_f(Total+x_right-PosC2-PosC1);
		lst_max = max(lst);
		pval_right = exp(log_f(PosC2)+log_f(Total-PosC2)+log_f(PosC1)+log_f(Total-PosC1)-log_f(Total)+lst_max+log(sum(exp(lst-lst_max))));
	end
end
end

function logfactoria = log_f(x)
% compile a table
persistent logftable

if isempty(logftable)
	logftable = log(0:100000);
	logftable(1) = 0;
	logftable = cumsum(logftable);
end

% refer to table
x(x<0) = 1;
logfactoria = logftable(x+1);
end
