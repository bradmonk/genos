function [pMin, pL, pR] = powerfish2(a, b, c, d)
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


sz = size(a,1);
tb = [a b c d];
FISHP  = zeros(sz,1);
FISHOR = zeros(sz,1);



logf    = log(0:100000);
logf(1) = 0;
logf    = cumsum(logf);

ab = a+b;
ac = a+c;
T  = a+b+c+d;

if a>=ab
	pL = 1;
else

	xL = 0:a;

    ff = xL;         ff(ff<0) = 1; ff=ff+1;
    gg = ac-xL;      gg(gg<0) = 1; gg=gg+1;
    hh = ab-xL;      hh(hh<0) = 1; hh=hh+1;
    ii = T+xL-ac-ab; ii(ii<0) = 1; ii=ii+1;
    jj = ac;         jj(jj<0) = 1; jj=jj+1;
    kk = T-ac;       kk(kk<0) = 1; kk=kk+1;
    mm = ab;         mm(mm<0) = 1; mm=mm+1;
    pp = T-ab;       pp(pp<0) = 1; pp=pp+1;
    qq = T;          qq(qq<0) = 1; qq=qq+1;
    %rr = xxx;        rr(rr<0) = 1; rr=rr+1;
    

	L = -logf(ff) - logf(gg) - logf(hh) - logf(ii);

	Lmax = max(L);

	pL = exp( logf(jj) + logf(kk) + logf(mm) + logf(pp) - logf(qq) +...
              Lmax + log(sum(exp(L-Lmax))) );
end





if a*T*ab*ac == 0
    pR = 1;
else
xR = a:(min(ab,ac)-a+1)+a-1;

ff = xR;         ff(ff<0) = 1; ff=ff+1;
gg = ac-xR;      gg(gg<0) = 1; gg=gg+1;
hh = ab-xR;      hh(hh<0) = 1; hh=hh+1;
ii = T+xR-ac-ab; ii(ii<0) = 1; ii=ii+1;
jj = ac;         jj(jj<0) = 1; jj=jj+1;
kk = T-ac;       kk(kk<0) = 1; kk=kk+1;
mm = ab;         mm(mm<0) = 1; mm=mm+1;
pp = T-ab;       pp(pp<0) = 1; pp=pp+1;
qq = T;          qq(qq<0) = 1; qq=qq+1;


    R = -logf(ff) - logf(gg) - logf(hh) - logf(ii);

    Rmax = max(R);

    pR = exp( logf(jj) + logf(kk) + logf(mm) + logf(pp) - logf(qq) +...
          Rmax + log(sum(exp(R-Rmax))) );
end



pMin = min([pL,pR]);

end

% function logfactoria = log_f(x)
% % compile a table
% persistent logftable
% 
% if isempty(logftable)
% 	logftable = log(0:100000);
% 	logftable(1) = 0;
% 	logftable = cumsum(logftable);
% end
% 
% % refer to table
% x(x<0) = 1;
% logfactoria = logftable(x+1);
% end
