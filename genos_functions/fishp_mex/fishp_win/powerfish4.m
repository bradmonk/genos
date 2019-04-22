function [FISHPM, FISHPL, FISHPR] = powerfish4(aa, bb, cc, dd)
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

sz = numel(aa);
% tb = [a b c d];
FISHPL  = zeros(sz,1);
FISHPR  = zeros(sz,1);
FISHPM  = zeros(sz,1);
% FISHOR = zeros(sz,1);



logf    = log(0:100000);
logf(1) = 0;
logf    = cumsum(logf);



for nn = 1:sz

    a = aa(nn);
    b = bb(nn);
    c = cc(nn);
    d = dd(nn);

    ab = a+b;
    ac = a+c;
    T  = a+b+c+d;


    jj = ac+1;         jj(jj<1) = 1;
    kk = T-ac+1;       kk(kk<1) = 1;
    mm = ab+1;         mm(mm<1) = 1;
    pp = T-ab+1;       pp(pp<1) = 1;
    qq = T+1;          qq(qq<1) = 1;


    if a>=ab
        pL = 1;
    else

        xL = 0:a;
        ffL = xL+1;         ffL(ffL<1) = 1;
        ggL = ac-xL+1;      ggL(ggL<1) = 1;
        hhL = ab-xL+1;      hhL(hhL<1) = 1;
        iiL = T+xL-ac-ab+1; iiL(iiL<1) = 1;

        L = -logf(ffL) - logf(ggL) - logf(hhL) - logf(iiL);
        Lmax = max(L);
        pL = exp( logf(jj) + logf(kk) + logf(mm) + logf(pp) - logf(qq) +...
                  Lmax + log(sum(exp(L-Lmax))) );

    end


    if a*T*ab*ac == 0
        pR = 1;
    else

        xR = a:(min(ab,ac)-a+1)+a-1;
        ffR = xR+1;         ffR(ffR<1) = 1;
        ggR = ac-xR+1;      ggR(ggR<1) = 1;
        hhR = ab-xR+1;      hhR(hhR<1) = 1;
        iiR = T+xR-ac-ab+1; iiR(iiR<1) = 1;

        R = -logf(ffR) - logf(ggR) - logf(hhR) - logf(iiR);	
        Rmax = max(R);
        pR = exp( logf(jj) + logf(kk) + logf(mm) + logf(pp) - logf(qq) +...
                  Rmax + log(sum(exp(R-Rmax))) );
    end


    FISHPL(nn)  = pL;
    FISHPR(nn)  = pR;
    FISHPM(nn)  = min([pL,pR]);

    if ~mod(nn,10000); disp(nn/sz); end
end


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
