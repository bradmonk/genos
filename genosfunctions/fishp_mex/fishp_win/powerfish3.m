function [FISHPM, FISHPL, FISHPR] = powerfish3(aa, bb, cc, dd)
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


    jj = ac;         jj(jj<0) = 1; jj=jj+1;
    kk = T-ac;       kk(kk<0) = 1; kk=kk+1;
    mm = ab;         mm(mm<0) = 1; mm=mm+1;
    pp = T-ab;       pp(pp<0) = 1; pp=pp+1;
    qq = T;          qq(qq<0) = 1; qq=qq+1;


    if a>=ab
        pL = 1;
    else

        xL = 0:a;
        ffL = xL;         ffL(ffL<0) = 1; ffL=ffL+1;
        ggL = ac-xL;      ggL(ggL<0) = 1; ggL=ggL+1;
        hhL = ab-xL;      hhL(hhL<0) = 1; hhL=hhL+1;
        iiL = T+xL-ac-ab; iiL(iiL<0) = 1; iiL=iiL+1;

        L = -logf(ffL) - logf(ggL) - logf(hhL) - logf(iiL);
        Lmax = max(L);
        pL = exp( logf(jj) + logf(kk) + logf(mm) + logf(pp) - logf(qq) +...
                  Lmax + log(sum(exp(L-Lmax))) );

    end


    if a*T*ab*ac == 0
        pR = 1;
    else

        xR = a:(min(ab,ac)-a+1)+a-1;
        ffR = xR;         ffR(ffR<0) = 1; ffR=ffR+1;
        ggR = ac-xR;      ggR(ggR<0) = 1; ggR=ggR+1;
        hhR = ab-xR;      hhR(hhR<0) = 1; hhR=hhR+1;
        iiR = T+xR-ac-ab; iiR(iiR<0) = 1; iiR=iiR+1;

        R = -logf(ffR) - logf(ggR) - logf(hhR) - logf(iiR);	
        Rmax = max(R);
        pR = exp( logf(jj) + logf(kk) + logf(mm) + logf(pp) - logf(qq) +...
                  Rmax + log(sum(exp(R-Rmax))) );
    end



    pMin = min([pL,pR]);


    if (numel(pL) ~=1) | (numel(pR) ~=1) | (numel(pMin) ~=1)
        keyboard
    end

    FISHPL(nn)  = pL;
    FISHPR(nn)  = pR;
    FISHPM(nn)  = pMin;
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
