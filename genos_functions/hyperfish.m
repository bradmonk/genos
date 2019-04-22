function [FISHPM, FISHPL, FISHPR, FISHOR] = hyperfish(aa, bb, cc, dd)
%% hyperfish(CASEREF,CTRLREF,CASEALT,CTRLALT)
% 
% a = LOCI.TRCASEREF(1:10000);
% b = LOCI.TRCTRLREF(1:10000);
% c = LOCI.TRCASEALT(1:10000);
% d = LOCI.TRCTRLALT(1:10000);
% 
% 
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

%--------------------------------------------------------------------------

sz = numel(aa);

% tb = [CASEREFS CASEALTS CTRLREFS CTRLALTS];
% tb = [aa cc bb dd];


FISHPL  = zeros(sz,1);
FISHPR  = zeros(sz,1);
FISHPM  = zeros(sz,1);

FISHOR = aa.*dd./cc./bb;


logf    = log(0:100000);
logf(1) = 0;
logf    = cumsum(logf);

a = aa;
b = bb;
c = cc;
d = dd;

b(a>=(a+b)) = b(a>=(a+b))+1;

ab = a+b;
ac = a+c;
T  = a+b+c+d;


jj = ac+1;         jj(jj<1) = 1;
kk = T-ac+1;       kk(kk<1) = 1;
mm = ab+1;         mm(mm<1) = 1;
pp = T-ab+1;       pp(pp<1) = 1;
qq = T+1;          qq(qq<1) = 1;

MINabac = min(ab,ac);



for nn = 1:sz

        xL = 0:a(nn);

        ffL = xL+1;                     ffL(ffL<1) = 1;
        ggL = ac(nn)-xL+1;              ggL(ggL<1) = 1;
        hhL = ab(nn)-xL+1;              hhL(hhL<1) = 1;
        iiL = T(nn)+xL-ac(nn)-ab(nn)+1; iiL(iiL<1) = 1;

        L = -logf(ffL) - logf(ggL) - logf(hhL) - logf(iiL);
        Lmax = max(L);
        pL = exp( logf(jj(nn)) + logf(kk(nn)) + logf(mm(nn)) + logf(pp(nn)) - logf(qq(nn)) +...
                  Lmax + log(sum(exp(L-Lmax))) );

        
        xR = a(nn):(MINabac(nn)-a(nn)+1)+a(nn)-1;

        ffR = xR+1;                         ffR(ffR<1) = 1;
        ggR = ac(nn)-xR+1;                  ggR(ggR<1) = 1;
        hhR = ab(nn)-xR+1;                  hhR(hhR<1) = 1;
        iiR = T(nn)+xR-ac(nn)-ab(nn)+1;     iiR(iiR<1) = 1;

        R = -logf(ffR) - logf(ggR) - logf(hhR) - logf(iiR);	
        Rmax = max(R);
        pR = exp( logf(jj(nn)) + logf(kk(nn)) + logf(mm(nn)) + logf(pp(nn)) - logf(qq(nn)) +...
                  Rmax + log(sum(exp(R-Rmax))) );



        

    FISHPL(nn)  = pL;
    FISHPR(nn)  = pR;
    FISHPM(nn)  = min([pL,pR]);

    if ~mod(nn,10000); disp(nn/sz); end
end


end

