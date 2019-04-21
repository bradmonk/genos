function [FISHP, FISHOR] = ffish(aa, bb, cc, dd)
%% fastfish(CASEREF,CTRLREF,CASEALT,CTRLALT)
% 
% a = LOCI.TRCASEREF(1:10000);
% b = LOCI.TRCTRLREF(1:10000);
% c = LOCI.TRCASEALT(1:10000);
% d = LOCI.TRCTRLALT(1:10000);
% 
% [FISHPM, FISHPL, FISHPR, FISHOR] = fastfish(a, b, c, d)
% 
% 
%--------------------------------------------------------------------------
%{
keyboard

sz = numel(aa);

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
MINI = (MINabac - a + 1) + a - 1;

fLR = logf(jj) + logf(kk) + logf(mm) + logf(pp) - logf(qq);


keyboard



for nn = 1:sz

    xL = 0:a(nn);
    ffL = xL+1;                     ffL(ffL<1) = 1;
    ggL = ac(nn)-xL+1;              ggL(ggL<1) = 1;
    hhL = ab(nn)-xL+1;              hhL(hhL<1) = 1;
    iiL = T(nn)+xL-ac(nn)-ab(nn)+1; iiL(iiL<1) = 1;

    L = -logf(ffL) - logf(ggL) - logf(hhL) - logf(iiL);
    Lmax = max(L);


    pL = exp( fLR(nn) + Lmax + log(sum(exp(L-Lmax))) );


    xR = a(nn):MINI(nn);
    ffR = xR+1;                         ffR(ffR<1) = 1;
    ggR = ac(nn)-xR+1;                  ggR(ggR<1) = 1;
    hhR = ab(nn)-xR+1;                  hhR(hhR<1) = 1;
    iiR = T(nn)+xR-ac(nn)-ab(nn)+1;     iiR(iiR<1) = 1;
    R = -logf(ffR) - logf(ggR) - logf(hhR) - logf(iiR);	
    Rmax = max(R);
    pR = exp( fLR(nn) + Rmax + log(sum(exp(R-Rmax))) );


    FISHPL(nn)  = pL;
    FISHPR(nn)  = pR;
    FISHPM(nn)  = min([pL,pR]);

    if ~mod(nn,10000); disp(nn/sz); end
end
%}

%%

    % fastfish(CASEREF,CTRLREF,CASEALT,CTRLALT)

    sz = numel(aa);

    %sz = size(CASEREFS,1);
    %tb = [aa cc bb dd];


    rr1 = sum([aa,cc],2);
    rr2 = sum([bb,dd],2);
    cc1 = sum([aa,bb],2);
    cc2 = sum([cc,dd],2);
    tot = sum([aa,bb,cc,dd],2);

    MINr1c1 = min(rr1,cc1);
    MINr2c2 = min(rr2,cc2);

    ifA = MINr1c1 <= MINr2c2;


    p1 = hygepdf(aa,tot,rr1,cc1);
    ep1 = 10.*eps(p1);
    e1 = p1+ep1;
        
    FISHP  = zeros(sz,1);
    FISHOR = zeros(sz,1);


    for nn = 1:sz

        if ifA(nn)
            V = (0 : MINr1c1(nn));
        else
            V = rr1(nn) - (cc2(nn) - (0 : MINr2c2(nn)));
        end       

        p2 = hygepdf(V,tot(nn),rr1(nn),cc1(nn));

        p = sum( p2(p2 < e1(nn)) );

        %or = x(1,1)*x(2,2)/x(1,2)/x(2,1);
        or = 1;

        FISHP(nn)  = p;
        FISHOR(nn) = or;

        if ~mod(nn,10000); disp(nn/sz); end
    end







%%
end

