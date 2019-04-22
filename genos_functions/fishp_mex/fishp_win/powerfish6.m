function [FISHPM, FISHPL, FISHPR] = powerfish6(aa, bb, cc, dd)
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

    a = aa;
    b = bb;
    c = cc;
    d = dd;

    b(a>=(a+b)) = b(a>=(a+b))+1;

    ab = a+b;
    ac = a+c;
    T  = a+b+c+d;


    %noR = a == 0;
    %noR = (a .* T .* ab .* ac) == 0;

    jj = ac+1;         jj(jj<1) = 1;
    kk = T-ac+1;       kk(kk<1) = 1;
    mm = ab+1;         mm(mm<1) = 1;
    pp = T-ab+1;       pp(pp<1) = 1;
    qq = T+1;          qq(qq<1) = 1;

    MINabac = min(ab,ac);


% keyboard

Lx = cell(sz,4);
Rx = cell(sz,4);

for nn = 1:sz

        Lx{nn,1} = 0:a(nn);
        Lx{nn,2} = ac(nn);
        Lx{nn,3} = ab(nn);
        Lx{nn,4} = T(nn);

        Rx{nn,1} = a(nn):(MINabac(nn)-a(nn)+1)+a(nn)-1;
        Rx{nn,2} = ac(nn);
        Rx{nn,3} = ab(nn);
        Rx{nn,4} = T(nn);

end






% eqOne = @(x) (x(x<1) = 1)

Fff = @(xx) xx+1;
Fgg = @(yy,xx) xx-yy+1;
Fhh = @(yy,xx) xx-yy+1;
Fii = @(xx,yy,zz,uu) xx+yy-zz-uu+1;

% FXlow = @(xx) xx<1;
% g = @(c) (integral(@(x) (x.^2 + c*x + 1),0,1));



A = cellfun(Fff, Lx(:,1),'UniformOutput',false);
B = cellfun(Fgg, Lx(:,1), Lx(:,2), 'UniformOutput',false);
C = cellfun(Fhh, Lx(:,1), Lx(:,2), 'UniformOutput',false);
D = cellfun(Fii, Lx(:,1), Lx(:,2), Lx(:,3), Lx(:,4),'UniformOutput',false);

E = cellfun(Fff, Rx(:,1),'UniformOutput',false);
F = cellfun(Fgg, Rx(:,1), Rx(:,2), 'UniformOutput',false);
G = cellfun(Fhh, Rx(:,1), Rx(:,2), 'UniformOutput',false);
H = cellfun(Fii, Rx(:,1), Rx(:,2), Rx(:,3), Rx(:,4),'UniformOutput',false);





for nn = 1:sz

    ffL = A{nn}; ffL(ffL<1) = 1; A{nn} = ffL;
    ggL = B{nn}; ggL(ggL<1) = 1; B{nn} = ggL;
    hhL = C{nn}; hhL(hhL<1) = 1; C{nn} = hhL;
    iiL = D{nn}; iiL(iiL<1) = 1; D{nn} = iiL;

    ffR = E{nn}; ffR(ffR<1) = 1; E{nn} = ffR;
    ggR = F{nn}; ggR(ggR<1) = 1; F{nn} = ggR;
    hhR = G{nn}; hhR(hhR<1) = 1; G{nn} = hhR;
    iiR = H{nn}; iiR(iiR<1) = 1; H{nn} = iiR;

end



FXLL = @(xx,yy,zz,uu) -logf(xx) - logf(yy) - logf(zz) - logf(uu);

FXpL = @(xx,yy,zz,uu) -logf(xx) - logf(yy) - logf(zz) - logf(uu);

LL = cellfun(FXLL, E, F, G, H,'UniformOutput',false);

FXLmax = cellfun(@max , LL,'UniformOutput',false);
Lmax = cellfun(@max , LL)';

FXLM = @(xx,yy) log(sum(exp(xx-yy)));

pLL = cellfun(FXLM , LL, FXLmax)';

exp( logf(jj) + logf(kk) + logf(mm) + logf(pp) - logf(qq) + Lmax + pLL)





    L = -logf(ffL) - logf(ggL) - logf(hhL) - logf(iiL);
    Lmax = max(L);
    pL = exp( logf(jj(nn)) + logf(kk(nn)) + logf(mm(nn)) + logf(pp(nn)) - logf(qq(nn)) +...
              Lmax + log(sum(exp(L-Lmax))) );


    R = -logf(ffR) - logf(ggR) - logf(hhR) - logf(iiR);	
    Rmax = max(R);
    pR = exp( logf(jj(nn)) + logf(kk(nn)) + logf(mm(nn)) + logf(pp(nn)) - logf(qq(nn)) +...
              Rmax + log(sum(exp(R-Rmax))) );


    FISHPL(nn)  = pL;
    FISHPR(nn)  = pR;
    FISHPM(nn)  = min([pL,pR]);

    if ~mod(nn,10000); disp(nn/sz); end

% end




for nn = 1:sz

        xL = Lx{nn};
        xR = Rx{nn};

        ffL = xL+1;                     ffL(ffL<1) = 1;
        ggL = ac(nn)-xL+1;              ggL(ggL<1) = 1;
        hhL = ab(nn)-xL+1;              hhL(hhL<1) = 1;
        iiL = T(nn)+xL-ac(nn)-ab(nn)+1; iiL(iiL<1) = 1;

        L = -logf(ffL) - logf(ggL) - logf(hhL) - logf(iiL);
        Lmax = max(L);
        pL = exp( logf(jj(nn)) + logf(kk(nn)) + logf(mm(nn)) + logf(pp(nn)) - logf(qq(nn)) +...
                  Lmax + log(sum(exp(L-Lmax))) );



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
