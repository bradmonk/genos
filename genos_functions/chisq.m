function [CHIP, CHIOR] = chisq(CASEREFS, CASEALTS, CTRLREFS, CTRLALTS)



ta = [CASEREFS CASEALTS CTRLREFS CTRLALTS];
tb = reshape(ta',2,2,[]);

sz = size(tb,3);

CHIP   = zeros(sz,1);
CHIOR  = zeros(sz,1);


for nn = 1:sz

    x = tb(:,:,nn);

    n = sum(x(:));

    P_ref = sum(x(1,:)) ./ n;
    P_alt = sum(x(2,:)) ./ n;
    P_ca = sum(x(:,1)) ./ n;
    P_co = sum(x(:,2)) ./ n;

    E_refca = P_ref * P_ca * n;
    E_refco = P_ref * P_co * n;
    E_altca = P_alt * P_ca * n;
    E_altco = P_alt * P_co * n;

    d1 = (x(1,1) - E_refca)^2 / E_refca;
    d2 = (x(1,2) - E_refco)^2 / E_refco;
    d3 = (x(2,1) - E_altca)^2 / E_altca;
    d4 = (x(2,2) - E_altco)^2 / E_altco;

    X = d1+d2+d3+d4;

    or = x(1,1)*x(2,2)/x(1,2)/x(2,1);

    CHIP(nn) = 1-chi2cdf(X,1);

    CHIOR(nn) = or;

    if ~mod(nn,10000); disp(nn/sz); end
end






% FISHP  = zeros(sz,1);
% FISHOR = zeros(sz,1);
% tic
% for nn = 1:sz
% 
%     x = tb(:,:,nn);
% 
%     r1 = sum(x(1,:));
%     r2 = sum(x(2,:));
%     c1 = sum(x(:,1));
%     c2 = sum(x(:,2));
%     n = sum(x(:));
% 
% 
%     if min(r1,c1)<= min(r2,c2)
%         x11 = (0 : min(r1,c1))';
%     else
%         x22 = (0 : min(r2,c2))';
%         x12 = c2 - x22;
%         x11 = r1 - x12;
%     end       
% 
%     p1 = hygepdf(x(1,1),n,r1,c1);
%     p2 = hygepdf(x11,n,r1,c1);
%     p = sum(p2(p2 < p1+10*eps(p1)));
% 
%     or = x(1,1)*x(2,2)/x(1,2)/x(2,1);
% 
%     FISHP(nn)  = p;
%     FISHOR(nn) = or;
% 
%     if ~mod(nn,10000); disp(nn/sz); end
% end
% toc
% [FISHP CHIP]



end