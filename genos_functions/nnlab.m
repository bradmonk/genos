function [yx] = nnlab(y,varargin)


sz = size(y);

ylen = length(y);

ywid = min(sz);

ismx = ywid > 1;

nlabs = numel(unique(y));





if ~ismx

    yx = zeros(ylen,nlabs);

    G = findgroups(y);

    for i = 1:nlabs

        yx(G==i,i) = 1;

    end


    if sz(2)>sz(1)
        yx = yx';
    end

end



if ismx

    if sz(1)>sz(2)

        yx = zeros(ylen,1);

        for i = 1:ywid
            yx(y(:,i)==1) = i-1;
        end

    else
        yx = zeros(1,ylen);

        for i = 1:ywid
            yx(y(i,:)==1) = i-1;
        end

    end

end








end










