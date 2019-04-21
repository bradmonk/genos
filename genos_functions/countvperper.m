function [CACOMX] = countvperper(SRR, AD, COHORTNUM, SNPCASE, SNPCTRL)
%% COUNT NUMBER OF VARIANTS PER PERSON


    si = size(SNPCASE,1);

    sz = size(SRR,1);

    CACOMX = zeros(sz,4);

    CACOMX(:,1) = SRR;

    CACOMX(:,2) = AD;

    CACOMX(:,3) = COHORTNUM;



    for nn = 1:si

        if ~isempty(SNPCASE{nn})
            ca = SNPCASE{nn}(:,1);
        else
            ca = 0;
        end

        if ~isempty(SNPCTRL{nn})
            co = SNPCTRL{nn}(:,1);
        else
            co = 0;
        end

        caco = [ca(:,1); co(:,1)];

        v = ismember(CACOMX(:,1),caco);

        CACOMX(:,4) = CACOMX(:,4) + v;

        if ~mod(nn,5000); disp(nn/si); end
    end


end