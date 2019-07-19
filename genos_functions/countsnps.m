function [CaN,CoN,CaU,CoU] = countsnps(...
          CASE,CTRL,USNP,PHECASE,PHECTRL)
%%

    sz = size(CASE,1);

    CaN = zeros(sz,1);
    CoN = zeros(sz,1);
    CaU = zeros(sz,1);
    CoU = zeros(sz,1);

    PHECASE = sort(PHECASE);
    PHECTRL = sort(PHECTRL);

    iNOTEMPTY_CASE = ~cellfun('isempty',CASE);
    iNOTEMPTY_CTRL = ~cellfun('isempty',CTRL);
    iNOTEMPTY_USNP = ~cellfun('isempty',USNP);

    disp(' ')
    disp('Counting variants at each loci...')

    for nn = 1:sz


        % COUNT CASE ALT VARIANTS
        if iNOTEMPTY_CASE(nn)

            v = CASE{nn}(:,1);
            x = CASE{nn}(:,2);


            % TRAINING DATASET
            ia = builtin('_ismemberhelper',v,PHECASE);
            CaN(nn) = sum( x(ia,:) );


        end



        % COUNT CTRL ALT VARIANTS
        if iNOTEMPTY_CTRL(nn)

            v = CTRL{nn}(:,1);
            x = CTRL{nn}(:,2);            

            % TRAINING DATASET
            ia = builtin('_ismemberhelper',v,PHECTRL);
            CoN(nn) = sum( x(ia,:) );

        end



        % COUNT UNCALLED INSTANCES
        if iNOTEMPTY_USNP(nn)

            v = USNP{nn}(:,1);

            CaU(nn) = sum( builtin('_ismemberhelper',v,PHECASE) );
            CoU(nn) = sum( builtin('_ismemberhelper',v,PHECTRL) );
        end


        if ~mod(nn,10000); disp(nn/sz); end
    end
    disp('done.'); disp(' ')





%----- THIS CODE WOULD SUFFICE TO COUNT UNCALLED INSTANCES -----
% ismem = @(v,r) sum(ismember(v, r ));
% CaU = cellfun(ismem,  USNP, repmat({PHECASE},size(USNP))  );
% CoU = cellfun(ismem,  USNP, repmat({PHECTRL},size(USNP))  );
%---------------------------------------------------------------
%%
end