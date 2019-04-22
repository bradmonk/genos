function [TRCaN,TRCoN,TECaN,TECoN,TRCaU,TRCoU,TECaU,TECoU] = ppsnpsum(...
          CASE,CTRL,USNP,TRCASESRR,TRCTRLSRR,TECASESRR,TECTRLSRR)
%%

    sz = size(CASE,1);

    TRCaN = zeros(sz,1);
    TRCoN = zeros(sz,1);
    TECaN = zeros(sz,1);
    TECoN = zeros(sz,1);
    TRCaU = zeros(sz,1);
    TRCoU = zeros(sz,1);
    TECaU = zeros(sz,1);
    TECoU = zeros(sz,1);

    TRCASESRR = sort(TRCASESRR);
    TRCTRLSRR = sort(TRCTRLSRR);
    TECASESRR = sort(TECASESRR);
    TECTRLSRR = sort(TECTRLSRR);



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
            ia = builtin('_ismemberhelper',v,TRCASESRR);
            TRCaN(nn) = sum( x(ia,:) >0 );


            % TESTING DATASET
            ia = builtin('_ismemberhelper',v,TECASESRR);
            TECaN(nn) = sum( x(ia,:) >0 );
        end



        % COUNT CTRL ALT VARIANTS
        if iNOTEMPTY_CTRL(nn)

            v = CTRL{nn}(:,1);
            x = CTRL{nn}(:,2);            

            % TRAINING DATASET
            ia = builtin('_ismemberhelper',v,TRCTRLSRR);
            TRCoN(nn) = sum( x(ia,:) >0 );

            % TESTING DATASET
            ia = builtin('_ismemberhelper',v,TECTRLSRR);
            TECoN(nn) = sum( x(ia,:) >0 );


        end



        % COUNT UNCALLED INSTANCES
        if iNOTEMPTY_USNP(nn)

            v = USNP{nn}(:,1);

            TRCaU(nn) = sum( builtin('_ismemberhelper',v,TRCASESRR) );
            TRCoU(nn) = sum( builtin('_ismemberhelper',v,TRCTRLSRR) );
            TECaU(nn) = sum( builtin('_ismemberhelper',v,TECASESRR) );
            TECoU(nn) = sum( builtin('_ismemberhelper',v,TECTRLSRR) );
        end


        if ~mod(nn,10000); disp(nn/sz); end
    end
    disp('done.'); disp(' ')





%{
    for nn = 1:sz


        % COUNT CASE ALT VARIANTS
        if iNOTEMPTY_CASE

            v = CASE{nn}(:,1);
            x = CASE{nn}(:,2);

            % TRAINING DATASET
            ai = ismember(v, TRCASESRR );
            TRCaN(nn) = sum( x(ai,:) );

            % TESTING DATASET
            ai = ismember(v, TECASESRR );
            TECaN(nn) = sum( x(ai,:) );
        end



        % COUNT CTRL ALT VARIANTS
        if iNOTEMPTY_CTRL

            v = CTRL{nn}(:,1);
            x = CTRL{nn}(:,2);            

            % TRAINING DATASET
            ai = ismember(v, TRCTRLSRR );
            TRCoN(nn) = sum( x(ai,:) );

            % TESTING DATASET
            ai = ismember(v, TECTRLSRR );
            TECoN(nn) = sum( x(ai,:) );


        end





        % COUNT UNCALLED INSTANCES
        if iNOTEMPTY_USNP(nn)

            v = USNP{nn}(:,1);

            TRCaU(nn) = numel(intersect(v, TRCASESRR ));
            TRCoU(nn) = numel(intersect(v, TRCTRLSRR ));
            TECaU(nn) = numel(intersect(v, TECASESRR ));
            TECoU(nn) = numel(intersect(v, TECTRLSRR ));
        end

        if ~mod(nn,10000); disp(nn/sz); end
    end
    %toc


%}

%----- THIS CODE WOULD SUFFICE TO COUNT UNCALLED INSTANCES -----
% ismem = @(v,r) sum(ismember(v, r ));
% TRCaU = cellfun(ismem,  USNP, repmat({TRCASESRR},size(USNP))  );
% TRCoU = cellfun(ismem,  USNP, repmat({TRCTRLSRR},size(USNP))  );
% TECaU = cellfun(ismem,  USNP, repmat({TECASESRR},size(USNP))  );
% TECoU = cellfun(ismem,  USNP, repmat({TECTRLSRR},size(USNP))  );
%---------------------------------------------------------------

%%
end