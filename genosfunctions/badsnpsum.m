function [TRCaN,TRCoN,TECaN,TECoN,TRCaU,TRCoU,TECaU,TECoU] = badsnpsum(...
          CASE,CTRL,USNP,TRCASESRR,TRCTRLSRR,TECASESRR,TECTRLSRR)
%%

% keyboard

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



    CACO = cellfun(@vertcat,CASE,CTRL,'UniformOutput',false);


    iNOTEMPTY_CACO = ~cellfun('isempty',CACO);
    iNOTEMPTY_USNP = ~cellfun('isempty',USNP);




    disp(' ')
    disp('Counting variants at each loci...')

    for nn = 1:sz


        ids = CACO{nn}(:,1);
        het = CACO{nn}(:,2);
        unk = USNP{nn};

        % COUNT UNCALLED INSTANCES
        if iNOTEMPTY_CACO(nn)

        trca = builtin('_ismemberhelper',ids,TRCASESRR);
        trco = builtin('_ismemberhelper',ids,TRCTRLSRR);
        teca = builtin('_ismemberhelper',ids,TECASESRR);
        teco = builtin('_ismemberhelper',ids,TECTRLSRR);

        TRCaN(nn) = sum( het(trca,:) );
        TRCoN(nn) = sum( het(trco,:) );
        TECaN(nn) = sum( het(teca,:) );
        TECoN(nn) = sum( het(teco,:) );

        end
        
        % COUNT UNCALLED INSTANCES
        if iNOTEMPTY_USNP(nn)

            TRCaU(nn) = sum( builtin('_ismemberhelper',unk,TRCASESRR) );
            TRCoU(nn) = sum( builtin('_ismemberhelper',unk,TRCTRLSRR) );
            TECaU(nn) = sum( builtin('_ismemberhelper',unk,TECASESRR) );
            TECoU(nn) = sum( builtin('_ismemberhelper',unk,TECTRLSRR) );

        end


        if ~mod(nn,10000); disp(nn/sz); end
    end
    disp('done.'); disp(' ')


%----- THIS CODE WOULD SUFFICE TO COUNT UNCALLED INSTANCES -----
% ismem = @(v,r) sum(ismember(v, r ));
% TRCaU = cellfun(ismem,  USNP, repmat({TRCASESRR},size(USNP))  );
% TRCoU = cellfun(ismem,  USNP, repmat({TRCTRLSRR},size(USNP))  );
% TECaU = cellfun(ismem,  USNP, repmat({TECASESRR},size(USNP))  );
% TECoU = cellfun(ismem,  USNP, repmat({TECTRLSRR},size(USNP))  );
%---------------------------------------------------------------

%%
end