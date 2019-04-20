function [] = gcounts(MX)

    format shortG

    [Gi,GNAMES] = findgroups( MX );

    GCOUNTS = splitapply(@numel, MX , Gi );

    if iscell(GNAMES)

        fprintf('\t Group\t\t  Count \n\t----------\t ---------\n')
        disp([repmat('    ',numel(GNAMES),1), char(GNAMES),...
              repmat('    ',numel(GNAMES),1), num2str(GCOUNTS)])

    else

        fprintf('\t Group\t\t  Count \n\t----------\t ---------\n')
        fprintf('%12.10g   |   %4.12g \n',[GNAMES GCOUNTS]')

    end

end