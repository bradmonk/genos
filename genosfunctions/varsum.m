function [G1num, G2num] = varsum(G1SNP, G1ID, G2SNP, G2ID)
%% varsum.m USAGE NOTES
%{
% 
% Syntax
% -----------------------------------------------------
%
%     [G1num, G2num] = varsum(G1SNP, G1ID, G2SNP, G2ID)
% 
% 
% Description
% -----------------------------------------------------
% 
%   varsum(G1SNP, G1ID, G2SNP, G2ID) 
%       takes G1SNP, G1ID and G2SNP, G2ID as arguments.
% 
%   SNP inputs are cells for each chr:pos when least 1 person has
%   a variant at that loci. Inside a cell is an Nx2 double array where 
%   column-1 are SRR IDs and column-2 are whether those people had alleles
%   that were HETEROZYGOUS ALTERNATE or HOMOZYGOUS ALTERNATE.
% 
%   ID arrays are Nx1 with SRR ID numbers
% 
%     
% 
% Example
% -----------------------------------------------------
% 
%     [G1num, G2num] = varsum(G1SNP, G1ID, G2SNP, G2ID)
% 
% 
% 
% See Also
% -----------------------------------------------------
%   http://bradleymonk.com/genos
%   http://bradleymonk.com/neuralnets
% 
% 
% Attribution
% -----------------------------------------------------
%   Created by: Bradley Monk
%   email: brad.monk@gmail.com
%   website: bradleymonk.com
%   2018.01.23
%
%}
%%


    %G1ID = uint64(G1ID);
    %G2ID = uint64(G2ID);

    sz = size(G1SNP,1);

    G1num = zeros(sz,1);

    G2num = zeros(sz,1);


    for nn = 1:sz

        if ~isempty(G1SNP{nn})

            v = G1SNP{nn}(:,1);
            x = G1SNP{nn}(:,2);

            [ai,~] = ismember(v, G1ID );

            G1num(nn) = sum( x(ai,:) );
        end



        if ~isempty(G2SNP{nn})

            v = G2SNP{nn}(:,1);
            x = G2SNP{nn}(:,2);

            [ai,~] = ismember(v, G2ID );

            G2num(nn) = sum( x(ai,:) );
        end


        if ~mod(nn,10000); disp(nn/sz); end
    end



end