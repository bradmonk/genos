function [snpmatrix] = makesnp(snpcells)
%% makesnp.m USAGE NOTES
%{
% 
% Syntax
% -----------------------------------------------------
%
%     [snpmatrix] = makesnp(snpcells)
% 
% 
% Description
% -----------------------------------------------------
% 
%   makesnp(snpcells) takes snpcells as an argument
% 
%   snpcells has a cell for each chr:pos that has at least 1 person with
%   a variant at that loci. Inside a cell is an Nx2 double array where 
%   column-1 are SRR IDs and column-2 are whether those people had alleles
%   that were HETEROZYGOUS ALTERNATE or HOMOZYGOUS ALTERNATE.
% 
%   This function converts snpcells to a cell array to an Nx3 matrix.
%   Column-1 is the SRR ID. Column-2 is whether het/homozygous. And
%   Column-3 is a number representing which cell the data was from
%   in the cell array.
% 
%     
% 
% Example
% -----------------------------------------------------
% 
%     [snpmatrix] = makesnp(snpcells)
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


    % uint32 handles values below to 2^32 (4,294,967,295)
    % preallocate 200,000,000 (for a study with 10k people)

    snpmx = uint32(zeros(200000000 ,3)); 


    szc = size(snpcells,1);
    mm=1;
    for nn = 1:szc

        v = snpcells{nn};

        if ~isempty(v)
            c = numel(v(:,1));
            vn = [v repmat(nn,c,1)];
            snpmx( mm:mm+c-1, :) = uint32(vn);
            mm=mm+c;
        end

        if ~mod(nn,50000); disp(nn/szc); end
    end

    ok = snpmx(:,1)>0;

    snpmatrix = snpmx(ok,:);
    % snpmx will end up closer to 100,00,000

end