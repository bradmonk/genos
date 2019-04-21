function [RCASE, RCTRL, RPHEN] = makerandgroups(PHEN, CASE, CTRL)


R = randperm(size(PHEN,1));

Rca = R(1:size(PHEN(PHEN.AD==1,:),1));

Rco = R((size(PHEN(PHEN.AD==1,:),1)+1):size(R,2));

RCASEPHEN = PHEN(Rca,:);

RCTRLPHEN = PHEN(Rco,:);

RCASEPHEN.AD(:) = 1;

RCTRLPHEN.AD(:) = 0;

RPHEN = [RCASEPHEN; RCTRLPHEN];

RCASEID = RCASEPHEN.SRR;

RCTRLID = RCTRLPHEN.SRR;





sz = size(CASE,1);

RCASE = cell(size(CASE));
RCTRL = cell(size(CTRL));

% nn=19
for nn = 1:size(CASE,1)

    if ~isempty(CASE{nn})
        v = CASE{nn}(:,1);
        x = CASE{nn}(:,2);

        [ai,aj] = ismember(v, RCASEID );

        [bi,bj] = ismember(v, RCTRLID );


        if ~isempty(v(ai))
        RCASE{nn}(:,1:2) = [v(ai) x(ai)];
        end

        if ~isempty(v(bi))
        RCTRL{nn}(:,1:2) = [v(bi) x(bi)];
        end
        
    end

    if ~isempty(CTRL{nn})
        v = CTRL{nn}(:,1);
        x = CTRL{nn}(:,2);

        [ai,aj] = ismember(v, RCASEID );

        [bi,bj] = ismember(v, RCTRLID );


        if ~isempty(v(ai))
        RCASE{nn}(end+1:end+numel(v(ai)),1:2) = [v(ai) x(ai)];
        end

        if ~isempty(v(bi))
        RCTRL{nn}(end+1:end+numel(v(bi)),1:2) = [v(bi) x(bi)];
        end
        
    end


    if ~mod(nn,10000); disp(nn/sz); end
end
