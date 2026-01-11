function Archive = UpdateArchive(Archive, NewPop, MaxSize)
% Update external archive for multimodal optimization
% Only keep non-dominated feasible solutions with good distribution

%% 1. Extract feasible solutions
if isempty(NewPop)
    return;
end
CV = sum(max(0, NewPop.cons), 2);
Feasible = NewPop(CV <= 1e-6);

if isempty(Feasible)
    return;
end

%% 2. Merge with existing archive
Combined = [Archive, Feasible];

%% 3. Duplicate elimination (based on decision space)
% Multi-modal optimization needs to keep solutions that are far in decision space
% Even if they have same objective values
PopDec = Combined.decs;
[~, uniIdx] = unique(round(PopDec, 4), 'rows');
Combined = Combined(uniIdx);

%% 4. Non-dominated sorting
% Only keep non-dominated solutions
FrontNo = NDSort(Combined.objs, 1);
Archive = Combined(FrontNo == 1);

%% 5. Capacity control using dual-space crowding distance
if length(Archive) > MaxSize
    Next = true(1, length(Archive));
    while sum(Next) > MaxSize
        RemainIdx = find(Next);
        TempPop = Archive(Next);
        
        % Normalize objectives and decisions
        Objs = TempPop.objs;
        Decs = TempPop.decs;
        f_min = min(Objs, [], 1);
        f_max = max(Objs, [], 1);
        objs_norm = (Objs - f_min) ./ (f_max - f_min + eps);
        
        d0_min = min(Decs, [], 1);
        d0_max = max(Decs, [], 1);
        decs_norm = (Decs - d0_min) ./ (d0_max - d0_min + eps);
        
        % Crowding distance in both spaces
        K = 3;
        dist_obj = pdist2(objs_norm, objs_norm);
        dist_obj = sort(dist_obj, 1);
        dobj = sum(dist_obj(1:min(K+1, size(dist_obj,1)), :), 1);
        
        dist_dec = pdist2(decs_norm, decs_norm);
        dist_dec = sort(dist_dec, 1);
        ddec = sum(dist_dec(1:min(K+1, size(dist_dec,1)), :), 1);
        
        % Combined score (smaller means more crowded)
        % We want to remove the most crowded one
        CrowdDis = dobj + ddec;
        [~, delIdx] = min(CrowdDis);
        
        Next(RemainIdx(delIdx)) = false;
    end
    Archive = Archive(Next);
end
end
