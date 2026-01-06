function Archive = UpdateArchive(Archive, NewPopulation, N)
    Union = [Archive, NewPopulation];
    
    % 宽松筛选：先保留前 1.5N 个最好的
    [FrontNo, ~] = NDSort(Union.objs, length(Union));
    [~, sort_idx] = sort(FrontNo);
    CandidateSize = min(length(Union), floor(1.5 * N));
    Candidate = Union(sort_idx(1:CandidateSize));
    
    % IDSS 精选：压回 N 个
    if length(Candidate) > N
        Archive = IDSS(Candidate, N);
    else
        Archive = Candidate;
    end
end