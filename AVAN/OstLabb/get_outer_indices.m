function outerIndices = get_outer_indices(NC,NR,NL,DIM)
%GET_OUTER_INDICES Summary of this function goes here
%   Detailed explanation goes here
    innerIndices = get_inner_indices(NC,NR,NL,DIM);
    
    if DIM == 2
        allIndices = [1:1:NR*NC];
        outerIndices = setdiff(allIndices,innerIndices);
    elseif DIM == 3
        allIndices = [1:1:NL*NR*NC];
        outerIndices = setdiff(allIndices,innerIndices);
    end
    
end

