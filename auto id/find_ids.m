function [ids, ids_prob, cols, rows] = find_ids(compartment, assignment_prob_ranks, assignments, method)
params = get_params();

if strcmp(method, 'mog')
    load([params.mog_folder, '/mog_model_', compartment, '.mat']);
else
    load([params.mog_folder, '/aligned_worms_', compartment, '.mat']);
end
neurons = N;

ncandidates = size(assignment_prob_ranks, 2);
ids = repmat({'NaN'}, size(assignments,1),1);
ids_prob = repmat({'NaN'}, size(assignments,1),ncandidates);

[rows,cols] = find(assignments);
ids(rows) = neurons(cols);

for i=1:size(assignment_prob_ranks,1)
    for j=1:ncandidates
        if(assignment_prob_ranks(i,j)>0)
            
            ids_prob{i,j}=neurons{assignment_prob_ranks(i,j)};
        end
    end
end

ids(cellfun(@isempty, ids)) = {'NaN'};
ids_prob(cellfun(@isempty, ids_prob)) = {'NaN'};