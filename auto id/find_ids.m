function sp_ids = find_ids(sp)


assignment_prob_ranks = sp.assignment_prob_ranks;
assignments = sp.assignments;

neurons = sp.id_neurons;

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

sp_ids = sp;
sp_ids.ids = ids;
sp_ids.ids_prob = ids_prob;
sp_ids.cols = cols;
sp_ids.rows = rows;
