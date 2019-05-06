function find_ids(sp)


assignment_prob_ranks = sp.get_meta_data('assignment_prob_ranks');
assignments = sp.get_meta_data('assignments');
neurons = sp.get_meta_data('ids');

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

sp.add_deterministic_ids(ids);
sp.add_probabilistic_ids(ids_prob);

sp.add_meta_data('cols', cols);
sp.add_meta_data('rows', rows);

end