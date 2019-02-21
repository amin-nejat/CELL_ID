function [assignments, assignment_det,  assignment_prob_ranks, assignment_prob_probs, rank_uncertain] = compute_assignments(LL);



%% Linear assignment

[assignments,~]=munkres(-LL);

ncandidates = 7;
niter=2000;
P = sinkhorn_logspace(LL, niter);


for i=1:size(P,1);
    prow = [P(i,:) 1-nansum(P(i,:))];
    total_prob(i)=nansum(P(i,:));
    [sort_prob,ind_prob] = sort(-prow);
    ind_prob(ind_prob==length(prow))=-1;

    assignment_prob_ranks(i,:)=[ind_prob(1:ncandidates)];
    assignment_prob_probs(i,:)=[-sort_prob(1:ncandidates)];

    det = find(assignments(i,:));

    if(isempty(det))
        assignment_det(i)=-1;
    else
        assignment_det(i)=det;
        
    end
end
for i=1:size(P,1)
  ind_det_in_probs = find(assignment_det(i)==assignment_prob_ranks(i,:));
  rest_indices = setdiff([1:ncandidates], ind_det_in_probs);
  ranks = assignment_prob_ranks(i,:);
  probs = assignment_prob_probs(i,:);
     
  assignment_prob_ranks(i, 1) = ranks(ind_det_in_probs);
  assignment_prob_ranks(i, 2:ncandidates)= ranks(rest_indices);
     
end


[~,rank_uncertain]=sort(assignment_prob_probs(:,1));

