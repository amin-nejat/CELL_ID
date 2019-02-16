function ids = auto_id(sp, compartment)
% Automatic labeling of previously identified neurons based on mixture of
% Gaussian model.
%
% sp is supervoxel data structure that is output by MP algorithm
%
% Erdem


params = get_params();

% load('/home/mn2822/Desktop/WormAutoID/data/MoG/ruoxi_MP_11_YAa.mat');
load([params.mog_folder, '/aligned_worms_', compartment, '.mat']); %% LOAD ALIGNED TRAINING DATA HEAD
sigma=inv(nancov(reshape(permute(M,[1 3 2]),[size(M,1)*size(M,3) size(M,2)])));

X=[sp.mean.*[1 1 4] sp.color(:,1:3)]; % Make sure that z-pixels are scaled accordingly
Phat=zeros(size(M,1),size(X,1));
for j=1:size(M,3)
    Y=M(:,:,j);
    neurons=N;
    idx=find(~any(isnan(Y),2));
    neurons(any(isnan(Y),2))=[];
    Y(any(isnan(Y),2),:)=[];
    [beta,P,Xhat,cost]=regression_woc(X,Y,sigma,20);
    Phat(idx,:)=Phat(idx,:)+P;
    Xmatch=P*Xhat;
    
    %% COMMENT THESE OUT IF YOU DONT WANT TO VISUALIZE
%     figure(1)
%     subplot(3,1,1)
%     cla
%     hold on
%     for i=1:size(Y,1)
%         try
%         plot3(Xmatch(i,1),Xmatch(i,2),Xmatch(i,3),'o','MarkerSize',8,'MarkerFaceColor',Xmatch(i,[4 5 6])./max(Xmatch(:,[4 5 6]),[],1),'MarkerEdgeColor','k');
%         plot3(Y(i,1),Y(i,2),Y(i,3),'o','MarkerSize',8,'MarkerEdgeColor',sp.color(i,1:3)/max(sp.color(i,1:3)));
%         catch
%             continue;
%         end
%     end
%     set(gca,'Color','k');
%     title('Testing alignment to Training');
%     drawnow
%     
%     subplot(3,1,2)
%     cla
%     hold on
%     for i=1:size(Y,1)
%         try
%         plot3(Xmatch(i,1),Xmatch(i,2),Xmatch(i,3),'o','MarkerSize',8,'MarkerFaceColor',Xmatch(i,[4 5 6])./max(Xmatch(:,[4 5 6]),[],1),'MarkerEdgeColor','k');
%         text(Xmatch(i,1),Xmatch(i,2),Xmatch(i,3),neurons{i},'Color','yellow','FontSize',8);
%         catch
%             continue;
%         end
%     end
%     set(gca,'Color','k');
%     title('Testing (Assignment)')
%     drawnow
%     
%     subplot(3,1,3)
%     cla
%     hold on
%     for i=1:size(Y,1)
%         plot3(Y(i,1),Y(i,2),Y(i,3),'o','MarkerSize',8,'MarkerEdgeColor',sp.color(i,1:3)/max(sp.color(i,1:3)));
%         text(Y(i,1),Y(i,2),Y(i,3),neurons{i},'Color','white','FontSize',8);
%     end
%     set(gca,'Color','k');
%     title('Training (Ground truth)')
%     drawnow
    
    
    
end

% niter=2000;
% P = sinkhorn_logspace(LL, niter);
% 
% ncandidates = 7;
% for i=1:size(P,1);
%     prow = [P(i,:) 1-nansum(P(i,:))];
%     total_prob(i)=nansum(P(i,:));
%     [sort_prob,ind_prob] = sort(-prow);
%     ind_prob(ind_prob==length(prow))=NaN;
% 
%     assignment_prob(i)=find(P(i,:)==max(P(i,:)));
% 
%     assignment_prob_rank(i,:)=[ind_prob(1:ncandidates)];
%     assignment_prob_probs(i,:)=[-sort_prob(1:ncandidates)];
% 
%     det = find(assignments(i,:));
%     if(isempty(find(a(i,:))))
%         det2(i)=NaN;
%     else
%         det2(i)=find(a(i,:));;
%     end
%     if(isempty(det))
%         assignment_det(i)=NaN;
%     else
%         assignment_det(i)=det;
%     end
% end

[m,idx]=max(Phat./sum(Phat,1),[],1);%%Assignments
idx(isnan(m))=length(N)+1;

N{end+1} = '';
ids = N(idx);

%% FINAL VISUALIZATION
figure
hold on
for i=1:size(X,1)
    plot3(X(i,1),X(i,2),X(i,3),'o','MarkerSize',8,'MarkerFaceColor',X(i,[4 5 6])./max(X(:,[4 5 6]),[],1),'MarkerEdgeColor','k');
    text(X(i,1),X(i,2),X(i,3),N(idx(i)),'Color','white','FontSize',8);
end
set(gca,'Color','k')
title('Testing --- final assignment')
drawnow



end