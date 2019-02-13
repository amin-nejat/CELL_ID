function ids = auto_id(sp, compartment)
% Automatic labeling of previously identified neurons based on mixture of
% Gaussian model.
%
% sp is supervoxel data structure that is output by MP algorithm
%
% Erdem

positions=sp.mean;
D=squareform(pdist(bsxfun(@times,positions,[1 1 4])));
for j=1:size(sp.color,1)
    worm_colors{j}=sp.color(:,[1 2 3]);
    worm_distances{j}=D(j,:);
    worm_positions{j}=positions(j,:);
end


%% Training ID likelihoods
params = get_params();
load([params.mog_folder, '/mog_model_', compartment, '.mat']);
neurons=N;
N=length(worm_distances);
for k=1:N
    k
    for j=1:length(mog)
        try
            LL(k,j)=nansum(log(pdf(mog{j},[worm_distances{k}((1:N)~=k)' worm_colors{k}((1:N)~=k,:)]))) + log(mvnpdf(worm_colors{k}(k,:),cd{j}.mean,cd{j}.cov));
        catch
        end
    end
end


%% Linear assignment

[assignments,cost]=munkres(bsxfun(@times,LL,1./nansum(LL,2)));

ids = repmat({''}, size(sp.mean,1),1);
[rows,cols] = find(assignments);
ids(rows) = neurons(cols);



%% Visualize ID's

% figure('units','normalized','outerposition',[0 0 1 1])
% imagesc(squeeze(max(v(:,:,:,[4 3 1])./max(vec(v(:,:,:,[4 3 1]))),[],3)))
% hold on
% for j=1:size(positions,1)
%     try
%         k=find(assignments(j,:));
%         text(positions(j,2),positions(j,1),neurons{k},'Color','white','FontSize',8)
%     catch
%         text(positions(j,2),positions(j,1),'X','Color','red','FontSize',8)
%     end
% end
% [a,b,c]=fileparts(filenames{i});
% [a,b,c]=fileparts(a);
% title(b)
% export_fig([b '.png']);

%% Color comparison QC
for k=1:length(rows)
    
    match(k,1)=cols(k);
    match(k,2)=rows(k);
    Icolor(k,1,[1 2 3])=cd{cols(k)}.mean;
    Icolor(k,2,[1 2 3])=sp.color(rows(k),[1,2,3]);
%     Icolor(k,3,[1 2 3])=colors(j,:);
end

figure
imagesc(Icolor./max(Icolor,[],1))
title('Column1: Training color mean, Column 2: matched neuron color readout, Column 3: matched neuron MP color')
figure
scatter(Icolor(:,1),Icolor(:,2),'.')
title('X: Training color mean, Y: matched neuron color')

end