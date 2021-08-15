%% Ensure neuron-type colors are identical
% (e.g., AWAL & AWAR should be the same color).

% Imports.
import Neurons.*;

% Initialize the models.
models = {
    'atlas_xx_rgb.mat'
    'atlas_xx_rgbw.mat'
    'atlas_xo_rgb.mat'
    'atlas_xo_rgbw.mat'
    };

% Initialize the body parts.
body = {
    'head'
    'tail'
    };

% Ensure neuron-type colors are identical.
for i = 1:length(models)
    
    % Load the atlas.
    m = models{i};
    atlas = [];
    load(m, 'atlas');
    
    % Fix each body part.
    for j = 1:length(body)
        
        % Sort the neurons.
        b = body{j};
        [atlas.(b).N, N_i] = sort(atlas.(b).N);
        atlas.(b).model.mu = atlas.(b).model.mu(N_i,:);
        atlas.(b).model.sigma = atlas.(b).model.sigma(:,:,N_i);
        
        % Replace identical neuron types with their mean color.
        neurons = atlas.(b).N;
        is_done = false(length(neurons), 1);
        for k = 1:length(neurons)
            
            % Did we already fix this neuron's color?
            if is_done(k)
                continue;
            end
            
            % Get the neuron name.
            neuron = neurons{k};
            [name, lr] = Hermaphrodite.stripLR(neuron);
            
            % Find identical neuron types.
            if isempty(lr) % there's only one neuron of this type
                continue;
            end
            type_i = find(startsWith(neurons, name));
            
            % Compute the mean neuron colors.
            mu = mean(atlas.(b).model.mu(type_i,4:end),1);
            sigma = mean(atlas.(b).model.sigma(4:end,4:end,type_i),3);
            
            % Fix the neuron colors.
            for l = 1:length(type_i)
                atlas.(b).model.mu(type_i(l),4:end) = mu;
                atlas.(b).model.sigma(4:end,4:end,type_i(l)) = sigma;
            end
            
            % We're done with this neuron type.
            is_done(type_i) = true;
        end
    end
    
    % Save the fixed model.
    save(m, 'atlas', '-append');
end
