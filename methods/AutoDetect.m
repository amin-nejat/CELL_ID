classdef AutoDetect < Singleton
   
   properties % Public Access
        options = []
        supervoxels = []
        
        fsize
        szext
        
        scale
   end
   
   methods(Access=private)
      function obj = AutoDetect()
        obj.options = optimoptions('fmincon', 'Display', 'off');
        obj.options.MaxFunEvals = 5000;
        obj.options.MaxIter = 5000;
        obj.options.StepTolerance = 1e-10;
      end
   end
   
   methods(Static)
      function obj = instance()
         persistent instance
         if isempty(instance)
            obj = AutoDetect();
            instance = obj;
         else
            obj = instance;
         end
      end
      
       function [c, ceq] = constrain_eigenvalues(x, lb, ub)
        % Nonlinear inequality constraints (eigenvalues) for fmincon.
        % Amin Nejat

            L = zeros(3, 3);
            L([1,2,3,5,6,9]) = x(4:9);
            variances = L*L';

            v = sort(eig(variances));
            c = -[v'-lb, ub-v'];
            
            ceq = [];
        end
        
        function resid = gaussian_cost(x, vol, norm_factor)
        % Calculating residual between a multi-color Gaussian function and the
        % original image.
        %
        % Amin Nejat

            x = x./norm_factor';
            
            mu = x(1:3);
            L = zeros(3, 3); L([1,2,3,5,6,9]) = x(4:9); cov = L*L';

            props = x(12:size(vol,4) + 11)*exp(x(10));
            baseline = x(size(vol,4) + 12:end-1)*x(11);
            sz = size(vol);
            recon = Utils.simulate_gaussian(sz(1:3), mu', cov, props', baseline', x(end));
            
            resid = sqrt(sum((vol(:) - recon(:)).^2));
        end
        
        function [shape, recon, mu, cov, colors] = get_gaussian(x, sz, norm)
        % Creating Gaussian function for a given set of parameters, used for
        % visualization.
        %
        % Amin Nejat

            x = x./norm;

            mu = x(1:3);
            L = zeros(3, 3); L([1,2,3,5,6,9]) = x(4:9); cov = L*L';

            colors = x(12:sz(4)+11)*exp(x(10));
            recon = Utils.simulate_gaussian(sz(1:3), mu, cov, colors, zeros(1,sz(4)), x(end));
            shape = recon(:,:,:,1)/colors(1);
        end
        
        
        function filter = get_filter(sz, sigma_factor, trunc)
        % Creating a Gaussian filter 2D, 3D, or 4D for smoothing an image.
        %
        % Amin Nejat

            fmu = (sz+1)/2;
            fsigma = diag(sz)./sigma_factor;
            pos = [];
            [pos(:, 1), pos(:, 2), pos(:, 3)] = ind2sub(sz, find(ones(sz(1:3))));
            filt = reshape(mvnpdf(pos, fmu, fsigma), sz);
            filt(filt < prctile(filt(:), trunc)) = 0;

            filter = filt/sqrt(sum(filt(:).*filt(:)));
        end
        
        function eval = evaluation(image, volume, mp_params, centers, scale, precision)
        % Evaluation of the segmentation method, the output contains shape
        % based evaluation (correlation and MSE) in the reconstruction
        % field and location based evaluation (TP, FP, TN, FN, Accuracy,
        % Precision, Recall, F1) in the centeres field.
        % Amin Nejat
            
            eval = [];
            eval.params.mp_params = mp_params;
            eval.params.precision = precision;
            
            sz = size(volume);
            reconstruction = image.get_3d_shape(sz,mp_params.hnsz);

            vec = @(x) x(:);


            % reconstruction based            
            eval.recon.mse = nanmean(vec((reconstruction - volume).^2));
            correlation = corrcoef(vec(reconstruction), vec(volume));
            eval.recon.cor = correlation(1,2);

            % shape based
%             segmentation = image.get_3d_segments(sz,mp_params.hnsz);

            % location based

            eval.centers.tp = [];
            eval.centers.fp = [];
            eval.centers.tn = [];
            eval.centers.fn = [];

            eval.centers.accuracy = [];
            eval.centers.f1 = [];
            eval.centers.precision = [];
            eval.centers.recall = [];
            eval.centers.mean_closest_distance = [];

            for i=1:size(image.get_positions(),1)
                positions = image.get_positions();
                positions = positions(1:i,:);

                distances = pdist2(positions.*scale', centers.*scale');
                distances(distances > precision) = Inf;

                if i == 1
                    tp = any(~isinf(distances));
                else
                    [~,~,assign] = munkres(distances);
                    [~,cols] = find(assign);
                    tp = length(cols);
                end

                fn = size(centers,1)-tp;
                fp = size(positions,1)-tp;
                tn = 0;

                eval.centers.tp(end+1) = tp;
                eval.centers.fp(end+1) = fp;
                eval.centers.tn(end+1) = tn;
                eval.centers.fn(end+1) = fn;

                eval.centers.accuracy(end+1) = (tp+tn)/(tp+fp+tn+fn);
                eval.centers.f1(end+1) = 2*tp/(2*tp+fp+fn);
                eval.centers.precision(end+1) = tp/(tp+fp);
                eval.centers.recall(end+1) = tp/(tp+fn);
                eval.centers.mean_closest_distance(end+1) = mean(min(pdist2(positions.*scale', centers.*scale')));
            end
        end
        
   end
   
   methods % Public Access
        function supervoxels = detect(obj, volume, filter, n_objects, cov_threshold, scale, exclusion)
        % Matching pursuit algorithm for finding the best mixture of Gaussians that
        % fits to input dataset. The algorithm runs greedy by subtracting off the
        % brightest Gaussian at each iteration.
        %
        % Amin Nejat
            obj.scale = scale;
            
            h = waitbar(0,'Initialize ...');
            
            obj.supervoxels = [];

            volume = double(volume);

            obj.szext = [size(volume), 1];
            obj.fsize = size(filter)-1;

            rho = Preprocess.filter_frame(volume, filter);

            N = 0;
            
            while N < n_objects && max(rho(:)) > 0.1
                try
                    waitbar((N+1)/n_objects,h,sprintf('%d%% completed ...',int16(100*(N+1)/n_objects)));
                catch
                    break;
                end
                
                rho_mag = max(rho, [], 4);
                [~, lmidx] = max(rho_mag(:));
                [x,y,z] = ind2sub(size(rho_mag), lmidx);
                
                bpatch = Utils.subcube(volume, [x,y,z], obj.fsize);

                [shape, sp] = obj.fit_gaussian(bpatch, squeeze(rho(x,y,z,:))', [x,y,z]);
                                
                fshape = imfilter(shape, filter, 'full');
                residual = Utils.placement(obj.szext(1:3), [x,y,z], fshape);

                for ch = 1: size(rho, 4)
                    rho(:,:,:,ch) = rho(:,:,:,ch)-residual*sp.color(ch);
                end
                
                
                exclusion_condition = isempty(obj.supervoxels) || min(sqrt(sum(((obj.supervoxels.mean-sp.mean).*obj.scale).^2, 2))) > exclusion;
                
                if min(eig(squeeze(sp.cov).*obj.scale)) > cov_threshold && exclusion_condition
                    obj.supervoxels = Utils.union_sp(obj.supervoxels, sp);
                    N = size(obj.supervoxels.mean, 1);
                end
            end
            try
                close(h);
            catch
                warning('The detection is canceled.');
            end
            supervoxels = obj.supervoxels;
        end


       
        function [shape, sp] = fit_gaussian(obj, bpatch, colorvec, absolute_position)

            norm = [obj.scale, ... % loc
                    sqrt([obj.scale(1),1,1,obj.scale(2),1,obj.scale(3)]), ... % covariance parameters
                    1, ... % color scale
                    1, ... % noise scale
                    ones(1, obj.szext(4)), ... % color
                    ones(1, obj.szext(4)), ... % noise
                    10000];
            
            eig_lb = [0.01,0.01,0.01];
            eig_ub = [5,5,5];
            

            nonlcon = @(x) AutoDetect.constrain_eigenvalues(x, eig_lb, eig_ub);

            colors = reshape(bpatch, [numel(bpatch)/size(bpatch, 4), size(bpatch, 4)]);
            colors(colors == 0) = eps;
            noisevec = prctile(colors, 10);

            colorvec(colorvec < 0) = eps;

            color_level = sum(colorvec(:));
            noise_level = abs(sum(noisevec(:)));

            colorvec = colorvec/color_level;
            noisevec = noisevec/noise_level;

            
            x0 =   double([(obj.fsize+1).*obj.scale, ... % loc
                sqrt([3,0,0,3,0,3]), ... % covariance parameters
                1.5*log(color_level), ... % color scale
                noise_level, ... % noise scale
                colorvec, ... % color
                noisevec, ... % noise
                0]);
            
            lb = [x0(1:3)-0.3*ones(1,3), zeros(1,6), 0, 0, zeros(1,obj.szext(4)), -ones(1,obj.szext(4)), 0];
            ub = [x0(1:3)+0.3*ones(1,3), 100*ones(1,6), 20, 0.1, ones(1,2*obj.szext(4)), 1000];
            
            
            A_eq = [zeros(1, 11), ones(1, obj.szext(4)), zeros(1, obj.szext(4)), 0; ...
                   [zeros(1, 11), zeros(1, obj.szext(4)), ones(1, obj.szext(4)), 1]];
            b_eq = [1; 1];

            
            f = @(x) AutoDetect.gaussian_cost(x', bpatch, norm);
            x_hat = fmincon(f, x0, [], [], A_eq, b_eq, lb, ub, nonlcon, obj.options);

            [shape, ~, ~, cov, col] = AutoDetect.get_gaussian(x_hat, [2*obj.fsize+1, obj.szext(4)], norm);
            bas = x_hat(12+obj.szext(4):11+2*obj.szext(4))*x_hat(11);

            relative_position = (x_hat(1:3)-x0(1:3))./norm(1:3);

            sp.mean = relative_position+absolute_position;
            sp.cov(1,:,:) = cov;
            sp.color = col;
            sp.baseline = bas;
            sp.truncation = x_hat(end)/norm(end);
        end
        
   end
   
end
