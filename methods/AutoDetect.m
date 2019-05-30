classdef AutoDetect < Singleton
   
   properties % Public Access
        options = []
        supervoxels = []
        
        fsize
        trunc
        szext
   end
   
   methods(Access=private)
      function obj = AutoDetect()
        obj.options = optimoptions('fmincon', 'Display', 'off');
        obj.options.MaxFunEvals = 5000;
        obj.options.MaxIter = 5000;
        obj.options.StepTolerance = 1e-8;
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
      
      
      function [c, ceq] = constrain_eigenvalues(x, norm, init, bound)
        % Nonlinear inequality constraints (eigenvalues) for fmincon.
        % Amin Nejat
            x = x./norm;

            L = zeros(3, 3);
            L([1,2,3,5,6,9]) = x(1: 6);
            variances = L'*L;

            v = sort(eig(variances));
            c = -[v'-init+bound, bound-v'+init];
            
            ceq = [];
        end
        
        function resid = fit_gaussian_fixed_cov(x, vol, norm, trunc, cov)
        % Calculating residual between a multi-color Gaussian function and the
        % original image where the covariance of the Gaussian is fixed.
        %
        % Amin Nejat

            x = x./norm';

            mu = x(1:3);

            props = x(6:size(vol,4) + 5)*exp(x(4));
            baseline = x(size(vol,4) + 6:end)*x(5);
            sz = size(vol);

            recon = Utils.simulate_gaussian(sz(1:3), ...
                mu', cov, props', baseline', zeros(1, size(vol, 4)), trunc);

            resid = sqrt(sum((vol(:) - recon(:)).^2));
        end
        
        function resid = fit_gaussian_cov(x, vol, norm, trunc, mu, color_props, baseline_props)
        % Calculating residual between a multi-color Gaussian for a given set of
        % parameters and the original image. In this function covariance of the
        % Gaussian is the only variable, rest of the parameters are fixed. Used by 
        % fmincon to find the best fit Gaussian function.
        %
        % Amin Nejat

            x = x./norm';

            L = zeros(3, 3);
            L([1,2,3,5,6,9]) = x(1:6);
            variances = L*L';

            props = color_props*exp(x(7));
            baseline = baseline_props*x(8);
            sz = [size(vol, 1), size(vol, 2), size(vol, 3), size(vol, 4)];

            recon = Utils.simulate_gaussian(sz(1: 3), ...
                mu', variances, props', baseline', zeros(1, size(vol, 4)), trunc);

            resid = sqrt(sum((vol(:) - recon(:)).^2));
        end
        
        function [shape, recon, mu, cov, colors] = get_gaussian_cov(x, sz, norm, trunc, mu, color_props)
        % Creating Gaussian function for a given set of parameters, used for
        % visualization.
        %
        % Amin Nejat

            x = x./norm;

            L = zeros(3, 3);
            L([1,2,3,5,6,9]) = x(1:6);
            cov = L*L';

            colors = color_props*exp(x(7));

            recon = Utils.simulate_gaussian(sz, ...
                mu, cov, colors, zeros(1, sz(4)), zeros(1, sz(4)), trunc);

            shape = recon(:, :, :, 1)/colors(1);
        end
        
        function filter = get_filter(sz, sigma_factor, truncate_percentage)
        % Creating a Gaussian filter 2D, 3D, or 4D for smoothing an image.
        %
        % Amin Nejat

            fmu = (sz+1)/2;
            fsigma = diag(sz)./sigma_factor;
            pos = [];
            [pos(:, 1), pos(:, 2), pos(:, 3)] = ind2sub(sz, find(ones(sz(1:3))));
            filt = reshape(mvnpdf(pos, fmu, fsigma), sz);
            filt(filt < prctile(filt(:), truncate_percentage)) = 0;

            filter = filt/sqrt(sum(filt(:).*filt(:)));
        end
        
   end
   
   methods % Public Access
        function supervoxels = detect(obj, volume, filter, n_objects, trunc, cov_threshold)
        % Matching pursuit algorithm for finding the best mixture of Gaussians that
        % fits to input dataset. The algorithm runs greedy by subtracting off the
        % brightest Gaussian at each iteration.
        %
        % Amin Nejat
            
            h = waitbar(0,'Initialize ...');
            
            obj.supervoxels = [];

            volume = double(volume);

            obj.szext = [size(volume), 1];
            obj.fsize = size(filter)-1;
            obj.trunc = trunc;

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

                [shape, sp, ~] = obj.fit_gaussian(bpatch, squeeze(rho(x,y,z,:))', [x,y,z]);

                fshape = imfilter(shape, filter, 'full');
                residual = Utils.placement(obj.szext(1:3), [x,y,z], fshape);

                for ch = 1: size(rho, 4)
                    rho(:,:,:,ch) = rho(:,:,:,ch)-residual*sp.color(ch);
                end

                if max(eig(squeeze(sp.cov))) > cov_threshold
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


       
        function [shape, sp, goodness] = fit_gaussian(obj, bpatch, colorvec, absolute_position)

            norm = [10./obj.szext(1: 3), ... % loc
                    1, ... % color scale
                    1, ... % noise scale
                    ones(1, obj.szext(4)), ... % color
                    ones(1, obj.szext(4)), ... % noise
                    ];

            init_eig = 2.5*sort(obj.fsize);
            init_bound = init_eig-1;

            nonlcon = @(x) AutoDetect.constrain_eigenvalues(x, ones(1,8), init_eig, init_bound);

            colors = reshape(bpatch, [numel(bpatch)/size(bpatch, 4), size(bpatch, 4)]);
            colors(colors == 0) = eps;
            noisevec = prctile(colors, 10);

            colorvec(colorvec < 0) = eps;

            color_level = sum(colorvec(:));
            noise_level = sum(noisevec(:));

            colorvec = colorvec/sum(colorvec(:));
            noisevec = noisevec/sum(noisevec(:));


            x0 =   double([obj.fsize+1, ... % loc
                log(color_level), ... % color scale
                noise_level, ... % noise scale
                colorvec, ... % color
                noisevec, ... % noise
                ]).*norm;


            if isnan(x0(4))
                x0(4) = 2;
            end

            bounds = [30*obj.fsize./(6*obj.szext(1:3)), ... % loc
                    2, ... % color scale
                    0.2, ... % noise scale
                    0.1*ones(1, obj.szext(4)), ... % color
                    0.5*ones(1, obj.szext(4)), ... % noise
                    ];

            x0(6: 5+2*obj.szext(4)) = max(0, x0(6: 5+2*obj.szext(4)));

            lb = x0 - bounds;
            ub = x0 + bounds;

            lb(6: 5+2*obj.szext(4)) = max(0, lb(6: 5+2*obj.szext(4)));
            lb(4) = max(1, lb(4));

            lb(5) = -1;
            ub(5) = 1;

            A_eq = [zeros(1, 5), ones(1, obj.szext(4)), zeros(1, obj.szext(4)); ...
                   [zeros(1, 5), zeros(1, obj.szext(4)), ones(1, obj.szext(4))]];
            b_eq = [1; 1];


            f = @(x) AutoDetect.fit_gaussian_fixed_cov(x', bpatch, norm, obj.trunc, diag(init_eig(end:-1:1)));
            res_fixed_cov = fmincon(f, x0, [], [], A_eq, b_eq, lb, ub, [], obj.options);


            f = @(x) AutoDetect.fit_gaussian_cov(x', bpatch, ones(1,8), obj.trunc, res_fixed_cov(1:3)'./norm(1:3)', res_fixed_cov(6:5+obj.szext(4)), res_fixed_cov(6+obj.szext(4):5+2*obj.szext(4)));
            res = fmincon(f, [sqrt([init_eig(3),0,init_eig(2),0,0,init_eig(1)]),res_fixed_cov(4:5)], [], [], [], [], [eps,eps,eps,eps,eps,eps,1,-1], [10,10,10,10,10,10,100,1], nonlcon, obj.options);

            [shape, rec, tr, cov, col] = AutoDetect.get_gaussian_cov(res, [2*obj.fsize+1, obj.szext(4)], ones(1,8), obj.trunc, res_fixed_cov(1:3)./norm(1:3), res_fixed_cov(6:5+obj.szext(4)));
            bas = res_fixed_cov(6+obj.szext(4):5+2*obj.szext(4))/res(end);

%             subplot(1,3,1)
%             image(squeeze(max(rec(:,:,:,[4,3,1]), [], 3))/20)
%             subplot(1,3,2)
%             image(squeeze(max(bpatch(:,:,:,[4,3,1]), [], 3))/20)
%             subplot(1,3,3)
%             image(squeeze(max(bpatch(:,:,:,[4,3,1]) - rec(:,:,:,[4,3,1]), [], 3))/20)
%             drawnow

            reccoef = corrcoef(rec(:), bpatch(:));
            goodness = reccoef(1,2);

            relative_position = (res_fixed_cov(1:3)-x0(1:3))./norm(1:3);

            sp.mean = relative_position+absolute_position;
            sp.cov(1,:,:) = cov;
            sp.color = col;
            sp.baseline = bas;
        end
        
   end
   
end
