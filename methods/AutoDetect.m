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
        
            sz = size(volume);
            reconstruction = image.get_3d_shape(sz,mp_params.hnsz);

            vec = @(x) x(:);


            % reconstruction based
            eval = [];
            eval.recon.mse = nanmean(vec((reconstruction - volume).^2));
            correlation = corrcoef(vec(reconstruction), vec(volume));
            eval.recon.cor = correlation(1,2);

            % shape based
            segmentation = image.get_3d_segments(sz,mp_params.hnsz);

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
                    [~,cols] = find(munkres(distances));
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
                min(eig(squeeze(sp.cov).*obj.scale))
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
            eig_ub = [4,4,4];
            

            nonlcon = @(x) AutoDetect.constrain_eigenvalues(x, eig_lb, eig_ub);

            colors = reshape(bpatch, [numel(bpatch)/size(bpatch, 4), size(bpatch, 4)]);
            colors(colors == 0) = eps;
            noisevec = prctile(colors, 10);

            colorvec(colorvec < 0) = eps;

            color_level = sum(colorvec(:));
            noise_level = sum(noisevec(:));

            colorvec = colorvec/sum(colorvec(:));
            noisevec = noisevec/sum(noisevec(:));

            
            x0 =   double([(obj.fsize+1).*obj.scale, ... % loc
                sqrt([3,0,0,3,0,3]), ... % covariance parameters
                1.5*log(color_level), ... % color scale
                noise_level, ... % noise scale
                colorvec, ... % color
                noisevec, ... % noise
                0]);
            
            lb = [x0(1:3)-0.3*ones(1,3), zeros(1,6), 0, -0.1, zeros(1,2*obj.szext(4)), 0];
            ub = [x0(1:3)+0.3*ones(1,3), 100*ones(1,6), 20, 10, ones(1,2*obj.szext(4)), 1];
            
            
            A_eq = [zeros(1, 11), ones(1, obj.szext(4)), zeros(1, obj.szext(4)), 0; ...
                   [zeros(1, 11), zeros(1, obj.szext(4)), ones(1, obj.szext(4)), 0]];
            b_eq = [1; 1];

            
            f = @(x) AutoDetect.gaussian_cost(x', bpatch, norm);
            x_hat = fmincon(f, x0, [], [], A_eq, b_eq, lb, ub, nonlcon, obj.options);

            [shape, rec, tr, cov, col] = AutoDetect.get_gaussian(x_hat, [2*obj.fsize+1, obj.szext(4)], norm);
            bas = x_hat(12+obj.szext(4):11+2*obj.szext(4))*x_hat(11);

            relative_position = (x_hat(1:3)-x0(1:3))./norm(1:3);

            sp.mean = relative_position+absolute_position;
            sp.cov(1,:,:) = cov;
            sp.color = col;
            sp.baseline = bas;
            sp.truncation = x_hat(end)/norm(end);
        end
        
        
        
                function [sp, output, neuron_collections, neuron_collections_fixed] = detect(Ds, Nneuron, file_name, fixed_ball, refine_shape_label, fitted_save_less)
        % running greedy algorithm -- matching pursuit (MP). 
        % Input: 
        %         Ds: z-scored raw data
        %         Nneuron: number of detections allowed
        %         file_name: to save MP results in super pixel format ("sp")
        %                    e.g. '38_YAaSP'
        %         fixed_ball (optional): neuron shape. default 3D gaussian
        %                                with coravance equal to [32,0,0;0,32,0;0,0,4]; 
        %         refine_shape_label (optional): Boolean. True (default): refine shape
        %                           after each MP iteration; False: not refine.
        %         fitted_save_less (optional): Boolean. True (default):
        %                                      not save 3D reconstruction
        %                                      for each iteration. False:
        %                                      save. 
        % Output: 
        %       sp: super pixel sp
        %       output:  x, y, z, r, g, b, w   
        %       neuron_collections:  collections of all neunons with fitted shape
        %       neuron_collections_fixed: collections of all neurons with fixed shape
        %       save "38_YAaSP_sp.mat" format for gui. 
        % Usage: 
        %   [sp, output, neuron_collections, neuron_collections_fixed] = detect(data, 10, '38_YAaSP'); 


            if nargin < 6
                fitted_save_less = true;
            end
            if nargin < 5
                refine_shape_label = true; 
            end
            if nargin < 4
                box= 10; 
                boxz = 4; 
                [Xmesh,Ymesh,Zmesh]=meshgrid(1:box+1,1:box+1,1:boxz+1); 
                sigma0 = [32,0,0;0,32,0;0,0,4]; 
                density = @(X) mvnpdf(X, [box/2+1,box/2+1,boxz/2+1], sigma0);  
                XX = [Ymesh(:), Xmesh(:), Zmesh(:)]; 
                fixed_ball = reshape(density(XX),[box+1,box+1,boxz+1]);
            end

            output =[];
            [xsize, ysize, zsize, ncolor] = size(Ds);
            box = size(fixed_ball,1)-1;
            boxz = size(fixed_ball,3)-1;

            
            w = [1,1,1,1/3] ; 

            % initial x0 (r, g, b)
            % should be adjusted according to data intensity level.
            % otherwise, optimization may fail. 
         
            x0 = [0, 4.3160e+03, 0, 4.3287e+03, 0, 4.3287e+03, 0, 4.1733e+03]; % normalize data
  
            rhos = zeros(size(Ds)); 

            for c = 1:ncolor
                rhos(:,:,:,c) = imfilter(Ds(:,:,:,c), fixed_ball, 'same'); 
            end

            % fixed ball window
            wz = (size(fixed_ball,3)-1)/2; 
            wx = (size(fixed_ball,1)-1)/2;
            wy=  (size(fixed_ball,2)-1)/2;


            neuron_collections = [];
            neuron_collections_fixed = []; 

            fixed_ball_l2_normalized = false % if fixed_ball is l2 normalized, sum(fixed_ball(:))==1, intensity can be calculated directly

            for i = 1:Nneuron
                i
                fixed = true;     
                rho = sum(rhos.^2, 4); %rho1.^2 + rho2.^2 + rho4.^2; 
                sign = sum(rhos, 4)>0;  % avoid negative
                [val, idx] = max(rho(:).*sign(:));
                [I_row, I_col, I_z] = ind2sub(size(rho),idx); 
                center = [I_row, I_col, I_z];
    
                % not fixed ball
                if refine_shape_label          
                    fixed = false; 
                    try                                      
                        fitted = fit_gaussian_colors4_box(Ds, center, w, x0, box, boxz);                
                    catch
                        fixed = true; 
                    end
  
                    % check if fitted value is correct or not
                    if ~fixed
                        if fitted.Sigma(1,1)>75 | fitted.Sigma(2,2)>75 | fitted.Sigma(3,3)> 7.5
                            disp(['fitted Sigma too large. i=',num2str(i)])
                            fixed = true; 
                        end
            
                        if fitted.Sigma(3,3)<1 
                            disp(['smallest threshold change to 1, not 2'])
                            disp(['fitted Sigma too small. i=',num2str(i)])
                            fixed = true;                         
                        end
  
                        fitted_x = fitted.box_size(1);
                        fitted_y = fitted.box_size(2);
                        fitted_z = fitted.box_size(3);

                        % fit center to close to edge
                        if fitted.box_center(1)>fitted_x-1 | fitted.box_center(1)<1
                            disp(['wrong box_center x. i=',num2str(i)])
                            fixed = true; 
                        end    
                        if fitted.box_center(2)> fitted_y-1 | fitted.box_center(2)<1
                             disp(['wrong box_center y. i=',num2str(i)])
                             fixed = true; 
                        end
             
                        if fitted.box_center(3)> fitted_z-1 | fitted.box_center(3)<1
                             disp(['wrong box_center z. i=',num2str(i)])
                             fixed = true; 
                        end
            
                    end % end refine shape
  
                 % varying-ball MP
                 if ~fixed     
                     fitted.box_vals(fitted.box_vals<0.1) = 0; %new
                     Ds  = Ds  - fitted.box_vals; 
                     for cc = 1:ncolor
                         this = fitted.box_vals(:,:,:,cc); 
                         rho = rhos(:,:,:,cc); 
                         tmp = imfilter(this, fixed_ball, 'same'); 
                         rho  = rho - tmp;
                         rhos(:,:,:,cc) =  rho; 
                     end
                     xis = cell2mat(fitted.scale); 
                     fitted.i = i; 
                     if fitted_save_less
                           fitted.box_vals = [];
                     end
                     neuron_collections = [neuron_collections fitted]; 
  
                 end
        
                end % end varying-ball MP
  
                % fixed ball MP
                if fixed
                    fitted.fixed = true; 
                    % fixed ball (need l2 normalization to make intensity correct)
                    if fixed_ball_l2_normalized %intensity
                        xis = zeros([1, ncolor]); 
                        fitted = {}; 
                        fitted.fixed = true; 
                        fitted.this = zeros(size(Ds)); 

                        for cc = 1:ncolor
       
                            D =  Ds(:,:,:,cc);  
                            rho = rhos(:,:,:,cc);%
                    
                            D_extend = zeros(size(D)+[2*wx, 2*wy, 2*wz]); 
                            D_extend(wx+1:size(D,1)+wx, wy+1:size(D,2)+wy, wz+1:size(D,3)+wz) = D; 
                            temp = D_extend(I_row:I_row+2*wx, I_col:I_col+2*wy, I_z:I_z+2*wz) .* fixed_ball(size(fixed_ball,1):-1:1,size(fixed_ball,2):-1:1,size(fixed_ball,3):-1:1); 
  
                            xi = sum(temp(:)); 
                            xis(cc) = xi;
  
                            this_extend = zeros(size(D)+[2*wx, 2*wy, 2*wz]);  
                            this_extend(I_row:I_row+2*wx, I_col:I_col+2*wy, I_z:I_z+2*wz) =xi*fixed_ball(size(fixed_ball,1):-1:1,size(fixed_ball,2):-1:1,size(fixed_ball,3):-1:1); 
                            this = this_extend(wx+1:size(D,1)+wx, wy+1:size(D,2)+wy, wz+1:size(D,3)+wz);
  
                            Ds(:,:,:,cc) = D - this;  
                            rhos(:,:,:,cc) =  rho - imfilter(this, fixed_ball, 'same');
 
                        end

                        fitted.this(:,:,:,cc) = this; 
                        fitted.scale = xis;            
                    else       
                
                        sigma0 = [32,0,0;0,32,0;0,0,4];
                        if ncolor==3
                            fitted = fit_intensity_colors3(Ds, center, sigma0, x0); 
                        elseif ncolor==4
                            fitted = fit_intensity_colors4(Ds, center, sigma0, w, x0); 
                        end
                        Ds  = Ds  - fitted.box_vals; 
            
  
                        for cc = 1:ncolor
                            this = fitted.box_vals(:,:,:,cc); 
                            rho = rhos(:,:,:,cc);      
                            rho  = rho - imfilter(this, fixed_ball, 'same');
                            rhos(:,:,:,cc) =  rho; 
                        end
                        xis = cell2mat(fitted.scale); 
  
                    end
                    fitted.i = i;
                    if fitted_save_less
                        fitted.box_vals = [];
                    end
                    neuron_collections_fixed = [neuron_collections_fixed fitted]; 
     
                end  % end fixed ball MP  
                output = [output;[fitted.center xis]];      
            end
       
             % save to super pixel format
            [mp_params, sp] = get_super_pixel(neuron_collections, neuron_collections_fixed);    
            save([char(file_name),'_sp.mat'], 'mp_params', 'sp'); 
      
        end
 
        
        function fitted = fit_gaussian_colors4_box(Ds, center, w, user_x0, box, boxz)
        % refine neuron's shape by fitting a gaussian, and extract center xyz, color, and background, etc    
        % Ruoxi Sun. 

            Xc = center(1);
            Yc = center(2);
            Zc = center(3);
    
            data1 = Ds(:,:,:,1); 
            data2 = Ds(:,:,:,2); 
            data3 = Ds(:,:,:,3); 
            data4 = Ds(:,:,:,4); 

            % extend border box/2, so that we can extract the window when the center is
            % at edge
            [xsize,ysize,zsize, ncolor] = size(Ds); 
            raw_data_extend = zeros([xsize+box+1, ysize+box+1, zsize+boxz+1, 4]); 
            raw_data = zeros([box+1,box+1,boxz+1,4]);
            b = zeros([1,4]); 
            for j = 1:4
                temp = Ds(:,:,:,j); 
                b(j) =  quantile(temp(:),.3);
                raw_data_extend(:,:,:,j) = ones([xsize+box+1, ysize+box+1, zsize+boxz+1])*b(j); 
                raw_data_extend(box/2+1:box/2+xsize, box/2+1:box/2+ysize, boxz/2+1:boxz/2+zsize,j) = temp;
                raw_data(:,:,:,j) = raw_data_extend(Xc:Xc+box, Yc:Yc+box, Zc:Zc+boxz,j); 
            end
            raw_data1 = raw_data(:,:,:,1);
            raw_data2 = raw_data(:,:,:,2);
            raw_data3 = raw_data(:,:,:,3);
            raw_data4 = raw_data(:,:,:,4);
 
            y1 = raw_data1(:);
            y2 = raw_data2(:);
            y3 = raw_data3(:);
            y4 = raw_data4(:);

            full_x = size(Ds,1); 
            full_y = size(Ds,2); 
            full_z = size(Ds,3); 

            [Xmesh,Ymesh,Zmesh]=meshgrid(1:box+1,1:box+1,1:boxz+1);
            X = [Ymesh(:),Xmesh(:),Zmesh(:)];

            f_mvn1 = @(x) abs(x(10))+abs(x(11))*mvnpdf(X,abs([x(1),x(2),x(3)]), [x(7),x(5),x(8);x(5),x(4),x(6);x(8),x(6),x(9)]); 
            f_mvn2 = @(x) abs(x(12))+abs(x(13))*mvnpdf(X,abs([x(1),x(2),x(3)]), [x(7),x(5),x(8);x(5),x(4),x(6);x(8),x(6),x(9)]); 
            f_mvn3 = @(x) abs(x(14))+abs(x(15))*mvnpdf(X,abs([x(1),x(2),x(3)]), [x(7),x(5),x(8);x(5),x(4),x(6);x(8),x(6),x(9)]); 
            f_mvn4 = @(x) abs(x(16))+abs(x(17))*mvnpdf(X,abs([x(1),x(2),x(3)]), [x(7),x(5),x(8);x(5),x(4),x(6);x(8),x(6),x(9)]); 

            f = @(x) sum( (w(1)*(y1-f_mvn1(x))).^2+(w(2)*(y2-f_mvn2(x))).^2+(w(3)*(y3-f_mvn3(x))).^2+(w(4)*(y4-f_mvn4(x))).^2); 
            x0_shape = [size(raw_data2,1)/2,size(raw_data2,2)/2, size(raw_data2,3)/2,  32, 0, 0, 32, 0, 4];
            x0 = [x0_shape  user_x0];  

            options = optimset('MaxFunEvals',10000,'MaxIter',10000);
            x = fminsearch(f,x0,options);

           %switch order
           [Xmesh1,Ymesh1,Zmesh1]=meshgrid(1-(Yc-box/2-1):full_y-(Yc-box/2-1),1-(Xc-box/2-1):full_x-(Xc-box/2-1),1-(Zc-boxz/2-1):full_z-(Zc-boxz/2-1));
           X1 = [Ymesh1(:),Xmesh1(:),Zmesh1(:)]; 

            % already subtract background
            f1 = @(X) abs(x(11))*mvnpdf(X,abs([x(1),x(2),x(3)]), [x(7),x(5),x(8);x(5),x(4),x(6);x(8),x(6),x(9)]);  
            f2 = @(X) abs(x(13))*mvnpdf(X,abs([x(1),x(2),x(3)]), [x(7),x(5),x(8);x(5),x(4),x(6);x(8),x(6),x(9)]);  
            f3 = @(X) abs(x(15))*mvnpdf(X,abs([x(1),x(2),x(3)]), [x(7),x(5),x(8);x(5),x(4),x(6);x(8),x(6),x(9)]);  
            f4 = @(X) abs(x(17))*mvnpdf(X,abs([x(1),x(2),x(3)]), [x(7),x(5),x(8);x(5),x(4),x(6);x(8),x(6),x(9)]);  

            vals1 = f1(X1); 
            vals2 = f2(X1);
            vals3 = f3(X1);
            vals4 = f4(X1);

            mvg_box1 = reshape(vals1,[full_x,full_y,full_z]); 
            mvg_box2 = reshape(vals2,[full_x,full_y,full_z]); 
            mvg_box3 = reshape(vals3,[full_x,full_y,full_z]); 
            mvg_box4 = reshape(vals4,[full_x,full_y,full_z]); 

            box_vals = zeros([size(mvg_box4),4]); 
            box_vals(:,:,:,1) = mvg_box1;  
            box_vals(:,:,:,2) = mvg_box2; 
            box_vals(:,:,:,3) = mvg_box3; 
            box_vals(:,:,:,4) = mvg_box4; 

            fitted = {};
            fitted.fixed = false;
            fitted.box_center = abs([x(1),x(2),x(3)]); 
            fitted.center = [Xc-box/2+abs(x(1))-1, Yc-box/2+abs(x(2))-1, Zc-boxz/2+abs(x(3))-1]; 
            fitted.Sigma = [x(7),x(5),x(8);x(5),x(4),x(6);x(8),x(6),x(9)]; 
            fitted.x0 = x0; 
            fitted.x = x; 
            fitted.funcs = {f_mvn1,f_mvn2,f_mvn3,f_mvn4}; 
            fitted.box_size = size(raw_data2); 

            % too large comments out (not save) !!!!!!!!!!!!!!!!!!
            fitted.box_vals = box_vals; 

            fitted.scale = {abs(x(11)),abs(x(13)),abs(x(15)),abs(x(17))};  
            fitted.background = {abs(x(10)),abs(x(12)),abs(x(14)),abs(x(16))}; 
            fitted.box_range = {[Xc-box/2:Xc+box/2], [Yc-box/2:Yc+box/2], [Zc-boxz/2:Zc+boxz/2]}; 

        end
        
        
        function [mp_params, sp] = get_super_pixel(neuron_collections, neuron_collections_fixed)
            % from MP format to super pixel format
            % Ruoxi Sun
    
            % not fixed
            nc_colors =  transpose(reshape(cell2mat([neuron_collections.scale]),[4, length(neuron_collections)])); 
            nc_centers =  transpose(reshape([neuron_collections.center],[3, length(neuron_collections)])); 
            nc_i =  [neuron_collections.i];
            nc_fixed =  [neuron_collections.fixed];  
            nc_background = transpose(reshape(cell2mat([neuron_collections.background]), 4, length(neuron_collections)));
            nc_sigma = reshape([neuron_collections.Sigma],3,3,length(neuron_collections)); 
            nc = [transpose(nc_i), nc_centers, nc_colors, transpose(nc_fixed), nc_background];

            % fixed
            if length(neuron_collections_fixed)~=0
                fixed_nc_colors =  transpose(reshape(cell2mat([neuron_collections_fixed.scale]),[4, length(neuron_collections_fixed)])); 
                fixed_nc_centers =  transpose(reshape([neuron_collections_fixed.center],[3, length(neuron_collections_fixed)])); 
                fixed_nc_i =  [neuron_collections_fixed.i]; 
                fixed_nc_fixed =  [neuron_collections_fixed.fixed]; 
                fixed_nc_background = transpose(reshape(cell2mat([neuron_collections_fixed.background]), 4, length(neuron_collections_fixed)));
                fixed_nc_sigma = reshape([neuron_collections_fixed.Sigma],3,3,length(neuron_collections_fixed)); 
                fixed_nc = [transpose(fixed_nc_i), fixed_nc_centers, fixed_nc_colors, transpose(fixed_nc_fixed),fixed_nc_background];
            else
                fixed_nc = [];
            end
    
            % all 
            neurons = [nc; fixed_nc]; 
            sigmas = cat(3, nc_sigma, fixed_nc_sigma); 
            [~,idx] = sort(neurons(:,1)); % sort just the first column
            sortedneurons = neurons(idx,:); 
            sortedsigmas = sigmas(:,:,idx); 

            sp.mean = sortedneurons(:,2:4);
            % remove duplicate
            [C,ii] = unique(round(sp.mean),'rows', 'stable'); 
            if length(ii)~=200
                length(ii)
            end
            sp.mean = sp.mean(ii,:); 
            sp.cov = permute(sortedsigmas,[3 1 2]);
            sp.cov = sp.cov(ii,:,:); 
            sp.color = sortedneurons(ii,5:8);
            sp.baseline = sortedneurons(ii,10:13);

             mp_params.hnsz = [10,10,4]; 
             mp_params.trunc = NaN;
             mp_params.k = length(ii); 
             mp_params.bca_iter = NaN; 
        end
        
               
        
        function [reconstructs, reconstructs_fixed] =  reconstruct2(full_size, ncolor, neuron_collections, neuron_collections_fixed)
            % reconstruct 3D volumn from MP results
            % Ruoxi Sun. 

            disp(['fitted neurons = ',num2str(length(neuron_collections))]);
            disp(['fixed neurons = ',num2str(length(neuron_collections_fixed))]); 


            % fitted
            full_x = full_size(1); 
            full_y = full_size(2); 
            full_z = full_size(3); 

            reconstructs = zeros([full_x,full_y,full_z,ncolor]);
 
            for i = 1:length(neuron_collections)%size(temp_this,4) 
                i
                fitted = neuron_collections(i); 

                if size(fitted.box_vals,1)==0
                    fitted = reconstruct_boxvals(full_size,fitted); 
                end 
    
                reconstructs = reconstructs + fitted.box_vals; 
 
            end

            % fixed

            reconstructs_fixed = zeros([full_x,full_y,full_z,ncolor]);
 

            for j = 1:length(neuron_collections_fixed)%size(temp_this,4) 
                j
    
                fitted = neuron_collections_fixed(j); 
    
                if size(fitted.box_vals,1)==0
                    fitted = reconstruct_boxvals(full_size,fitted); 
                end 
    
 
                reconstructs_fixed = reconstructs_fixed + fitted.box_vals;  %(:,:,:,cc) = temp_this; %sum_this; 
 
            end
    
        end

        
        
        
        
        
        
   end
   
end
