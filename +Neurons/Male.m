classdef Male
    %MALE Male related data.

    methods (Static)
        %% Utility functions.
        function is_neuron = isNeuron(name)
            %ISNEURON Is this a neuron name?
            is_neuron = ...
                ~isempty(find(strcmp(name, Neurons.Male.getNeurons()),1));
        end

        function is_cell = isCell(name)
            %ISCELL Is this a cell name?

            % Is this AWC?
            %is_cell = Neurons.Male.isAWC(name);
            is_cell = false;

            % Is this another neuron's name?
            if ~is_cell
                is_cell = ...
                    ~isempty(find(strcmp(name, Neurons.Male.getNeurons()),1));
            end

            % Is this a non-neuronal name?
            if ~is_cell
                is_cell = ...
                    ~isempty(find(strcmp(name, Neurons.Male.non_neuronal_cells),1));
            end
        end

        function is_AWC = isAWC(name)
            %ISAWC Is this neuron's name an acceptable form of AWC?

            % Does the neuron's name start with AWC?
            is_AWC = false;
            if length(name) > 3 && strcmp(name(1:3), 'AWC')
                if name(4) == 'L' || name(4) == 'R'

                    % AWCL/R.
                    if length(name) == 4
                        is_AWC = true;

                    % AWCL/R-ON/OFF
                    elseif strcmp(name(5:end), '-OFF') || ...
                           strcmp(name(5:end), '-ON')
                    end
                end
            end
        end

        function [name, LR] = stripLR(name)
            %STRIPLR Strip the L/R info from the neuron's name.
            
            % Is there a name?
            LR = [];
            if isempty(name)
                return;
            end
            
            % Is the neuron a left right neuron?
            name = upper(name);
            if (name(end) ~= 'L' && name(end) ~= 'R') || ...
                    any(strcmp(name, Neurons.Male.non_LR_neurons))
                return;
            end
            
            % Strip the L/R info from the neuron's name.
            LR = name(end);
            name = name(1:(end-1));
        end
        
        function name = flipLR(name)
            %FLIPLR Flip the L/R info in the neuron's name.
            
            % Is there a name?
            if isempty(name)
                return;
            end
            
            % Is the neuron a left right neuron?
            name = upper(name);
            if (name(end) ~= 'L' && name(end) ~= 'R') || ...
                    any(strcmp(name, Neurons.Male.non_LR_neurons))
                return;
            end
            
            % Flip the L/R info in the neuron's name.
            switch name(end)
                case 'L'
                    name(end) = 'R';
                case 'R'
                    name(end) = 'L';
            end
        end
        
        function [name, is_on] = stripOnOff(name)
            %STRIPONOFF Strip the ON/OFF info from the neuron's name.
            is_on = nan;
            if length(name) > 4 && strcmp(name(1:3), 'AWC')

                % Is the neuron ON or OFF?
                switch name(5:end)
                    case '-OFF'
                        is_on = false;
                    case '-ON'
                        is_on = true;
                end

                % Strip the name of ON/OF.
                name = name(1:4);
            end
        end

        function names = neuronStartsWith(str)
            %NEURONSTARTSWITH Which neurons start with this string?
            names = [];
            neurons = Neurons.Male.getNeurons();
            i = startsWith(neurons, str);
            if ~isempty(i)
                names = neurons(i);
            end
        end

        
        %% Neuron color data.
        function color = getNeuronColor(name)
            %GETNEURONCOLOR Get the neuron's NeuroPAL RGB color.
            color = [];
            [names, colors] = Neurons.Male.getColors();
            i = find(strcmp(names, name), 1);
            if ~isempty(i)
                color = colors{i};
            end
        end

        
        function name = getNeuroPALColorClass(name)
            %GETNEUROPALCOLORCLASS Get the NeuroPAL-limited color class for
            % the neuron. Most left/right and dorsal/ventral neurons share
            % equivalent NeuroPAL colors and thus have the same color class.
            if ~Neurons.Male.isNeuron(name)
                name = [];
                return;
            end
            
            % Is the neuron deterministically asymetrically colored?
            if any(strcmp(name, Neurons.Male.asym_deterministic_LR_neurons))
                return;
            end
            
            % Strip the neuron of its L/R designation.
            name = Neurons.Male.stripLR(name);
            
            % Does the neuron have 2, 4, or 6 fold symmetry?
            if name(end) == 'D' || name(end) == 'V'
                nameNoDV = name(1:(end-1));
                if any(strcmp(nameNoDV, Neurons.Male.sym_2_fold_DV_neurons)) || ...
                        any(strcmp(nameNoDV, Neurons.Male.sym_4_fold_DV_neurons)) || ...
                        any(strcmp(nameNoDV, Neurons.Male.sym_6_fold_DV_neurons))
                    name = [nameNoDV 'DV'];
                end
            end
        end
        

        %% Neuron data.
        function [names, colors] = getColors()
            %GETCOLORS Get a list of NeuroPAL neuron colors.
            %   names = neuron names
            %   colors = corresponding (R,G,B) color vector
            persistent names_data;
            persistent colors_data;
            if isempty(names_data)
                load('NeuroPAL_herm_data.mat', 'colors');
                names_data = {colors.name}';
                colors_data = {colors.RGB}';
            end
            names = names_data;
            colors = colors_data;
        end

        function neurons = getNeurons()
            %GETNEURONS Get a list of neurons.
            persistent neurons_data;
            if isempty(neurons_data)
                load('male_data.mat', 'neurons');
                neurons_data = neurons;
            end
            neurons = neurons_data;
        end

        function neurons = getNeuronsOrderedByGanglia()
            %GETNEURONSORDEREDBYGANGLIA Get a list of neurons ordered by
            % ganglia combining left & right ganglia.
            persistent neurons_data;
            if isempty(neurons_data)
                load('male_data.mat', 'neuronsOrderedByGanglia');
                neurons_data = neuronsOrderedByGanglia;
            end
            neurons = neurons_data;
        end
        
        function neurons = getNeuronsOrderedByGangliaLR()
            %GETNEURONSORDEREDBYGANGLIALR Get a list of neurons ordered by
            % ganglia separating left & right ganglia.
            persistent neurons_data;
            if isempty(neurons_data)
                load('male_data.mat', 'neuronsOrderedByGangliaLR');
                neurons_data = neuronsOrderedByGangliaLR;
            end
            neurons = neurons_data;
        end
        
        function classes = getClasses()
            %GETCLASSES Get a list of neuron classes.
            persistent classes_data;
            if isempty(classes_data)
                load('male_data.mat', 'classes');
                classes_data = classes;
            end
            classes = classes_data;
        end
        
        function neuron_class = getNeuronClass(neuron)
            % GETNEURONCLASS Get the class for this neuron.
            
            % Load the data.
            persistent neurons2classes_data;
            if isempty(neurons2classes_data)
                load('male_data.mat', 'neurons2classes');
                neurons2classes_data = neurons2classes;
            end
            
            % Find the neuron.
            neuron_class = [];
            neuron_i = strcmpi(neuron, Neurons.Male.getNeurons());
            if all(neuron_i == 0)
                return;
            end
            
            % Find the neuron's class.
            class_i = neurons2classes_data(neuron_i);
            classes = Neurons.Male.getClasses();
            neuron_class = classes{class_i};
        end
        
        function neuron = getNeuroPALName(neuron)
            % GETNEUROPALNAME Get the NeuroPAL-limited name for the neuron.
            % Left/right neurons, within the same ganglia, cannot be
            % distinguished. Therefore,they have a degenerate neuron name.
            % For example, RIGL & RIGR have a NeuroPAL-limited ID of RIG.
            
            % Is this a neuron?
            if ~Neurons.Male.isNeuron(neuron)
                neuron = [];
                return;
            end
            
            % Translate the neuron to its NeuroPAL-limited name.
            switch neuron
                % Anterior ganglion.
                %case {'IL1DL', 'IL1DR'}
                %    neuron = 'IL1D';
                %case {'IL1VL', 'IL1VR'}
                %    neuron = 'IL1V';
                %case {'IL2DL', 'IL2DR'}
                %    neuron = 'IL2D';
                %case {'IL2VL', 'IL2VR'}
                %    neuron = 'IL2V';
                
                % Ventral ganglion.
                case {'AIAL', 'AIAR'}
                    neuron = 'AIA';
                case {'RMFL', 'RMFR'}
                    neuron = 'RMF';
                case {'RMHL', 'RMHR'}
                    neuron = 'RMH';
                %case {'AIML', 'AIMR'}
                %    neuron = 'AIM';
                %case {'AIYL', 'AIYR'}
                %    neuron = 'AIY';
                %case {'AVKL', 'AVKR'}
                %    neuron = 'AVK';
                
                % Retro-vesicular ganglion.
                case {'AVFL', 'AVFR'}
                    neuron = 'AVF';
                case {'RIFL', 'RIFR'}
                    neuron = 'RIF';
                case {'RIGL', 'RIGR'}
                    neuron = 'RIG';
                case {'SABVL', 'SABVR'}
                    neuron = 'SABV';
                %case {'VD1', 'VD2'}
                    %neuron = 'VD1/2';
                
                % Pre-anal ganglion.
                case {'PVPL', 'PVPR'}
                    neuron = 'PVP';
                %case {'VD12', 'VD13'}
                    %neuron = 'VD12/13';
            end
        end
        

        %% Neuron counts.
        function num_neurons = numNeurons(body)
            %NUMNEURONS Get the number of neurons for this body part.
            num_neurons = [];
            switch lower(strtrim(body))
                case 'whole worm'
                    num_neurons = 300;
                case 'head'
                    num_neurons = Neurons.Male.numHeadNeurons();
                case 'midbody'
                    num_neurons = Neurons.Male.numMidbodyNeurons();
                case 'anterior midbody'
                    num_neurons = Neurons.Male.numAnteriorMidbodyNeurons();
                case 'central midbody'
                    num_neurons = Neurons.Male.numCentralMidbodyNeurons();
                case 'posterior midbody'
                    num_neurons = Neurons.Male.numPosteriorMidbodyNeurons();
                case 'tail'
                    num_neurons = Neurons.Male.numTailNeurons();
            end
        end
        
        function num_neurons = numHeadNeurons()
            %NUMHEADNEURONS Get the number of head neurons.
            import Neurons.*;
            persistent num;
            if isempty(num)
                num = length(Male.getAnteriorPharynx()) + ...
                    length(Male.getPosteriorPharynx()) + ...
                    length(Male.getLeftAnteriorGanglion()) + ...
                    length(Male.getRightAnteriorGanglion()) + ...
                    length(Male.getDorsalGanglion()) + ...
                    length(Male.getLeftLateralGanglion()) + ...
                    length(Male.getRightLateralGanglion()) + ...
                    length(Male.getVentralGanglion()) + ...
                    length(Male.getRetroVesicularGanglion());
            end
            num_neurons = num;
        end
        
        function num_neurons = numMidbodyNeurons()
            %NUMMIDBODYNEURONS Get the number of midbody neurons.
            import Neurons.*;
            persistent num;
            if isempty(num)
                num = length(Male.getAnteriorMidbody()) + ...
                    length(Male.getCentralMidbody()) + ...
                    length(Male.getPosteriorMidbody()) + ...
                    length(Male.getVentralNerveCord());
            end
            num_neurons = num;
        end
        
        function num_neurons = numAnteriorMidbodyNeurons()
            %NUMANTERIORMIDBODYNEURONS Get the number of anterior midbody
            % neurons.
            import Neurons.*;
            persistent num;
            if isempty(num)
                num = length(Male.getAnteriorMidbody()) + ...
                    round(length(Male.getVentralNerveCord()) / 3);
            end
            num_neurons = num;
        end
        
        function num_neurons = numCentralMidbodyNeurons()
            %NUMCENTRALMIDBODYNEURONS Get the number of central midbody
            % neurons.
            import Neurons.*;
            persistent num;
            if isempty(num)
                num = length(Male.getCentralMidbody()) + ...
                    round(length(Male.getVentralNerveCord()) / 3);
            end
            num_neurons = num;
        end
        
        function num_neurons = numPosteriorMidbodyNeurons()
            %NUMPOSTERIORMIDBODYNEURONS Get the number of posterior midbody
            % neurons.
            import Neurons.*;
            persistent num;
            if isempty(num)
                num = length(Male.getPosteriorMidbody()) + ...
                    round(length(Male.getVentralNerveCord()) / 3);
            end
            num_neurons = num;
        end
        
        function num_neurons = numTailNeurons()
            %NUMTAILNEURONS Get the number of tail neurons.
            import Neurons.*;
            persistent num;
            if isempty(num)
                num = length(Male.getPreAnalGanglion()) + ...
                    length(Male.getDorsoRectalGanglion()) + ...
                    length(Male.getLeftLumbarGanglion()) + ...
                    length(Male.getRightLumbarGanglion()) + ...
                    length(Male.getLeftRays()) + ...
                    length(Male.getRightRays()) + ...
                    length(Male.getLeftCloacalGanglion()) + ...
                    length(Male.getRightCloacalGanglion());
            end
            num_neurons = num;
        end

        
        %% Ganglia data.
        function ganglia = getGanglia()
            %GETGANGLIA Get a list of ganglia info.
            %   A stuct array where:
            %   name = ganglion name
            %   neurons = ganglion neurons
            persistent ganglia_data;
            if isempty(ganglia_data)
                load('male_data.mat', 'ganglia');
                ganglia_data = ganglia;
            end
            ganglia = ganglia_data;
        end

        function names = getGanglionNames()
            %GETGANGLIONNAMES Get a list of ganglion names.
            persistent ganglion_names;
            if isempty(ganglion_names)
                ganglia = Neurons.Male.getGanglia();
                ganglion_names = {ganglia.name}';
            end
            names = ganglion_names;
        end


        %% Pharyngeal ganglia data.
        function names = getAnteriorPharynx()
            %GETANTERIORPHARYNX A list of anterior pharynx neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.Male.getGanglia();
                i = find(strcmp(Neurons.Male.getGanglionNames(), ...
                    'Anterior Pharyngeal Bulb'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end

        function names = getPosteriorPharynx()
            %GETPOSTERIORPHARYNX A list of posterior pharynx neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.Male.getGanglia();
                i = find(strcmp(Neurons.Male.getGanglionNames(), ...
                    'Posterior Pharyngeal Bulb'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end


        %% Head ganglia data.
       function names = getLeftAnteriorGanglion()
            %GETLEFTANTERIORGANGLION A list of left anterior ganglion neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.Male.getGanglia();
                i = find(strcmp(Neurons.Male.getGanglionNames(), ...
                    'Anterior Ganglion (Left)'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end

        function names = getRightAnteriorGanglion()
            %GETRIGHTANTERIORGANGLION A list of right anterior ganglion neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.Male.getGanglia();
                i = find(strcmp(Neurons.Male.getGanglionNames(), ...
                    'Anterior Ganglion (Right)'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end

        function names = getDorsalGanglion()
            %GETDORSALGANGLION A list of dorsal ganglion neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.Male.getGanglia();
                i = find(strcmp(Neurons.Male.getGanglionNames(), ...
                    'Dorsal Ganglion'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end

        function names = getLeftLateralGanglion()
            %GETLEFTLATERALGANGLION A list of left lateral ganglion neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.Male.getGanglia();
                i = find(strcmp(Neurons.Male.getGanglionNames(), ...
                    'Lateral Ganglion (Left)'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end

        function names = getRightLateralGanglion()
            %GETRIGHTLATERALGANGLION A list of right lateral ganglion neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.Male.getGanglia();
                i = find(strcmp(Neurons.Male.getGanglionNames(), ...
                    'Lateral Ganglion (Right)'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end

        function names = getVentralGanglion()
            %GETVENTRALGANGLION A list of ventral ganglion neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.Male.getGanglia();
                i = find(strcmp(Neurons.Male.getGanglionNames(), ...
                    'Ventral Ganglion'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end
        function names = getRetroVesicularGanglion()
            %GETRETROVESICULARGANGLION A list of retro-vesicular ganglion neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.Male.getGanglia();
                i = find(strcmp(Neurons.Male.getGanglionNames(), ...
                    'Retro-Vesicular Ganglion'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end


        %% Midbody ganglia data.
        function names = getAnteriorMidbody()
            %GETANTERIORMIDBODY A list of anterior midbody neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.Male.getGanglia();
                i = find(strcmp(Neurons.Male.getGanglionNames(), ...
                    'Anterior Midbody'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end

        function names = getCentralMidbody()
            %GETCENTRALRMIDBODY A list of central midbody neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.Male.getGanglia();
                i = find(strcmp(Neurons.Male.getGanglionNames(), ...
                    'Central Midbody'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end

        function names = getPosteriorMidbody()
            %GETPOSTERIORMIDBODY A list of posterior midbody neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.Male.getGanglia();
                i = find(strcmp(Neurons.Male.getGanglionNames(), ...
                    'Posterior Midbody'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end

        function names = getVentralNerveCord()
            %GETVENTRALNERVECORD A list of ventral nerve cord neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.Male.getGanglia();
                i = find(strcmp(Neurons.Male.getGanglionNames(), ...
                    'Ventral Nerve Cord'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end


        %% Tail ganglia data.
        function names = getPreAnalGanglion()
            %GETPREANALGANGLION A list of pre-anal ganglion neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.Male.getGanglia();
                i = find(strcmp(Neurons.Male.getGanglionNames(), ...
                    'Pre-Anal Ganglion'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end

        function names = getDorsoRectalGanglion()
            %GETDORSORECTALGANGLION A list of dorso-rectal ganglion neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.Male.getGanglia();
                i = find(strcmp(Neurons.Male.getGanglionNames(), ...
                    'Dorso-Rectal Ganglion'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end

        function names = getLeftLumbarGanglion()
            %GETLEFTLUMBARGANGLION A list of left lumbar ganglion neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.Male.getGanglia();
                i = find(strcmp(Neurons.Male.getGanglionNames(), ...
                    'Lumbar Ganglion (Left)'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end

        function names = getRightLumbarGanglion()
            %GETRIGHTLUMBARGANGLION A list of right lumbar ganglion neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.Male.getGanglia();
                i = find(strcmp(Neurons.Male.getGanglionNames(), ...
                    'Lumbar Ganglion (Right)'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end
        
        function names = getLeftRays()
            %GETLEFTRAYS A list of left ray neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.Male.getGanglia();
                i = find(strcmp(Neurons.Male.getGanglionNames(), ...
                    'Rays (Left)'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end
        
        function names = getRightRays()
            %GETRIGHTRAYS A list of right ray neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.Male.getGanglia();
                i = find(strcmp(Neurons.Male.getGanglionNames(), ...
                    'Rays (Right)'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end
        
        function names = getLeftCloacalGanglion()
            %GETLEFTCLOACALGANGLION A list of left cloacal ganglion neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.Male.getGanglia();
                i = find(strcmp(Neurons.Male.getGanglionNames(), ...
                    'Cloacal Ganglion (Left)'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end
        
        function names = getRightCloacalGanglion()
            %GETRIGHTCLOACALGANGLION A list of right cloacal ganglion neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.Male.getGanglia();
                i = find(strcmp(Neurons.Male.getGanglionNames(), ...
                    'Cloacal Ganglion (Right)'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end
    end


    %% Constant properties.
    properties (Constant)
        
        % Non-neuronal cells.
        non_neuronal_cells = { ...
            'AMSOL'
            'AMSOR'
            'HMC'
            'PHSO1L'
            'PHSO1R'
            };
        
        % Non-left/right neurons.
        non_LR_neurons = {
            'ADL'
            'AQR'
            'AVL'
            'OLL'
            'PQR'
            'PVR'
            'RIR'
            'SDQL'
            'SDQR'
            };
        
        % Asymmetrically deterministic neurons.
        asym_deterministic_LR_neurons = {
            'ASEL'
            'ASER'
            };
        
        % Asymmetrically stochastic neurons.
        asym_stochastic_LR_neurons = {
            'AWCL'
            'AWCR'
            };
        
        % 2-fold symmetric neurons.
        sym_2_fold_DV_neurons = {
            'RME'
            };
        
        % 4-fold symmetric neurons.
        % Note: the SIBs & SMDs are excluded because their coloring is not
        % always symmetric for the dorsal and ventral cells.
        sym_4_fold_DV_neurons = {
            'CEP'
            'OLQ'
            'SAA'
            'SIA'
            'SMB'
            'URA'
            'URY'
            };

        % 6-fold symmetric neurons.
        sym_6_fold_DV_neurons = {
            'IL1'
            'IL2'
            'RMD'
            };
    end
end
