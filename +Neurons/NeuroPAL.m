classdef NeuroPAL
    %NEUROPAL NeuroPAL and worm related data.

    methods (Static)
        %% Utility functions.
        function is_neuron = isNeuron(name)
            %ISNEURON Is this a neuron name?
            is_neuron = ...
                ~isempty(find(strcmp(name, Neurons.NeuroPAL.getNeurons()),1));
        end

        function is_cell = isCell(name)
            %ISCELL Is this a cell name?

            % Is this AWC?
            %is_cell = Neurons.NeuroPAL.isAWC(name);
            is_cell = false;

            % Is this another neuron's name?
            if ~is_cell
                is_cell = ...
                    ~isempty(find(strcmp(name, Neurons.NeuroPAL.getNeurons()),1));
            end

            % Is this a non-neuronal name?
            if ~is_cell
                is_cell = ...
                    ~isempty(find(strcmp(name, Neurons.NeuroPAL.non_neuronal_cells),1));
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
            neurons = Neurons.NeuroPAL.getNeurons();
            i = startsWith(neurons, str);
            if ~isempty(i)
                names = neurons(i);
            end
        end

        function color = getNeuronColor(name)
            %GETNEURONCOLOR Get the neuron's NeuroPAL RGB color.
            color = [];
            [names, colors] = Neurons.NeuroPAL.getColors();
            i = find(strcmp(names, name), 1);
            if ~isempty(i)
                color = colors{i};
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
                load('herm_data.mat', 'neurons');
                neurons_data = neurons;
            end
            neurons = neurons_data;
        end

        function neurons = getNeuronsOrderedByGanglia()
            %GETNEURONSORDEREDBYGANGLIA Get a list of neurons ordered by
            % ganglia combining left & right ganglia.
            persistent neurons_data;
            if isempty(neurons_data)
                load('herm_data.mat', 'neuronsOrderedByGanglia');
                neurons_data = neuronsOrderedByGanglia;
            end
            neurons = neurons_data;
        end
        
        function neurons = getNeuronsOrderedByGangliaLR()
            %GETNEURONSORDEREDBYGANGLIALR Get a list of neurons ordered by
            % ganglia separating left & right ganglia.
            persistent neurons_data;
            if isempty(neurons_data)
                load('herm_data.mat', 'neuronsOrderedByGangliaLR');
                neurons_data = neuronsOrderedByGangliaLR;
            end
            neurons = neurons_data;
        end
        
        function classes = getClasses()
            %GETCLASSES Get a list of neuron classes.
            persistent classes_data;
            if isempty(classes_data)
                load('herm_data.mat', 'classes');
                classes_data = classes;
            end
            classes = classes_data;
        end
        
        function neuron_class = getNeuronClass(neuron)
            % GETNEURONCLASS Get the class for this neuron.
            
            % Load the data.
            persistent neurons2classes_data;
            if isempty(neurons2classes_data)
                load('herm_data.mat', 'neurons2classes');
                neurons2classes_data = neurons2classes;
            end
            
            % Find the neuron.
            neuron_class = [];
            neuron_i = strcmpi(neuron, Neurons.NeuroPAL.getNeurons());
            if all(neuron_i == 0)
                return;
            end
            
            % Find the neuron's class.
            class_i = neurons2classes_data(neuron_i);
            classes = Neurons.NeuroPAL.getClasses();
            neuron_class = classes{class_i};
        end

        
        function neuron = getNeuroPALName(neuron)
            % GETNEUROPALNAME Get the NeuroPAL-limited name for the neuron.
            % Left/right neurons, within the same ganglia, cannot be
            % distinguished. Therefore,they have a degenrate neuron name.
            % For example, RIGL & RIGR have a NeuroPAL-limited ID of RIG.
            
            % Is this a neuron?
            if ~Neurons.NeuroPAL.isNeuron(neuron)
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
                    num_neurons = Neurons.NeuroPAL.numHeadNeurons();
                case 'midbody'
                    num_neurons = Neurons.NeuroPAL.numMidbodyNeurons();
                case 'anterior midbody'
                    num_neurons = Neurons.NeuroPAL.numAnteriorMidbodyNeurons();
                case 'central midbody'
                    num_neurons = Neurons.NeuroPAL.numCentralMidbodyNeurons();
                case 'posterior midbody'
                    num_neurons = Neurons.NeuroPAL.numPosteriorMidbodyNeurons();
                case 'tail'
                    num_neurons = Neurons.NeuroPAL.numTailNeurons();
            end
        end
        
        function num_neurons = numHeadNeurons()
            %NUMHEADNEURONS Get the number of head neurons.
            import Neurons.*;
            persistent num;
            if isempty(num)
                num = length(NeuroPAL.getAnteriorPharynx()) + ...
                    length(NeuroPAL.getPosteriorPharynx()) + ...
                    length(NeuroPAL.getLeftAnteriorGanglion()) + ...
                    length(NeuroPAL.getRightAnteriorGanglion()) + ...
                    length(NeuroPAL.getDorsalGanglion()) + ...
                    length(NeuroPAL.getLeftLateralGanglion()) + ...
                    length(NeuroPAL.getRightLateralGanglion()) + ...
                    length(NeuroPAL.getVentralGanglion()) + ...
                    length(NeuroPAL.getRetroVesicularGanglion());
            end
            num_neurons = num;
        end
        
        function num_neurons = numMidbodyNeurons()
            %NUMMIDBODYNEURONS Get the number of midbody neurons.
            import Neurons.*;
            persistent num;
            if isempty(num)
                num = length(NeuroPAL.getAnteriorMidbody()) + ...
                    length(NeuroPAL.getCentralMidbody()) + ...
                    length(NeuroPAL.getPosteriorMidbody()) + ...
                    length(NeuroPAL.getVentralNerveCord());
            end
            num_neurons = num;
        end
        
        function num_neurons = numAnteriorMidbodyNeurons()
            %NUMANTERIORMIDBODYNEURONS Get the number of anterior midbody
            % neurons.
            import Neurons.*;
            persistent num;
            if isempty(num)
                num = length(NeuroPAL.getAnteriorMidbody()) + ...
                    round(length(NeuroPAL.getVentralNerveCord()) / 3);
            end
            num_neurons = num;
        end
        
        function num_neurons = numCentralMidbodyNeurons()
            %NUMCENTRALMIDBODYNEURONS Get the number of central midbody
            % neurons.
            import Neurons.*;
            persistent num;
            if isempty(num)
                num = length(NeuroPAL.getCentralMidbody()) + ...
                    round(length(NeuroPAL.getVentralNerveCord()) / 3);
            end
            num_neurons = num;
        end
        
        function num_neurons = numPosteriorMidbodyNeurons()
            %NUMPOSTERIORMIDBODYNEURONS Get the number of posterior midbody
            % neurons.
            import Neurons.*;
            persistent num;
            if isempty(num)
                num = length(NeuroPAL.getPosteriorMidbody()) + ...
                    round(length(NeuroPAL.getVentralNerveCord()) / 3);
            end
            num_neurons = num;
        end
        
        function num_neurons = numTailNeurons()
            %NUMTAILNEURONS Get the number of tail neurons.
            import Neurons.*;
            persistent num;
            if isempty(num)
                num = length(NeuroPAL.getPreAnalGanglion()) + ...
                    length(NeuroPAL.getDorsoRectalGanglion()) + ...
                    length(NeuroPAL.getLeftLumbarGanglion()) + ...
                    length(NeuroPAL.getRightLumbarGanglion());
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
                load('herm_data.mat', 'ganglia');
                ganglia_data = ganglia;
            end
            ganglia = ganglia_data;
        end

        function names = getGanglionNames()
            %GETGANGLIONNAMES Get a list of ganglion names.
            persistent ganglion_names;
            if isempty(ganglion_names)
                ganglia = Neurons.NeuroPAL.getGanglia();
                ganglion_names = {ganglia.name}';
            end
            names = ganglion_names;
        end


        %% Pharyngeal ganglia data.
        function names = getAnteriorPharynx()
            %GETANTERIORPHARYNX A list of anterior pharynx neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.NeuroPAL.getGanglia();
                i = find(strcmp(Neurons.NeuroPAL.getGanglionNames(), ...
                    'Anterior Pharyngeal Bulb'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end

        function names = getPosteriorPharynx()
            %GETPOSTERIORPHARYNX A list of posterior pharynx neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.NeuroPAL.getGanglia();
                i = find(strcmp(Neurons.NeuroPAL.getGanglionNames(), ...
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
                ganglia = Neurons.NeuroPAL.getGanglia();
                i = find(strcmp(Neurons.NeuroPAL.getGanglionNames(), ...
                    'Anterior Ganglion (Left)'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end

        function names = getRightAnteriorGanglion()
            %GETRIGHTANTERIORGANGLION A list of right anterior ganglion neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.NeuroPAL.getGanglia();
                i = find(strcmp(Neurons.NeuroPAL.getGanglionNames(), ...
                    'Anterior Ganglion (Right)'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end

        function names = getDorsalGanglion()
            %GETDORSALGANGLION A list of dorsal ganglion neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.NeuroPAL.getGanglia();
                i = find(strcmp(Neurons.NeuroPAL.getGanglionNames(), ...
                    'Dorsal Ganglion'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end

        function names = getLeftLateralGanglion()
            %GETLEFTLATERALGANGLION A list of left lateral ganglion neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.NeuroPAL.getGanglia();
                i = find(strcmp(Neurons.NeuroPAL.getGanglionNames(), ...
                    'Lateral Ganglion (Left)'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end

        function names = getRightLateralGanglion()
            %GETRIGHTLATERALGANGLION A list of right lateral ganglion neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.NeuroPAL.getGanglia();
                i = find(strcmp(Neurons.NeuroPAL.getGanglionNames(), ...
                    'Lateral Ganglion (Right)'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end

        function names = getVentralGanglion()
            %GETVENTRALGANGLION A list of ventral ganglion neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.NeuroPAL.getGanglia();
                i = find(strcmp(Neurons.NeuroPAL.getGanglionNames(), ...
                    'Ventral Ganglion'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end
        function names = getRetroVesicularGanglion()
            %GETRETROVESICULARGANGLION A list of retro-vesicular ganglion neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.NeuroPAL.getGanglia();
                i = find(strcmp(Neurons.NeuroPAL.getGanglionNames(), ...
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
                ganglia = Neurons.NeuroPAL.getGanglia();
                i = find(strcmp(Neurons.NeuroPAL.getGanglionNames(), ...
                    'Anterior Midbody'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end

        function names = getCentralMidbody()
            %GETCENTRALRMIDBODY A list of central midbody neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.NeuroPAL.getGanglia();
                i = find(strcmp(Neurons.NeuroPAL.getGanglionNames(), ...
                    'Central Midbody'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end

        function names = getPosteriorMidbody()
            %GETPOSTERIORMIDBODY A list of posterior midbody neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.NeuroPAL.getGanglia();
                i = find(strcmp(Neurons.NeuroPAL.getGanglionNames(), ...
                    'Posterior Midbody'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end

        function names = getVentralNerveCord()
            %GETVENTRALNERVECORD A list of ventral nerve cord neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.NeuroPAL.getGanglia();
                i = find(strcmp(Neurons.NeuroPAL.getGanglionNames(), ...
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
                ganglia = Neurons.NeuroPAL.getGanglia();
                i = find(strcmp(Neurons.NeuroPAL.getGanglionNames(), ...
                    'Pre-Anal Ganglion'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end

        function names = getDorsoRectalGanglion()
            %GETDORSORECTALGANGLION A list of dorso-rectal ganglion neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.NeuroPAL.getGanglia();
                i = find(strcmp(Neurons.NeuroPAL.getGanglionNames(), ...
                    'Dorso-Rectal Ganglion'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end

        function names = getLeftLumbarGanglion()
            %GETLEFTLUMBARGANGLION A list of left lumbar ganglion neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.NeuroPAL.getGanglia();
                i = find(strcmp(Neurons.NeuroPAL.getGanglionNames(), ...
                    'Lumbar Ganglion (Left)'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end

        function names = getRightLumbarGanglion()
            %GETRIGHTLUMBARGANGLION A list of right lumbar ganglion neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = Neurons.NeuroPAL.getGanglia();
                i = find(strcmp(Neurons.NeuroPAL.getGanglionNames(), ...
                    'Lumbar Ganglion (Right)'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end
    end


    %% Constant properties.
    properties (Constant)
        non_neuronal_cells = { ...
            'AMSOL'
            'AMSOR'
            'HMC'
            'PHSO1L'
            'PHSO1R'
            };
    end
end
