classdef NeuroPAL
    %NEUROPAL NeuroPAL and worm related data.

    methods (Static)
        %% Utility functions.
        function is_neuron = isNeuron(name)
            %ISNEURON Is this a neuron name?
            is_neuron = ...
                ~isempty(find(strcmp(name, NeuroPAL.getNeurons()),1));
        end
                
        function name = stripOnOff(name)
            %STRIPONOFF Strip the ON/OFF info from the neuron's name.
            if length(name) > 4 && strcmp(name(1:3), 'AWC')
                name = name(1:4);
            end
        end
        
        function names = neuronStartsWith(str)
            %NEURONSTARTSWITH Which neurons start with this string?
            names = [];
            neurons = NeuroPAL.getNeurons();
            i = startsWith(neurons, str);
            if ~isempty(i)
                names = neurons(i);
            end
        end
        
        function color = getNeuronColor(name)
            %GETNEURONCOLOR Get the neuron's NeuroPAL RGB color.
            color = [];
            [names, colors] = getColors();
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
        
        function classes = getClasses()
            %GETCLASSES Get a list of neuron classes.
            persistent classes_data;
            if isempty(classes_data)
                load('herm_data.mat', 'classes');
                classes_data = classes;
            end
            classes = neurons_data;
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
                ganglia = NeuroPAL.getGanglia();
                ganglion_names = {ganglia.name}';
            end
            names = ganglion_names;
        end
        
        
        %% Pharyngeal ganglia data.
        function names = getAnteriorPharynxNeurons()
            %GETANTERIORPHARYNXNEURONS A list of anterior pharynx neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = NeuroPAL.getGanglia();
                i = find(strcmp(NeuroPAL.getGanglionNames(), ...
                    'Anterior Pharyngeal Bulb'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end
        
        function names = getPosteriorPharynxNeurons()
            %GETPOSTERIORPHARYNXNEURONS A list of posterior pharynx neurons.
            persistent neurons;
            if isempty(neurons)
                ganglia = NeuroPAL.getGanglia();
                i = find(strcmp(NeuroPAL.getGanglionNames(), ...
                    'Posterior Pharyngeal Bulb'), 1);
                neurons = ganglia(i).neurons;
            end
            names = neurons;
        end
        
        
        %% Head ganglia data.
        
        
        %% Midbody ganglia data.
        
        
        %% Tail ganglia data.
    end
end

