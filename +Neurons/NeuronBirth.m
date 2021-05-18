%% Construct the neuron birth data.
classdef NeuronBirth
    methods (Static)
    
        %% Neuron text data.
        
        % Get the hermaphrodite neuron birth stages as text.
        function [neurons, stages, neuron_stages] = getHermaphroditeBirthText()
            %GETHERMAPHRODITEBIRTHTEXT Get a list of hermaphrodite neuron
            % birth times as a text string.
            %
            % Output:
            %   neurons = neuron names
            %   stages = birth time text
            %   neuron_stages = neuron name and birth time text
            import Neurons.*;
            
            % Generate the persisitent data.
            persistent neuron_data;
            persistent stage_data;
            persistent neuron_stage_data;
            if isempty(neuron_data)
                 [neuron_data, ~, stages, times] = ...
                     NeuronBirth.getHermaphroditeBirthTimes();
                 stage_data = cellfun(@(x,y) sprintf('%s %s', x, y), ...
                     times, stages, 'UniformOutput', false);
                 neuron_stage_data = cellfun(@(x,y,z) ...
                     sprintf('%s (%s %s)', x, y, z), ...
                     neuron_data, times, stages, 'UniformOutput', false);
            end
            
            % Return the data.
            neurons = neuron_data;
            stages = stage_data;
            neuron_stages = neuron_stage_data;
        end
        
        % Get the male neuron birth stages as text.
        function [neurons, stages, neuron_stages] = getMaleBirthText()
            %GETMALEBIRTHTEXT Get a list of male neuron birth times as a
            % text string.
            %
            % Output:
            %   neurons = neuron names
            %   stages = birth time text
            %   neuron_stages = neuron name and birth time text
            import Neurons.*;
            
            % Generate the persisitent data.
            persistent neuron_data;
            persistent stage_data;
            persistent neuron_stage_data;
            if isempty(neuron_data)
                 [neuron_data, ~, stages, times] = ...
                     NeuronBirth.getMaleBirthTimes();
                 stage_data = cellfun(@(x,y) sprintf('%s %s', x, y), ...
                     times, stages, 'UniformOutput', false);
                 neuron_stage_data = cellfun(@(x,y,z) ...
                     sprintf('%s (%s %s)', x, y, z), ...
                     neuron_data, times, stages, 'UniformOutput', false);
            end
            
            % Return the data.
            neurons = neuron_data;
            stages = stage_data;
            neuron_stages = neuron_stage_data;
        end
        
        
        %% Neuron birth data.
        
        % Get all hermaphrodite neuron birth data.
        function [names, minutes, stages, times] = getHermaphroditeBirthTimes()
            %GETHERMAPHRODITEBIRTHTIMES Get a list of hermaphrodite neuron
            % birth times.
            %
            % Output:
            %   names = neuron names
            %   minutes = neuron birth times in minutes
            %   stages = larval stage
            %   times = stage-specific birth time
            import Neurons.*;
            
            % Generate the persisitent data.
            persistent name_data;
            persistent minute_data;
            persistent stage_data;
            persistent time_data;
            if isempty(name_data)
                
                % Load the birth times.
                load('birth_times.mat', 'NeuronBirthTimes');
                name_data = NeuronBirthTimes.neuron;
                minute_data = NeuronBirthTimes.minutes;
                
                % Generate the stages and times.
                [stage_data, time_data] = ...
                    NeuronBirth.getHermaphroditeStages(name_data, minute_data);
                
                % Sort the names.
                [name_data, sort_i] = sort(name_data);
                minute_data = minute_data(sort_i);
                stage_data = stage_data(sort_i);
                time_data = time_data(sort_i);
            end
            
            % Return the data.
            names = name_data;
            minutes = minute_data;
            stages = stage_data;
            times = time_data;
        end
        
        % Get all male neuron birth data.
        function [names, minutes, stages, times] = getMaleBirthTimes()
            %GETMALEBIRTHTIMES Get a list of male neuron
            % birth times.
            %
            % Output:
            %   names = neuron names
            %   minutes = neuron birth times in minutes
            %   stages = larval stage
            %   times = stage-specific birth time
            import Neurons.*;
            
            % Generate the persisitent data.
            persistent name_data;
            persistent minute_data;
            persistent stage_data;
            persistent time_data;
            if isempty(name_data)
                
                % Get the hermaphrodite data.
                [name_data, minute_data, stage_data, time_data] = ...
                    NeuronBirth.getHermaphroditeBirthTimes();
                
                % Load the male neurons.
                load('male_data.mat', 'neurons');

                % Add the male neurons.
                male_neurons = setdiff(neurons, name_data);
                num_male = length(male_neurons);
                name_data = [name_data; male_neurons];
                minute_data(end+1:end+num_male) = nan;
                stage_data(end+1:end+num_male) = {'L4'};
                time_data(end+1:end+num_male) = {'Late'};
                
                % Sort the names.
                [name_data, sort_i] = sort(name_data);
                minute_data = minute_data(sort_i);
                stage_data = stage_data(sort_i);
                time_data = time_data(sort_i);
            end
            
            % Return the data.
            names = name_data;
            minutes = minute_data;
            stages = stage_data;
            times = time_data;
        end
        
        
        %% Neuron birth stages and stage-specific birth times.
        
        % Get the hermaphrodite neuron birth stages and stage-specific times.
        function [stages, times] = getHermaphroditeStages(neurons, minutes)
            %GETHERMAPHRODITESTAGES Get a list of hermaphrodite neuron
            % birth stages and stage-specific times.
            %
            %  Input:
            %   names = neuron names
            %   minutes = neuron birth times in minutes
            %
            %  Output:
            %   stages = larval stage
            %   times = stage-specific birth time
            import Neurons.*;
            
            % Initialize the data.
            stages = cell(length(minutes),1);
            times = cell(length(minutes),1);
            
            % Generate the embryonic times.
            stages(minutes < 1160) = {'Embryo'};
            times(minutes < 375) = {'Early'};
            times(minutes >= 375 & minutes <= 435) = {'Mid'};
            times(minutes > 435 & minutes < 1160) = {'Late'};
            
            % Generate the L1 times.
            stages(minutes >= 1160 & minutes < 1730) = {'L1'};
            times(minutes >= 1160 & minutes <= 1490) = {'Mid'};
            times(minutes > 1490 & minutes < 1730) = {'Late'};
            
            % Generate the L2 times.
            stages(minutes >= 1730 & minutes <= 2240) = {'L2'};
            times(minutes >= 1730 & minutes < 1880) = {'Early'};
            times(minutes >= 1880 & minutes <= 2000) = {'Mid'};
            times(minutes > 2000 & minutes <= 2240) = {'Late'};
            
            % Correct the late-born neurons.
            late_neurons = NeuronBirth.late_herm_neurons;
            for i = 1:length(late_neurons)
                j = startsWith(neurons, late_neurons{i});
                stages(j) = NeuronBirth.late_herm_neuron_stages(i);
                times(j) = NeuronBirth.late_herm_neuron_times(i);
            end
        end
    end
    
%% Constant properties.
    properties (Constant)
        
        % Hermaphrodite neurons that differentiate long after they are born.
        late_herm_neurons = {
            'PDA'
            'HSN'
            'VC'
            };
        late_herm_neuron_stages = {
            'L3'
            'Adult'
            'Adult'
            };
        late_herm_neuron_times = {
            'Early'
            'Early'
            'Early'
            };
    end
end
