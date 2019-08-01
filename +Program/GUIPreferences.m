classdef GUIPreferences < handle
    %GUIPreferences GUI preferences.
    
    % GUI properties.
    properties (Access = public)
        
        % Version.
        %version = Program.ProgramInfo.version;
        
        % Preferences file.
        prefs_file = '.NeuroPAL_ID_prefs.mat';
        
        % Image directory info.
        image_dir = [];
        
        % Display info.
        GFP_color = [1,1,0]; % GFP reporter color for GUI
        neuron_dot; % neuron dot info: (un)selected marker + line sizes
    end
    
    % GUI public static methods.
    methods (Static, Access = public)
        function obj = instance()
            %INSTANCE get the GUIPreferences singelton.
             persistent instance;
             if isempty(instance)
                 
                 % Add the image directory.
                 obj = Program.GUIPreferences();
                 user_dir = what('~/');
                 obj.image_dir = user_dir.path;
                 
                 % Initilaize the neurons dots.
                 obj.neuron_dot.unselected.marker = 40; % unselected neuron marker size
                 obj.neuron_dot.unselected.line = 1; % unselected neuron line size
                 obj.neuron_dot.selected.marker = 200; % selected neuron marker size
                 obj.neuron_dot.selected.line = 4; % selected neuron line size
                 
                 % Save the instantiation.
                 instance = obj;
             else
               obj = instance;
             end
        end
        
        function save()
            %SAVE save the GUI preferences to their file.
            
            % Instantiate the class.
            obj = Program.GUIPreferences.instance();
            
            % Save the preferences.
            version = Program.ProgramInfo.version;
            GUI_prefs = obj;
            save(obj.prefs_file, 'version', 'GUI_prefs');
        end
        
        function is_loaded = load()
            %LOAD load the GUI preferences from their file.
            
            % Instantiate the class.
            obj = Program.GUIPreferences.instance();
            
            % Load the preferences.
            is_loaded = false;
            if exist(obj.prefs_file, 'file')
                is_loaded = true;
                prefs = load(obj.prefs_file);
                Program.GUIPreferences.read(prefs.GUI_prefs);
            end
        end
        
        function read(prefs)
            %READ read and store the GUI preferences.
            
            % Instantiate the class.
            obj = Program.GUIPreferences.instance();
            
            % Read the preferences.
            obj.prefs_file = prefs.prefs_file;
            obj.image_dir = prefs.image_dir;
            obj.GFP_color = prefs.GFP_color;
            obj.neuron_dot = prefs.neuron_dot;
        end
        
        function dot = inputNeuronDot()
            %INPUTNEURONDOT get the neuron dot info from the user.
            
            % Instantiate the class.
            obj = Program.GUIPreferences.instance();
            
            % Initialize the input dialog.
            title = 'Neuron Dots';
            prompt = {'Selected neuron dot size:', ...
                'Unselected neuron dot size:', ...
                'Selected neuron line size:', ...
                'Unselected neuron line size:'};
            definput = {num2str(obj.neuron_dot.selected.marker), ...
                num2str(obj.neuron_dot.unselected.marker), ...
                num2str(obj.neuron_dot.selected.line), ...
                num2str(obj.neuron_dot.unselected.line)};
            dims = [1 35];
            
            % Get the user input.
            dot = [];
            answer = inputdlg(prompt, title, dims, definput);
            if isempty(answer)
                return;
            end
            
            % Parse the user input.
            dot.selected.marker = round(str2double(answer{1}));
            dot.unselected.marker = round(str2double(answer{2}));
            dot.selected.line = round(str2double(answer{3}));
            dot.unselected.line = round(str2double(answer{4}));
            
            % Check the user input.
            if dot.selected.marker < 1
                dot.selected.marker = obj.neuron_dot.selected.marker;
            end
            if dot.unselected.marker < 1
                dot.unselected.marker = obj.neuron_dot.unselected.marker;
            end
            if dot.selected.line < 1
                dot.selected.line = obj.neuron_dot.selected.line;
            end
            if dot.unselected.line < 1
                dot.unselected.line = obj.neuron_dot.unselected.line;
            end
            
            % Store the user input.
            obj.neuron_dot.selected.marker = dot.selected.marker;
            obj.neuron_dot.unselected.marker = dot.unselected.marker;
            obj.neuron_dot.selected.line = dot.selected.line;
            obj.neuron_dot.unselected.line = dot.unselected.line;
        end
    end
    
    % GUI private methods.
    methods (Static, Access = private)
        % Hide the constructor.
        function obj = GUIPreferences()
        end
    end
end
