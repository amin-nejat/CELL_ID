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
        unselected_neuron_size = 40; % unselected neuron marker size
        selected_neuron_size = 200; % selected neuron marker size
        unselected_neuron_line = 1; % unselected neuron line size
        selected_neuron_line = 4; % selected neuron line size
    end
    
    % GUI public methods.
    methods (Static, Access = public)
        function obj = instance()
            %INSTANCE get the GUIPreferences singelton.
             persistent instance;
             if isempty(instance)
                 
                 % Add the image directory.
                 obj = Program.GUIPreferences();
                 user_dir = what('~/');
                 obj.image_dir = user_dir.path;
                 
                 % Save the instantiation.
                 instance = obj;
             else
               obj = instance;
             end
        end
    end
    
    % GUI private methods.
    methods (Static, Access = private)
        % Hide the constructor.
        function obj = GUIPreferences()
        end
    end
end
