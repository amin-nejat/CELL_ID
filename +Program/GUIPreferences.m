classdef GUIPreferences < handle
    %GUIPreferences GUI preferences.
    
    % GUI properties.
    properties (Access = public)
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
                obj = Program.GUIPreferences();
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
