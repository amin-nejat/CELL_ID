classdef Neuron < handle
    %NEURON A neuron.
    %
    %   A neuron has a position within an image, color, list of potential
    %   IDs, probabilities associated with these IDs, and methods to draw
    %   it within a GUI.
    
    properties (Access = private)
        brain
    end
    
    properties
        position % neuron pixel position (x,y,z)
        color % neuron color (R,G,B,W), W = white channel, values=[0-255]
        ids % neuron IDs listed by descending probability
        id_probs % neuron ID probabilities
    end
    
    methods
        function addToBrain(obj, b)
            %ADDTOBRAIN Add this neuron to a brain image.
            obj.brain = b;
        end

        function remove(obj)
            %REMOVE Remove this neuron.
            if ~isempty(obj.brain)
                obj.brain.remove(obj);
            end
        end
    end
end
