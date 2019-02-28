function names = stripNeuronLR(names, varargin)
%STRIPNEURONLR Strip a neuron name of the left/right character.
%
%   STRIPNEURONLR(NAMES)
%   STRIPNEURONLR(NAMES, KEEP_ASYMMETRIC_NAMES)
%
%   Input:
%   name - the names of the neuron to strip of the left/right character.
%   keep_asymmetric - should we keep keep asymmetric left/right neuron names?
%       default: false
%
%   Outputs:
%   name - the neuron name without the left/right character.

% Are we permitting left/right designations in asymmetric neurons?
keep_asym = false;
if ~isempty(varargin)
    keep_asym = varargin{1};
end

% Setup the don't strip list.
dont_strip_list = [];
if keep_asym
    dont_strip_list = {
        'ADL'
        'AQR'
        'ASEL'
        'ASER'
        'AVL'
        'AWC(OFF)'
        'AWC(ON)'
        'OLL'
        'PQR'
        'PVR'
        'RIR'
        'SDQL'
        'SDQR'
        };
    
else
        dont_strip_list = {
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
end

% Are we stripping a single neuron name?
is_char = ischar(names);
if is_char
    names = {names};
end

% Strip the neurons of their left/right character?
for i = 1:length(names)
    
    % Clean up AWC ON & OFF.
    name = names{i};
    if strncmp('AWC', name, 3)
        if keep_asym && length(name) > 3
            lr_char = lower(name(4));
            if lr_char == 'l' || lr_char == 'r'
                names{i} = [name(1:3) name(5:end)];
            end
        else
            names{i} = 'AWC';
        end
    end
    
    % Should we strip the neuron of its left/right character?
    if ~any(strcmp(dont_strip_list, name))
        lr_char = lower(name(end));
        if lr_char == 'l' || lr_char == 'r'
            names{i} = name(1:end-1);
        end
    end
end

% Are we stripping a single neuron name?
if is_char
    names = names{1};
end
end
