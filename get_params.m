function params = get_params
% Constant parameters of the software
%
% Amin Nejat
params = [];
params.mp.default_window = '[13,13,5]';
params.mp.default_trunc = '60';
params.mp.default_k = '180';
params.mp.max_eig_thresh = '10';

% !!! ADDED FOR MY OWN CONVENIENCE GET RID OF THIS LATER :)
if strcmpi('ev', getenv('USER'))
    params.default_data_directory = '/Users/ev/Documents/Github';
    params.fiji_lib_folder = '/Applications/Fiji.app';
elseif strcmpi('Erdem',getenv('USERNAME'))
    params.default_data_directory = 'C:\Users\Erdem\Documents\GitHub';
    params.fiji_lib_folder = 'C:\Users\Erdem\Dropbox\Projects\worm\Fiji.app';
else % Amin's setup.
    params.default_data_directory = '/home/mn2822/Desktop/WormAutoID/data';
    params.fiji_lib_folder = '/home/mn2822/Desktop/WormTracking/libraries/Fiji.app';
end

end

