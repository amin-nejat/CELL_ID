function params = get_params
% Constant parameters of the software
%
% Amin Nejat
params = [];
params.mp.default_window = '[15,15,9]';
params.mp.default_trunc = '90';
params.mp.default_k = '200';

% !!! ADDED FOR MY OWN CONVENIENCE GET RID OF THIS LATER :)
if strcmpi('ev', getenv('USER'))
    params.default_data_directory = '/Users/ev/Documents';
    params.fiji_lib_folder = '/Applications/Fiji.app';
    params.mij_path = '/Users/ev/Documents/Github/CELL_ID/data handling/mij.jar';
else % Amin's setup.
    params.default_data_directory = '/home/mn2822/Desktop/WormAutoID/data';
    params.fiji_lib_folder = '/home/mn2822/Desktop/WormTracking/libraries/Fiji.app';
    params.mij_path = '/home/mn2822/Desktop/WormTracking/libraries/Fiji.app/mij.jar';
end

end

