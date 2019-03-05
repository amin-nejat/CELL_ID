function    sad = folder_list( varargin )
%   folder_list creates list of folders. Transposed dir-format. No '.' or '..'.
%   Purpose: remove clutter in code
%
%   See also: dir, folder_file_list, findfiles
 
%   2013-02-01, TODO: "_list" should be used for <1xn cell> exclusively(?)  
%
% FIXME: folder_list(varargin) with varargin={'c:\MyData\Vk\Pennfaktaren\txt','2012-05-07'}
%   returns empty, and with varargin = {'c:\MyData\Vk\Pennfaktaren\txt','2012-05-0*'}
%   returns 'c:\MyData\Vk\Pennfaktaren\txt\2012-05-07'. Example from Pennfaktaren to X

    sad = dir( fullfile( varargin{:} ) );
    sad( not( [sad.isdir] ) ) = [];
    for ii = numel( sad ) : -1 : 1
        if sad(ii).name(1) == '.'
            sad(ii) = [];
        end
    end
    sad( not( [sad.isdir] ) ) = [];
    sad = transpose( sad );
end