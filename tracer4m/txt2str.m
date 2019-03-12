function    str = txt2str( filespec, nbytes )
% txt2str returns the content of text file as a string row vector.
%
%   str = txt2str( filespec )    returns the content of the text file filespec
%   str = txt2str( filespec, n ) returns the first n bytes of the content of filespec
%
%   The code is based on Matlab fileread. The file is opened in text mode, 'rt', thus 
%   '\r\n' is automatically replaced by '\n'. Force: simplify regular expressions. 
%
%   See also fread, fileread

    % Validate input args
    narginchk(1,2);
    
    if nargin == 1
       nbytes = inf; 
    end

    % get  filespec 
    if ~ischar( filespec ), 
        error(  'pia:txt2str:filenameNotString'                 ...
            ,   'The specifier, %s, is not a string', filespec  ); 
    end

    % do some validation
    if isempty( filespec ), 
        error( 'pia:txt2str:emptyFilename', 'Empty filespec' ); 
    end

    % open the file
    [fid, msg] = fopen( filespec, 'rt' );
    if fid == (-1)
        error( 'pia:txt2str:cannotOpenFile', 'Cannot open %s. %s', filespec , msg );
    end

    try
        % read file
        str = permute( fread( fid, nbytes, '*char' ), [2,1] );
    catch me
        % close file
        fclose(fid);
        throw(me);
    end

    % close file
    fclose(fid);
end