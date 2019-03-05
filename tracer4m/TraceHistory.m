classdef ( Sealed )  TraceHistory < handle    
% Store and displays a trace of calls (/messages) to functions and methods
%
%   Public properties   none
%
%   Public methods          
%       add             -   used by the function, tracer, to add info on calls   
%       clearHistory    -   remove trace of calls
%       disp            -   overloaded 
%       get             -       
%       Instance        -   returns handle to the TraceHistory object (Singleton)
%       setup           -   registers file(s) for tracing, i.e. sets breakpoints  
%   
%   Singleton pattern  
%   **NOTE** the methods are not static, i.e log = TraceHistory.Instance, etc.
%
%   See also: tracer, tracer_test, tracer4m_demo, Messenger, TestLogger

%{
%   Example: 
    log = TraceHistory.Instance;
    log.setup( { 'testfile4tracer' } )  
    testfile4tracer;
    disp(  log )
    report = log.get;
%}
%   Bugs
%   2012-02-14, poi: Classes being traced cannot have a property named "ID".
%   2011-04-30, poi: Functions must have a closing "end"

%   HISTORY
%   2016-08-17, poi: QuickFix: R2016a, replaced
%               sam = transpose( mlint( '-calls', char( filespec ) ) ); by
%               str = evalc('mlint( ''-calls'', char( filespec ) )');   etc.

%   2016-09-07, poi: getcallinfo, Returns called functions and their first and last lines.
%               This function is unsupported and might change or be removed without notice
%               in a future version.

%   profile( ..., '-history' )
%{
%   profile( ..., '-history' ) when did it appear. Search comp.soft-sys.matlab 
%   Subject: Function call tree, from: Etta Makwetta, Date: 13 Sep, 2002 20:26:48
%   Etta, While it won't give you a tree, the DEPFUN .../Bob Gilmore, The MathWorks
%   ... take your code as an input and then outputs the structure i.e who's calling who
%   ... Try this: ... profile on -history, surf(peaks), stats = profile('info');
%   Steve Lord, 2003-06-10

%   2016-09-08, http://undocumentedmatlab.com/blog/undocumented-profiler-options-part-3
%   The history data is actually a numeric matrix, where the first row contains the values
%   0 (=function entry) or 1 (=function exit), and the second row is the corresponding
%   index into profData.FunctionTable, indicating the called function. We can easily
%   convert this matrix into human-readable form using the following code snippet:

%     offset = cumsum(1-2*history(1,:)) - 1;  % calling depth
%     entryIdx = history(1,:)==1;     % history items of function entries
%     funcIdx = history(2,entryIdx);  % indexes of relevant functions
%     funcNames = {profData.FunctionTable(funcIdx).FunctionName};
%     for idx = 1: length(funcNames);
%        disp([repmat(' ',1,offset(idx)) funcNames{idx}]);
%     end
% 
%   which generates the following calling list in the MATLAB Command Window:
% 
%     isempty
%      isempty
%       transpose
%      meshgrid
%       peaks
%      nargchk
%       error
%      ishg2parent
%     ...
%}

%   author:     per isakson
%   e-mail:     per-ola-isakson(at)gmail-com
%   created:    2008-11-11
%   modified:   2016-09-11

%#ok<*PRTCAL> allow output to command screen
%#ok<*AGROW>  allow grow inside loop

    properties  ( Access = private )        %
        History
    end   
    methods     ( Access = private )        %
        function  this = TraceHistory()                      
            this.History  = cell( 5, 0 );    
        end
    end
    methods     ( Access = public  )        %
        
        function    add( this, caller_name, called_name, ID, begin_end )  
%           add is called by tracer            
            time_stamp      = datestr( now, 'HH:MM:SS,FFF' );
            this.History    = cat( 2, this.History          ...
                            ,   {   caller_name             ... from the stack      1
                                ;   called_name             ... from the stack      2
                                ;   ID                      ... #nnn or C/M/F/S/N   3
                                ;   time_stamp              ... of call to tracer   4
                                ;   begin_end   }   );      ... 'begin' or 'end'    5
%           pause( 0.001 )
        end
        %
        function    clearHistory            ( this           )  
            % the name "clear" is taken - overload - no
            this.History  = cell( 5, 0 );
        end
        function    disp                    ( this           )  
% FIXME: Overloading to show the trace of calls might be nice when debugging, but ...
            this.flip_ID_of_constructors();
            intent  = '';
            fprintf( 1, '%s\n', '--- tracer4m ---' )
            for ca = this.History

                if isempty( intent )
                    fprintf( 1, '%s\n', ca{1} )                         % caller
                end
                if strcmp( ca{5}, 'begin' )
                    intent( end+1 : end+4 )  = ' ';
                    fprintf( 1, '%s%s\n', intent, ca{2} )               % called

                elseif  strcmp( ca{5}, 'end' )
                    intent  = intent( 1 : end-4 );
                else
                    error( 'TraceHistory:disp: ') 
                end
                if  strcmp( ca{3}, 'A' )            
                    intent  = intent( 1 : end-4 );
                end
            end
        end
        function    log = get               ( this           )  
            this.flip_ID_of_constructors();
            log = this.History;
        end
        function    setup                   ( this, caFileList )  
            
%   ----            
%   2012-11-27, poi: Added wildcard functionality

            is_wildcard = any( not( cat( 2                                      ...
                ,   cellfun( @(str) isempty(strfind(str,'*')), caFileList )     ...
                ,   cellfun( @(str) isempty(strfind(str,'?')), caFileList ) ) ) );
                
            if is_wildcard
                ffs = cell( 2, 0 );
                for cac = caFileList
                    [ folder_name, file_name, ext ] = fileparts( cac{:} );
                    ffs = cat( 2, ffs, { folder_name ; [ file_name, ext ] } );
                end
                
                ffs = TrenddataShuffler.expand_folder_file_spec_list            ...
                                                (   'FolderFileSpecList', ffs   );
                            
                is_invalid  = cellfun( @(str) strncmp(str,'_',1), ffs(2,:) );
                ffs( :, is_invalid ) = [];
                
                file_name_list  = cellfun(  @(p,f) fullfile(p,f)            ...
                                        ,   ffs(1,:), ffs(2,:), 'uni', false );
                            
            else
                file_name_list  = caFileList;
            end
            
%           Is it possible to scan for Code Analyzer errors?
%             for cac = file_name_list
%                if not( isempty( checkcode(cac{:}) ) )
%                    17;
%                end
%
%           In summary, you could simply do:
%           errMsgs = mlint('-m2',srcFileNames); % m2 = errors only
%           m1Msgs  = mlint('-m1',srcFileNames); % m1 = errors and severe warnings only
%           allMsgs = mlint('-m0',srcFileNames); % m0 = all errors and warnings
%           Note that mlint returns the data in struct format, while mlintmex returns 
%           a string that you must parse.
%           Yair Altman
%   ----            

            for caFL = file_name_list
                %   this.curStatus = dbstatus( caFL{:} ) and
                %   dbstop( this.curStatus ) at the "end" is more professional.  
                %   However, for now I want to keep the trace-breakpoints to have a 
                %   chance to inspect them. 
                saStatus = transpose( dbstatus( caFL{:} ) );
                for sa = saStatus
                    for ii = 1 : numel( sa.expression )
                        if not( isempty( sa.expression{ii} ) )
                            dbclear( 'in', caFL{:}, 'at', num2str( sa.line(ii) ) )
                        end
                    end
                end
                calls = MethodInfo( caFL{:} ); 

                for s = calls
                    switch s.FunctionType
                        case { 'Constructor' }
                            
                            dbstop( 'in',   s.Name{1}                       ...
                                ,   'at',   s.Range{1,1}                    ...
                                ,   'if',   'tracer( ''C'', ''begin'' )'    )

                            dbstop( 'in',   s.Name{1}                       ...
                                ,   'at',   s.Range{1,2}                    ...
                                ,   'if',   break_condition_for_class_( s ) )
                            
                        case { 'Method', 'Property' }
                            
                            dbstop( 'in',   s.Name{1}                       ...
                                ,   'at',   s.Range{1,1}                    ...
                                ,   'if',   break_condition_for_class_( s ) )
                            
                            dbstop( 'in',   s.Name{1}                       ...
                                ,   'at',   s.Range{1,2}                    ...
                                ,   'if',   'tracer( ''M'', ''end''   )'    )
                            
                        case { 'Main', 'Sub', 'Nested' }

                            dbstop( 'in',   s.Name{1}                       ...
                                ,   'at',   s.Range{1,1}                    ...
                                ,   'if',   'tracer( ''F'', ''begin'' )'    )

%   FIXME:  You cannot set a breakpoint past the start of the last expression in the
%            file. Error in TraceHistory/setup (line 306), dbstop( 'in',   s.Name{1}

                            dbstop( 'in',   s.Name{1}                       ...
                                ,   'at',   s.Range{1,2}                    ...
                                ,   'if',   'tracer( ''F'',  ''end''   )'   )
                            
%FIXME: dbstop fails with methods. See support answer; [THREAD ID: 1-DXEKZW]

                        case { 'Static' }
                            
                            dbstop( 'in',   s.Name{1}                       ...
                                ,   'at',   s.Range{1,1}                    ...
                                ,   'if',   'tracer( ''S'', ''begin''  )'   )
                            
                            dbstop( 'in',   s.Name{1}                       ...
                                ,   'at',   s.Range{1,2}                    ...
                                ,   'if',   'tracer( ''S'', ''end''   )'    )
                            
                        case { 'Anonymous' }
                            
% FIXME: Anonymous function on continuation lines causes trouble. mlint returns  
%        the continuation line number, which when used by dbstop puts the break  
%        on first executable following the continuation line. Skip anonymous 
%        function for now!                               
%{
                            dbstop( 'in',   s.Name{1}                       ...
                                ,   'at',   [s.Range{1,1},'@']              ...
                                ,   'if',   'tracer( ''A'', ''begin'' )'    )
%}                        
                        otherwise
                            warning('TraceHistory:setup:FunctionType'   ...
                                ,   'Unknown function type: "%s"'       ...
                                ,   s.FunctionType                      )
                    end
                end
            end
            this.clearHistory()
        end
    end
    methods     ( Access = private )        %
        function    flip_ID_of_constructors ( this          )   %
            this.History = flip_ID_of_constructors_( this.History );
        end
    end
    methods     ( Static )                  %
        function    this = Instance()
            persistent  Instance
            if isempty( Instance ) || not( isvalid( Instance ) )
                Instance = TraceHistory;
            end
            this = Instance;
        end
    end
end
function    calls = MethodInfo              ( file_name )   %
    
    filespec = which( file_name, '-all' );   
    filespec = filespec(1,1);   % Q&D
    
    assert( not( isempty( filespec ) )              ...
        ,   'TraceHistory:MethodInfo:FileNotFound'  ...
        ,   'File not found: "%s"'                  ...
        ,   filespec                                )

    assert( numel( filespec ) == 1                      ...
        ,   'TraceHistory:MethodInfo:ShadowFile'        ...
        ,   'There are more than one file named: "%s"'  ...
        ,   file_name                                   )
    
%   2016-08-17, poi: Code broken by R2016a(?). The mlint-call returns a scalar structure
%       instead of an expected array. The scalar structure is the first item of the 
%       expected structure array. 
%
%       sam = mlint( '-calls', char( filespec ) )
%           sam = 
%               message: 'S0 29 28 testclass…'
%               ...  
%   
%   I failed to use checkcode (which quitely returns empty for some files) and  
%   getcallinfo (which returns data on a totally different format)
%
%   QuickFix:
%
    str = evalc('mlint( ''-calls'', char( filespec ) )');
    cac = textscan( str, '%s', 'Delimiter', '\n' );
    sam = cell2struct( cac{1}, {'message'}, 2 );
    ise = cellfun( @isempty, {sam.message}, 'uni',true );
    sam(ise) = []; % empty rows
    
%   FunctionType /http://undocumentedmatlab.com/blog/function-definition-meta-info 
%   M   main (top-level) function
%   S   sub-function
%   N   nested function
%   U   out-of-scope (external/built-in) function
%   A   anonymous function
%   E   end-of-function indication

    calls = struct( 'Name'          ,   {}  ...        
                ,   'MlintType'     ,   {}  ...
                ,   'FunctionType'  ,   {}  ...      
                ,   'MlintLevel'    ,   {}  ...
                ,   'Range'         ,   {}  );
    
    for ii = 1 : numel( sam )
        s1  = regexp( sam(ii).message                                               ...
            , '^(?<type>\w)(?<level>\d) (?<line>\d+) (?<column>\d+) (?<name>.*+)$'  ...
            , 'names'                                                               );

        if  any( strcmp( s1.type, {'E','U'} ) ),    continue
        end
        
        if not( strcmp( s1.type, 'A' ) )    % Increment the number; 
            s1.line = sprintf( '%.0f', sscanf( s1.line, '%u' ) + 1 );
        else
            % do nothing
        end
        
        s2  = regexp( sam(ii+1).message                                             ...
            , '^(?<type>\w)(?<level>\d) (?<line>\d+) (?<column>\d+) (?<name>.*+)$'  ...
            , 'names'                                                               );
        
        if strcmp( s1.type, 'A' ),  calls(end+1).Name = { file_name; 'anonymous'      };
        else                        calls(end+1).Name = { file_name; strtrim(s1.name) };
        end                 
        calls(end).MlintType    = s1.type;
        calls(end).MlintLevel   = s1.level;
        calls(end).Range        = { s1.line, s2.line ; s1.column, s2.column };   
    end
    
    cac         = NameParts( file_name );
    metaobject  = meta.class.fromName( cac{end} );
    
    if not( isempty( metaobject ) )
        MethodNames   = cellfun( @(obj) obj.Name , metaobject.Methods   , 'uni', false );
        PropertyNames = cellfun( @(obj) obj.Name , metaobject.Properties, 'uni', false );
        SetGetNames   = [ cellfun( @(c) ['get.',c], PropertyNames, 'uni', false )
                          cellfun( @(c) ['set.',c], PropertyNames, 'uni', false ) ];
    else
        MethodNames   = {''};
        SetGetNames   = {''};
    end
    
    for ii = 1 : numel( calls )
        ism = strcmp( calls(ii).Name{2}, MethodNames );
        isp = strcmp( calls(ii).Name{2}, SetGetNames );
        if any( ism )
            if metaobject.Methods{ism}.Static
                calls(ii).FunctionType              = 'Static';
            else
                if  strcmp( calls(ii).Name{1} ...  
                        ,   calls(ii).Name{2} )      
                    calls(ii).FunctionType          = 'Constructor';
                else
                    calls(ii).FunctionType          = 'Method';
                end
            end
        elseif any( isp )
            calls(ii).FunctionType                  = 'Property';
        else
            switch calls(ii).MlintType
                case 'M',   calls(ii).FunctionType  = 'Main';   
                case 'S',   calls(ii).FunctionType  = 'Sub';   
                case 'N',   calls(ii).FunctionType  = 'Nested';   
                case 'A',   calls(ii).FunctionType  = 'Anonymous';   
                otherwise
                    warning('TraceHistory:MethodInfo:UnknownMlintType'  ...
                        ,   'Unknown mlint function type: "%s"'         ...
                        ,   calls(ii).MlintType                         )
            end
        end
    end
end
function    cac = NameParts                 ( fullname  )   %

%   Doc says: Packages are special folders that can contain class folders, function   
%   and class definition files, and other packages. ...  mypack.mysubpack.myfcn 
%   Assumption: 
%   1.  "set" and "get" are not used as package names! 
%   2.  No packages in packages

    cac = textscan( fullname, '%s', 'delimiter', '.', 'whitespace', '' );
    cac = cac{:};
    ism = ismember( cac, { 'set', 'get' } );
    
    if any( ism )
        ix  = find( ism );
        assert( numel( ix == 1 )                                ... 
            ,   'TraceHistory:NameParts:SetGetTrouble'          ...
            ,   'More than on occurance of "set/get" in "%s"'   ...
            ,   fullname                                        )
        if ix == numel( ism )
            % fine - do nothing
        else
            cac{ ix+1 } = cat( 2, cac{ix}, '.', cac{ix+1} );
            cac( ix   ) = [];
        end
    end
end
function    out = break_condition_for_class_( sas       )   %
    
    ffs = which( [sas.Name{1},'.m'] );
    str = txt2str( ffs );
    
    num = regexp( str, '(?<=(^|\n)[ ]*)classdef', 'once' );
    assert( not(isempty( num )), 'tracer4m:TraceHistory:NotClassdef' ...
        ,   'Cannot find the string, "classdef", in "%s"', ffs      )
    
    regstr2 = regexptranslate( 'escape', sas.Name{2} );
    if strcmp( sas.Name(1), sas.Name(2) )
        % constructor
        xpr = ['(?<=function\s*)\S+(?=\s*\=\s*',regstr2,')'];
        tag = 'end';    % Q%D: the object is not available at the "begin"
    else
        xpr = ['(?<=',regstr2,'\s*\(\s*)\S+?(?=[ ,\)])'];
        tag = 'begin';
    end
    obn = regexp( str, xpr, 'match', 'once' );
    out = sprintf( 'tracer( %s, ''%s'' )', obn, tag ); 
end
function    cac = flip_ID_of_constructors_  ( cac       )   %

    ix_constructor = find( ismember( cac(3,:), 'C'     ) ...
                         & ismember( cac(5,:), 'begin' ) );
    
    if isempty( ix_constructor )
        return                          %   RETURN
    end
    
    out = cac;

    for ixc = ix_constructor
        class_spec  = cac( 1, ixc );
        constructor = cac( 2, ixc );

        is_ends = strcmp( cac(1,:), class_spec )  ...
                & strcmp( cac(2,:), constructor ) ...
                & not( strcmp( cac(3,:), 'C' ) )  ...
                & strcmp( cac(5,:), 'end' )       ;

        ix_ends = find( is_ends );

        ixe = ix_ends( find( ix_ends > ixc, 1, 'first' ) );

        out(3,[ixc,ixe]) = cac(3,[ixe,ixc]);
    end
    cac = out;
end

