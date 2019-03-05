function isOut = tracer( varargin )
% tracer is used in conditional breaks to transfer information to a log 
% 
%   Inputs
%       arg1    <calling object> or the type of the caller. (This is too complicated!)
%       arg2    'begin'/'end' indicating whether call is from begin or end of function
%                           
%   Output  
%       isOut   <false>, the call to tracer will not cause the conditional break to fire
%
% See also: TraceHistory, tracer_test 

%   author:     per isakson
%   e-mail:     per-ola-isakson(at)gmail-com
%   created:    2008-11-11
%   modified:   2010-10-03

%{
  2009-10-28
  Thread Subject: How to get debug to catch errors in listener callbacks?
  Subject: How to get debug to catch errors in listener callbacks?
  Date: 6 Dec, 2008 01:14:17

  Say I add a listener to an event called 'hello', defined for some object, 'obj':
      addlistener(obj, 'hello', @(src, evnt) MyLazyCallback(src, evnt, data));
  Now, if there is an error in MyLazyCallback, my Matlab only catches the error at 
  the line where obj is notified:
      ->notify(obj, 'hello')
  My question:  Is there any way to get Matlab to error inside MyLazyCallback? 
  It's a bit of a pain right now as I have to go through every callback listening 
  to that event to find the error...
  ....
  ....
  From: Ryan Ollos
  Date: 28 Oct, 2009 04:46:01
  I believe this has been fixed as of r2009b. poi: No, it is not fixed

2014-02-24, poi: What arguments are passed to the method is relevant in testing. 
    name    = evalin( 'caller', 'whos' );
    value   = evalin( 'caller', 'name' );
%} 
    assert( nargin == 2                                             ...
        ,   'tracer:WrongNumberInputArguments'                      ...
        ,   'Wrong, "%u", number of input arguments. Must be two.'  ...
        ,   nargin                                                  )
     
    stack   = dbstack(1);
    if numel( stack ) >= 2
        called_name = stack(1).name; 
        caller_name = stack(2).name; 
    else
        called_name = stack(1).name; 
        caller_name = 'base'; 
    end
    
    if isobject( varargin{1} )
        obj     = varargin{1};
        mobj    = metaclass( obj );
        mprops  = mobj.Properties;
        for ii = 1 : numel( mprops )
            if strcmp( 'ID_', mprops{ii}.Name )
                ID  = obj.ID_;
                break
            end
        end
    else
        ID = varargin{1};   
    end
%   2012-02-14, poi: There was never a variable named, create?    
%   if not( exist( 'created', 'var' ) )     
%       ID = '----';
%   end
%   2012-10-07, poi: Undefined function or variable "ID" for Hdf5Adapter.Hdf5Adapter 
    if not( exist( 'ID', 'var' ) )     
        ID = '---';
    else
%       2012-12-01, poi: this is special for the state pattern of test4ida       
        called  = regexp( called_name, '\.', 'split' );
        caller  = regexp( caller_name, '\.', 'split' );
        if strcmp( varargin{2}, 'end' )             ...
            && strcmp( caller{1}, 'IdaTester' )     ...
            && strcmp( called{1}(end-4:end), 'State' )
            
            obj = evalin('caller', 'this' );
            ID  = obj.context.ID;
        end
    end
    log = TraceHistory.Instance;
    log.add( caller_name, called_name, ID, varargin{2} )
    
    isOut = false;
end