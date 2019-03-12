function    str = value2short( val )
%   value2short converts value to a short string that is suitable to display
%
%   See also: mat2str

%   all( cellfun( 'isempty', this.SpringContainer ) )
%   ans =
%        1     0     1
%   value2short( all( cellfun( 'isempty', this.SpringContainer ) ) )
%   ans =
%   <1x3 logical>
%
%   value2short( {'a'} )
%   ans =
%   <1x1 cell>

%   2011-08-21, poi:  Wish, value2short( [ 1:2 ] ) should return '[1 2]'
%               workspacefunc( 'getshortvalue', [1:11] )
%               ans =
%               <1x11 double>
%               workspacefunc( 'getshortvalue', [1:10] )
%               ans =
%               [1 2 3 4 5 6 7 8 9 10]
%   2015-10-05, poi: value2short([1:2]) returns '[1 2]'. The limit len<=10 is hardcoded

%/  -------------------------------------------------
%   author:     per isakson
%   e-mail:     poi(at)bim-kth-se
%   created:    2008-06-10
%   modified:   2008-06-10
%/  -------------------------------------------------
    
%   2014-02-14, poi: rewrite
    if nargin > 0
        str     = workspacefunc( 'getshortvalue', val );
        max_len = 48;
        if length( str ) >= max_len
            str = [ str(1:max_len-4 ), ' ...' ];
        end
    else
        str = 'NIL';
    end
%{
%     if nargin > 0
%         
%         if numel( val ) >= 2  && not( ischar( val ) )           %   2011-04-11, poi: 
%             varargout = { workspacefunc( 'getshortvalue', val(1) ) };
%         else
%             varargout = { workspacefunc( 'getshortvalue', val    ) };
%         end
%     else
%         varargout = {'NIL'};
%     end
%}
end