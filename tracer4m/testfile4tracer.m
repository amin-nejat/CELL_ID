function    testfile4tracer( varargin )
%testfile4tracer is used by unit tests in tracer_test

%   tracer_test/DispNestedFoo_1_test shall display
%{
--- tracer ---
tracer_test.DispNestedFoo_1_test
    testfile4tracer
        testfile4tracer/nested_01
            testfile4tracer/nested_01/nested_11
            testfile4tracer/nested_02
        testfile4tracer/nested_02
        subfunction_01
            subfunction_02
        subfunction_02
%}

%#ok<*SETNU>
%#ok<*DEFNU>
%#ok<*NASGU>

    abc = 1;
    def = 2; 
    ghi = 3;
    nested_01
    nested_02
    subfunction_01
    subfunction_02
    
    function    nested_01
        abc = 1;
        nested_11
        def = 2; 
        ghi = 3;
        function    nested_11
            abc = 1;
            def = 2; 
            ghi = 3;
        end
        nested_02
    end
    function    nested_02
        abc = 1; 
        def = 2; 
        ghi = 3;
    end
end
function    subfunction_01
    abc = 1;  
    subfunction_02
    def = 2; 
    ghi = 3;
end
function    subfunction_02
    abc = 1;
    def = 2; 
    ghi = 3;
end
