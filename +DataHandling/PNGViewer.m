classdef PNGViewer
    %PNGVIEWER View PNG files.
    
    methods (Static)
        function show(file, name)
            %SHOW Show the PNG file.
            
            % Show the PNG.
            hFig = figure('Name', name, ...
                'Toolbar', 'none', 'Menubar', 'none');
            hIm = imshow(file);
            
            % Create a scroll panel for the PNG.
            imscrollpanel(hFig,hIm);
        end
    end
end
