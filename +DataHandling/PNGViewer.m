classdef PNGViewer
    %PNGVIEWER View PNG files.
    
    methods (Static)
        function show(file, name)
            %SHOW Show the PNG file.
            
            % Show the PNG.
            hFig = figure('Name', name, ...
                'Menubar', 'none', 'Toolbar', 'none');
            hIm = imshow(file);
            
            % Create a scroll panel for the PNG.
            imscrollpanel(hFig,hIm);
        end
    end
end
