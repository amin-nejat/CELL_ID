# CELL_ID
C-elegans Cell Identifier

- Open the get_params.m function first and change the directories to your local directory.
- For opening czi and tiff files both Fiji software and mij.jar are required. If any of these are missing the data must be in the .mat format and the dimensions must follow this format (x,y,z,c).
- Once updating get_params.m function is finished open visualize_light.mlapp using Matlab 2018a or 2018b.
- Pess "Load" and select your dataset.
- Select "Detect Cells" and select the parameters.
- Once cell detection is finished use left click to remove the objects that are not cells. The software will remove them from the data structure.
- For adding cells simply use the middle click in the image.
- After the changes are done click on "Save" button.
- Select "ID Cells" and once it finished by clicking on each cell in the "Preprocess" tab the ID that's assigned to it by the algorithm will be presented. The user can change the ID in the edit field next to it and by pressing "Annotate" the ID will be updated to user selection.