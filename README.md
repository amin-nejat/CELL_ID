# CELL_ID
C-elegans Cell Identifier

- Open the get_params.m function first and change the directories to your local directory.
- For opening czi and tiff files Fiji software is required which is publibcly available here: https://imagej.net/Fiji/Downloads. If any of these are missing the data must be in the .mat format and the dimensions must follow this format (x,y,z,c).
- Once updating get_params.m function is finished open visualize_light.mlapp using Matlab 2018b.
- Pess "Load Data" and select your dataset. The pre-processing menu contains some basic pre-precessing steps that you can apply before analyzing the image. These pre-processing steps include spatial decimation, thresholding, z-scoring, and denoising (it is recommended to run the cell detection procedure on z-scored image).
- Select "Detect Neurons" and select the parameters. If you've ran the cell detection before you can reload the pre-stored data structure that contains the information about the cells (location, covariance, signal and noise color).
- You can reconstruct the neuron shapes and color by pressing "Reconstruct" button in the Reconstruction tab.
- Once cell detection is finished use double click to remove the objects that are not cells or add the cells that are not detected by the automatic cell detection procedure.
- After the changes are done click on "Save" button.
- Select "ID Neurons" and choose the parameters for automatic cell identification. Once it's finished by left-clicking on each cell the ID that is assigned to it by the algorithm will be displayed. The user can change the ID in the edit field next to it and by pressing "Annotate" the ID will be updated to user selection.
