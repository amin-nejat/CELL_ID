# CELL_ID
C-elegans Cell Identifier

- Open the get_params.m function first and change the directories to your local directory.
- For opening czi and tiff files both Fiji software and mij.jar are required. If any of these are missing the data must be in the .mat format and the dimensions must follow this format (x,y,z,c).
- Once updating get_params.m function is finished open visualize_light.mlapp using Matlab 2018a or 2018b.
- Type the dataset address and press "Load".
- In the "MP" tab select "Run MP" and select the parameters.
- Once running MP is finished go back to the first tab, check "Remove" and select the objects that are not cells. The software will remove them from the data structure. After removing all the non-cell objects uncheck "Remove".
- Go to "MP" tab again and select "Save MP".
- In the "MP" tab Select "Run AutoID" and once it finished by clicking on each cell in the "Preprocess" tab the ID that's assigned to it by the algorithm will be presented. The user can change the ID in the edit field next to it and by pressing "Go" the ID will be updated to user selection. Once annotation is finished Go to "MP" tab and select "Save AutoID". Note that for running AutoID two .mat files are required: "mog_model_tail.mat", "mog_model_head.mat".
