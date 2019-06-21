# CELL_ID
C-elegans Cell Identifier

**GUI environment** 


![GUI](https://dl.dropboxusercontent.com/s/d3lysgcr3dk8g5d/GUI-labeled.png)

(a) open and load a NeuroPAL image, (b) menu for basic preprocessing and improving the visualization, (c) automatically detect the neurons in NeuroPAL images, (d) automatically identify and label the detected neurons (e) clicking on each neurons shows the potential names that our method has found for that neuron and their associated uncertainties, (f) confirm the neuron identity, (g,h) manually annotate neurons, (i) selecting color channels used for visualizing NeuroPAL image, (j) select double click action: manually add a neuron/automatically find a neuron near the clicked location/remove artefacts, (k) list of neuron names sorted by their uncertainty according to the result of our method, (l) true neuron names and colors categorized by their corresponding ganglia, (m) z-stacks of NeuroPAL image and text labels for annotated or confirmed neurons, (n) controller for moving through different z-stacks, (o) maximum intensity projection of the image.

**Class, Package, and Usecase Diagrams**


To build a fast and scalable system we followed object oriented design practices. We identified internally coherent and externally independent classes in our environment and established their properties and behavior, as well as their interactions with other modules. We have separate packages for data handling, methods, logging, biological objects, and user interface. Each of these packages contains relevant classes and the functionalities of the system is built upon the interactions between the classes. For reusability purposes we tried to develop a modular system with independent modules. This allows the users of the system to reuse different compartments of the system for other purposes.

![Class](https://dl.dropboxusercontent.com/s/43l6x3vw55c8s5j/Class.png)
![Package](https://dl.dropboxusercontent.com/s/f26jeky6mzdlz3c/Package.png)
![Usecase](https://dl.dropboxusercontent.com/s/clgk3w8ju7d1i8n/Usecase.png)

**Related article** 

Here is a link to the draft of the paper: [CELL_ID](https://dl.dropboxusercontent.com/s/vxcq10ffdsgfa88/worm-neural-ID.pdf)

# Installation

To open NeuroPAL images the following package needs to be installed:

https://docs.openmicroscopy.org/bio-formats/5.3.4/users/matlab/index.html

We also use the following MATLAB functions but we have included them in the software:

https://github.com/tminka/lightspeed/blob/master/logsumexp.m
https://www.mathworks.com/matlabcentral/fileexchange/20328-munkres-assignment-algorithm
