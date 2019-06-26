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

For more details about the experimental setup and the methods look at the following paper:

```
@article {Yemini676312,
	author = {Yemini, Eviatar and Lin, Albert and Nejatbakhsh, Amin and Varol, Erdem and Sun, Ruoxi and Mena, Gonzalo E. and Samuel, Aravinthan DT and Paninski, Liam and Venkatachalam, Vivek and Hobert, Oliver},
	title = {NeuroPAL: A Neuronal Polychromatic Atlas of Landmarks for Whole-Brain Imaging in C. elegans},
	elocation-id = {676312},
	year = {2019},
	doi = {10.1101/676312},
	publisher = {Cold Spring Harbor Laboratory},
	abstract = {Resolving whole-brain images of neuronal gene expression or neuronal activity patterns to the level of single neuron types with defined identities represents a major challenge. We describe here the development and use of a multicolor Caenorhabditis elegans transgene, called NeuroPAL (a Neuronal Polychromatic Atlas of Landmarks), to resolve unique neural identities in whole-brain images. NeuroPAL worms share a stereotypical multicolor map, permitting complete, unambiguous and automated determination of individual neuron identities in conjunction with GCaMP-based neuronal activity reporters and GFP/YFP/CFP gene-expression reporters. To demonstrate the method and its potential, we use NeuroPAL and GFP-based reporters to map expression for the whole family of metabotropic acetylcholine, glutamate, and GABA neurotransmitter receptors encoded in the C. elegans genome, revealing a vast number of potential molecular connections that go far beyond the anatomically-defined connectome. We then expand the technique to whole-brain activity, employing NeuroPAL and a panneuronal neural-activity sensor (GCaMP6s) for functional analysis. Using this tool we delineate extensive nervous system activity patterns in response to several stimuli with single, identified neuron resolution. We find that attractive odors sensed by the same neuron class exhibit dissimilar activity patterns implying that, despite their shared valence and stimulus modality, these odors drive distinct neural circuitry. Our results also indicate that the connectome is a poor predictor of functional activity, emphasizing the need for comprehensive brain-activity recordings that delineate behavior-specific circuitry. Lastly, we illustrate the NeuroPAL as an unbiased analysis tool for investigating neuronal cell fate in specific mutant backgrounds. With these applications in mind, we establish a high-throughput software pipeline for automated and semi-automated cell identification using NeuroPAL. In conclusion, we demonstrate the power of the NeuroPAL as a tool for decoding whole-brain gene expression and maps of functional activity.},
	URL = {https://www.biorxiv.org/content/early/2019/06/20/676312},
	eprint = {https://www.biorxiv.org/content/early/2019/06/20/676312.full.pdf},
	journal = {bioRxiv}
}
```

Here is a link to biorxiv the paper: [CELL_ID](https://www.biorxiv.org/content/10.1101/676312v1)

# Installation

To open NeuroPAL images the following packages needs to be installed:

https://docs.openmicroscopy.org/bio-formats/5.3.4/users/matlab/index.html
https://github.com/StackOverflowMATLABchat/mlapptools

We also use the following MATLAB functions but we have included them in the software:

https://github.com/tminka/lightspeed/blob/master/logsumexp.m
https://www.mathworks.com/matlabcentral/fileexchange/20328-munkres-assignment-algorithm
