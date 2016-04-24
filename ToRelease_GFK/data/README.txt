This package contains the data and (GFK) code used in [1] and [2]. 

10 common categories are picked out from 4 datasets of object images: 
caltech-256 [3], amazon, webcam, and dslr [4]. We then follow the 
previously reported protocols for preparing features [4], i.e., extracting 
SURF features and quantizing them into an 800-bin histogram with codebooks 
computed via K-means on a subset of images from amazon. 

Data:
*_SURF_L10.mat:    features and labels
*_SURF_L10_imgs:   image indices of the data instances

Code:
GFK.m:             code for the geodesic flow kernel

Generally, data preprocessing helps:
load('amazon_SURF_L10.mat');
fts = fts ./ repmat(sum(fts,2),1,size(fts,2)); 
fts = zscore(fts,1);  

------------
[1] B. Gong, Y. Shi, F. Sha, and K. Grauman. Geodesic Flow Kernel for 
	Unsupervised Domain Adaptation. In CVPR 2012.
[2] B. Gong, K. Grauman, and F. Sha. Connecting the Dots with Landmarks: 
	Discriminatively Learning Domain-Invariant Features for Unsupervised 
	Domain Adaptation. In ICML 2013.
[3] Griffin, G., Holub, A., and Perona, P. Caltech-256 object category 
	dataset. Technical report, California Institute of Technology, 2007.
[4] Saenko, K., Kulis, B., Fritz, M., and Darrell, T. Adapting visual 
	category models to new domains. In ECCV, 2010.
		
-----------
Contact: Boqing Gong (boqinggo@usc.edu)