# Source Code for Multiple Sclerosis Lesion Change Detection using Multi-Scale Radiomic Features of MRI Images

The code is copied from the original repository to preserve anonymity, so history/changes are not visible.

### Prerequisites

The predictions and initial data are saved as NIfTI (Neuroimaging Informatics Technology Initiative) files which are loaded for each run. 21x21 patches are also in Hierarchical Data Format 5 (HDF5) files to be loaded rather than extracted on the fly in each experiment. They are not uploaded due to size limitations.

Run on Matlab R2013b. Other Matlab scripts necessary: Statistics and Machine Learning Toolbox, Tools for NIfTI and ANALYZE image

### Compared Methods in /methods

Method 1: detectionBaseline.m

Method 2: neighbors.m

Method 3: texturefeatures.m

Multi-scale Method: notensemble.m

### Other Files of Interest in /other

Creating 21x21 patches and texture-based feature extraction: createtexturepatches.m

Statistical validation: temp.m

Format ground truth: positivegtvalidate.m

### Future Work in Progress (To Be Continued) in /future

Automated prediction of future lesion change based on quantitative MR imaging features using cross-sectional image data: predictionBaseline.m

Pipeline could continue to be extended to evaluating MR images based on larger patches than 21x21 using supervised deep learning: original_cnn.ipynb, createpatches.ipynb
