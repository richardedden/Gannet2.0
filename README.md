# Gannet 3.0

Open-source, MATLAB-based software for automated data processing and quantification of edited magnetic resonance spectroscopy (MRS) data.

Visit our [website](http://www.gabamrs.com/) for the latest news on Gannet and our developments in edited MRS methodology.

## Getting Started

### Prerequisites

Gannet runs in MATLAB (we recommend using the latest release if possible). Additionally, Gannet requires that the following MATLAB toolboxes are installed:

* Image Processing
* Optimization
* Signal Processing
* Statistics and Machine Learning

To run the voxel co-registration and structural image segmentation modules, [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) must be installed.

### Installing

The simplest way to install Gannet 3.0 is to download the zipped [master folder](https://github.com/richardedden/Gannet3.0/archive/master.zip), unzip it and then move the Gannet3.0-master folder to your MATLAB directory. To add the folder to your search path, run the following line of code in the command window of MATLAB:

macOS
```
addpath(genpath('~/Documents/MATLAB/Gannet3.0-master'));
```

Alternatively, open the Set Path dialog box from the MATLAB menu, click 'Add with Subfolders', find the Gannet3.0-master folder and then select it. When done, press 'Save' to permanently save the Gannet folder to MATLAB's search path.

## Built With

* [MATLAB](https://www.mathworks.com/products/matlab.html)

## Bugs, Contributions, Queries

For bug reporting, contribution requests and/or queries, please contact us at: gabamrs@gmail.com

## Versioning

Versioning is conducted on a module-specific basis using the style YYMMDD. That is, each Gannet module has its own release version. Major changes or addition of new functionalities are reflected in the version updates, while minor bug fixes are not.

## Authors

* **Richard Edden** (The Johns Hopkins University)
* **C. John Evans** (Cardiff University)
* **Ashley Harris** (University of Calgary)
* **Nicolaas Puts** (The Johns Hopkins University)
* **Georg Oeltzschner** (The Johns Hopkins University)
* **Muhammad Saleh** (The Johns Hopkins University)
* **Mark Mikkelsen** (The Johns Hopkins University)

## License

This software is open-source and does not have a specific license, but should you publish material that made use of Gannet in some way please cite the following publications:

* Edden RAE, Puts NAJ, Harris AD, Barker PB, Evans CJ. Gannet: A batch-processing tool for the quantitative analysis of gamma-aminobutyric acid-edited MR spectroscopy spectra. J. Magn. Reson. Imaging 2014;40:1445–1452. doi: [10.1002/jmri.24478](http://doi.wiley.com/10.1002/jmri.24478)
* Harris AD, Puts NAJ, Edden RAE. Tissue correction for GABA-edited MRS: Considerations of voxel composition, tissue segmentation, and tissue relaxations. J. Magn. Reson. Imaging 2015;42:1431–1440. doi: [10.1002/jmri.24903](http://doi.wiley.com/10.1002/jmri.24903) (if you report tissue-corrected, water-referenced measurements)

## Acknowledgments

We wish to thank the following individuals for their contributions to the development of Gannet and shared data processing code:

* Ralph Noeske (GE Berlin)
* Peter Barker (The Johns Hopkins University)
* Robin de Graaf (Yale School of Medicine)
* Philipp Ehses (Max Planck Institute for Biological Cybernetics, Tübingen)
* Wouter Potters (UMC Amsterdam)
