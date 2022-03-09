# Confidence Composition: Case Studies
The data and analysis scripts for case studies on Confidence Composition.

Full paper: [arXiv](https://arxiv.org/abs/2111.03782). Feel free to contact the first author with any questions. 

### Folders

* coco-case-studies/: source code for the analysis of the case studies. 
* coco-case-studies/mountaincar-logs: the mountain car case study from OpenAI Gym. 
* coco-case-studies/uuv-logs: the unmanned underwater vehicle (UUV) pipeline-scanning case study based on the ROS UUV simulator. 

### Installing Wolfram Script

The goal of the analysis is to reproduce the performance measurements of confidence monitors as reported in Tables 1 and 2 in the paper references above.

To run the analysis scripts, you would need to dowload and install a free version of the Wolfram Engine: 

1) Download and install the [Wolfram Engine](https://www.wolfram.com/engine/) for your operating system. 
2) Create a Wolfram ID account (if you do not have one). 
3) Accept the terms of the [free license](https://www.wolfram.com/engine/free-license)
4) Run "wolframscript" in your command line
5) Input your Wolfram ID and password into the prompt

### Executing the analysis via scripts 

You can follow the Docker-based repeatability instructions in the [reproduction-instruction.pdf](./reproduction-instruction.pdf). Follow the build-image-from-scratch route or download the built [docker image](https://drive.google.com/drive/folders/1YY6M_aFl7YsUhBvdzMZTTsf-8UrVAN2e).

To run the analysis scripts: 
* UUV case study
   - Neutral calibration: wolframscript <path_to>/uuv-analysis-neutral.wls
   - Conservative calibration: wolframscript <path_to>/uuv-analysis-conservative.wls
* Mountain car case study 
   - Neutral calibration: wolframscript <path_to>/mountaincar-analysis-neutral.wls
   - Conservative calibration: wolframscript <path_to>/mountaincar-analysis-conservative.wls

"wolframscript" is equivalent to and can be replaced with "wolfram -script". 

### Interactive analysis via notebooks 

The notebook-based analysis allows interactive processing and visualization of the data. You would need to install [Wolfram Mathematica](https://www.wolfram.com/mathematica/), version 12.1+. You can obtain a [free trial](https://www.wolfram.com/mathematica/trial/) or you may be at an institution with an [academic license](https://www.wolfram.com/mathematica/pricing/colleges-universities/). For read-only access, you can download a free [Wolfram Player](https://www.wolfram.com/player/).

The following notebooks are provided: 
* composition-utils.nb: utility functions to analyze compositions of confidence monitors. It contains the functionality analogous to the composition-utils-package.wl package. 
* mountaincar-analysis.nb: a notebook with the analysis and exploration of the data from the mountain car case study. 
* uuv-analysis.nb: a notebook with the analysis and exploration of the data from the underwater vehicle case study. 

