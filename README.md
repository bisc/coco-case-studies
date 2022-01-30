# Confidence Composition: Case Studies
The data and analysis scripts for case studies on Confidence Composition.

Full paper: [arXiv](https://arxiv.org/abs/2111.03782)

### Folders

* mountaincar-logs: the mountain car case study from OpenAI Gym. 
* uuv-logs: the unmanned underwater vehicle (UUV) pipeline-scanning case study based on the ROS UUV simulator. 

### Installing Wolfram Script

The goal of the analysis is to reproduce the performance measurements of confidence monitors as reported in Tables 1 and 2 in the paper references above.

To run the analysis scripts, you would need to dowload and install a free version of the Wolfram Engine: 

1) Download and install the [Wolfram Engine](https://www.wolfram.com/engine/) for your operating system. 
2) Create a Wolfram ID account (if you do not have one). 
3) Accept the terms of the [free license](https://www.wolfram.com/engine/free-license)
4) Run "wolframscript" in your command line
5) Input your Wolfram ID and password into the prompt

### Executing the analysis

To run the analysis scripts: 
* UUV case study
   - Neutral calibration: wolframscript uuv-analysis-neutral.wls
   - Conservative calibration: wolframscript uuv-analysis-conservative.wls
* Mountain car case study , 
   - Neutral calibration: wolframscript mountaincar-analysis-neutral.wls
   - Conservative calibration: wolframscript mountaincar-analysis-conservative.wls

"wolframscript" is equivalent to and can be replaced with "wolfram -script". 
