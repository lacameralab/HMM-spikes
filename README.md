# HMM for spike data

Author: Tianshu Li

The HMM-spikingdata toolbox contains MATLAB source code to train and decode several HMM algorithms to spike data. 
Full details about the algorithms are in the article:

T. Li and G. La Camera, *A sticky Poisson Hidden Markov Model for solving the problem of over-segmentation and rapid state switching in cortical datasets*, bioRxiv 2024.08.07.606969; doi: https://doi.org/10.1101/2024.08.07.606969

MATLAB version: MATLAB 2023a

## Folders and files

- src/ contains all the source code.

- src/demo/ contains four demos:
  standard Poisson HMM (standardPHMMdemo.m)\
  sticky-Poisson HMM (stickypHMMdemo.m)\ 
  Poisson HMM with Dirichlet prior (DPHMM, DirichletPHMMdemo.m)\
  Multinoulli HMM (mHMMdemo.m).\ 
  In each case, the demo code shows an example of traning and decoding the corresponding HMM algorithm on a dataset 
  generated by a clustered spiking network. The same algorithm can be applied to other datasets by changing 
  the data folder and the relevant parameters (see demo file headers).
    
- src/pHMM/ contains functions for Poisson HMM.
    
- src/mHMM/ contains functions for Multinoulli-HMM.
    
- src/utils/ contains other functions used in the demo, including functions used for plotting. 
  AIC.m and BIC.m are provided in this folder. They are not used in demo, but are needed for model selection.

- data/ contains the example data used in the demos.
    

### Comments
    
- Detailed instructions are in the header of the demo files.
    
- PHMM does not require the length of each trial to be the same. 
    The demo uses the data with trials having the same length just for simplicity. 
    If you want to apply the toolbox to trials with different lengths, use 
    ./src/utils/plotpHMMtrials_varT.m instead of ./src/utils/plotpHMMtrials.m to plot the decoding 
    of all trials in one figure.
