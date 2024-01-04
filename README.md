# HMM-spikingdata
---
    author: Tianshu Li
    topic: description of the HMM-spikingdata toolbox
---

# Notes for HMM toolbox

The HMM-spikingdata toolbox contains the source code used in the article:

*T. Li and G. La Camera, "Hidden Markov Modeling of spike data"*


## Folders and files

    src/ contains all the source code.

    src/demo/ contains four demos: standard Poisson HMM (standardPHMMdemo.m), sticky-Poisson HMM (stickypHMMdemo.m), 
    Poisson HMM with Dirichlet prior (DPHMM, DirichletPHMMdemo.m), Multinoulli HMM (mHMMdemo.m). 
    The demos showed an HMM procedure of a dataset generated from a clustered spiking network, 
    including generating initial parameters, training and decoding. One can apply the algorithms to other datasets 
    by changing the directory of data and parameters according to the instructions.
    
    src/pHMM/ contains functions for Poisson HMM.
    
    src/mHMM/ contains functions for Multinoulli-HMM.
    
    src/utils/ contains other functions used in the demo, including functions used for plotting. 
    AIC.m and BIC.m are provided in this folder. They are not used in demo, but are needed for model selection.

    data/ contains the example data used in the demo.
    

## Comments:
    
    Detailed instructions are in demo files.
    
    PHMM does not require the length of each trial to be the same. 
    The demo uses the data with trials having the same length just for simplicity. 
    If one wants to apply the toolbox to trials with different lengths, please use 
    ./src/utils/plotpHMMtrials_varT.m instead of ./src/utils/plotpHMMtrials.m to plot the decoding of all trials in one figure.
