
# sPCF-EM: Switching Poisson Cubature Filter - Expectation Maximization <br/>


Learn parameters for a switching dynamical system with Poisson observations using a switch EM framework with a Poisson Cubature Filter (PCF) embedded.


# Publication
For the derivation of PCF and how it was embedded in both stationary and switching EM frameworks, see the paper below. Validation of the method on both simulated and experimental data is also shown in the paper with analysis of one of the sessions from the experimental dataset (publicly available from the Sabes lab) given as an example here.

Song C.Y., Shanechi M.M., "Unsupervised learning of stationary and switching dynamical system models from Poisson observations", Journal of Neural Engineering, Dec. 2023

Link to the paper: https://doi.org/10.1088/1741-2552/ad038d

#

For details on inference of switching dynamical systems (estimating behavior and regimes from neural data given learned model parameters) see the below publication: 

Song C.Y. et al, "Modeling and Inference Methods for Switching Regime-Dependent Dynamical Systems with Multiscale Neural Observations",  Journal of Neural Engineering, Oct. 2022

Link to paper: https://doi.org/10.1088/1741-2552/ac9b94

Original preprint: https://doi.org/10.1101/2022.06.09.494416

Summary: https://twitter.com/MaryamShanechi/status/1536343291347668992



# Usage guide
## Dependencies
Code was developed and tested on MATLAB R2021a and R2023a. The following packages are also used:
- Optimization Toolbox
- Statistics and Machine Learning Toolbox

## Initialization
Add the source directory and its subdirectories to the path. Adding the experimental directory is optional unless you want to run the example script. You can run init.m to do this. 

## Main learning function
The main learning function [source/fitSwitchEM.m](source/fitSwitchEM.m). A complete usage guide is available in the function. Also see [source/example scripts/scriptLearnSabes.m](source/example%20scripts/scriptLearnSabes.m). The following shows an example case: 
```
trnS = struct('fs',fs,'obsNt',nt);
prmS = struct('dimXtEst',nx,'dimStEst',ns);
mdlS = fitSwitchEM(trnS,prmS,nIter);
```
Inputs:
- nt is the input to trnS.obsNt which is a dimension x time matrix with spiking activity observation data. Can also be a cell array of dimension x varying time length matrices for trial based data.
- fs is the input to trnS.fs which is the sampling frequency in Hz.
- nx is the input to prmS.dimXtEst which is the assumed latent state xt dimension.
- ns is the input to prmS.dimStEst which is the assumed number of regimes.
- nIter is the number of iterations to run EM.

Output:
- mdlS: a structure containing learned model parameters and any tracked metrics. The final parameters can be obtained by mdlS.thetaCell{end}. See code for details.

## Estimating latent and regime states in test data
Once a model is learned, you can apply the model to new data to estimate latent and regime states. The latent states can also be projected onto some desired behavior to estimate behavior as well. See [source/example scripts/scriptLearnSabes.m](source/example%20scripts/scriptLearnSabes.m) for a more detailed example and [source/filterStates.m](source/filterStates.m) for more details on outputs. For acausal estimation (smoothing), see [source/smoothStates.m](source/smoothStates.m).
```
theta = mdlS.thetaCell{end};
mngr = filterStates(theta,ntTest);
PStDec = mngr.PStEst;
```
Input:
- theta: model parameters
- ntTest: neural observations dimension x time.

Outputs:
- mngr: structure containing filtered states and estimation covariances
- mngr.PStDec: decoded probability of regimes. Choose regime with highest probability to get decoded regimes.
- mngr.xDec: decoded latent states xt. You can fit a projection L onto behavior in training and then apply it on test L*mngr.xDec to get decoded behavior. See the example script.


# Usage examples
An example script analyzing a publicly available session from the Sabes lab which was also performed in the above paper is provided in [source/example scripts/scriptLearnSabes.m](source/example%20scripts/scriptLearnSabes.m). The script learns a model from Poisson observations and then performs causal decoding of a held-out test set to estimate regimes as in Figure 8.


# Licence
Copyright (c) 2023 University of Southern California  
See full notice in [LICENSE.md](LICENSE.md)  
Christian Y. Song and Maryam M. Shanechi  
Shanechi Lab, University of Southern California