This repository contains the MATLAB scripts used for the paper Qin, S., Li, Q., Tang, C., & Tu, Y. (2019). Optimal compressed sensing strategies for an array of nonlinear olfactory receptor neurons with and without spontaneous activity. Proceedings of the National Academy of Sciences, 116(41), 20286-20295.

Given the statistics of an odor space, characterized by sparse odorants and wide range of concentration, the goal is to understand the optimal strategies of the odor receptor neurons. 

CMA-ES was used to find the optimal sensitivity matrix of the odor receptors for neurons with and without spontaneous firing. To numerically estimate the mutual information, we used 

The code is tested on MATLAB 2018a, 2020b

# Dependence
- The `CMA-ES` package can be downloaded [here](https://cma-es.github.io/cmaes_sourcecode_page.html)
- `Information theoretical estimators toolbox` can be downloaded [here](https://bitbucket.org/szzoli/ite/src/master/)
- `Gaussian copula mutual information` can be downloaded [here](https://github.com/robince/gcmi)

# Simulation
## Log-normal concentration distribution
`optMatrixCMA_v7.m` is the main script used to find the optimal sensitivity matrix. Basic usage: `[wmin,fmin] = optMatrixCMA_v7(numOdor,numRecp,sparsity,sig)` 
- `numOdor` is the total number of possible odorants, default value 100
- `numRecp` is the number of OSN, default 30. 
- `sparsity` the number of average odorants appear in an odor, default 2.
- `Cons`  a boolean value indicating whether constraints are applied to the input, for example the correlation
- `r0`    the relative spontaneous firing, 0 ~ 1, by default it is 0
- `sig` the standard deviation of the log concentration. 
- `wmin` the optimal values of the sensitivity matrix
- `fmin` the value of the objective function
