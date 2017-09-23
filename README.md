## Synopsis

Projected Power Method (PPM) tries to infer the values of variables given pairwise noisy modulo measurements. It has applications to 3D object alignments that can provide better oriented training examples for convolutional neural networks.
Even though this framework works for modulo measurements, it does not provide the correct initialization for extending the algorithm to community detection. In this project, we have mathematically determined the correct initialization that can help the algorithm start from the basin of attraction and eventually reach a global optima, even though the overall problem is non-convex. 

## Motivation

The problem of robust community detection in censored block model is important in itself. Moreover, we have some concrete applications in genomics where this framework can be used. 

## Robustness

Apart from the community detection in friendship-foe networks, our algorithm also incorporates robustness, where outlier nodes with random connections can be detected and eliminated in subsequent iterations. The robust algorithm is 
motivated from robust regression where there is a significant fraction of outliers (<= 0.4), yet theoretical guarantees exist for identifying the right set of points.

## Installation

Download the git folder and run `main_CBM_Adverse.m` in MATLAB. 


## Tests

An example run is as follows :

```
Running Robust Alg -->
PPM without Robustness , Error = 0.645000   False positive = (45 / 300 )
Iteration 1 , Error = 0.385000   False positive = (13 / 300 )
Iteration 2 , Error = 0.030000   False positive = (6 / 300 )
Iteration 3 , Error = 0.030000   False positive = (6 / 300 )
Iteration 4 , Error = 0.035000   False positive = (6 / 300 )
Iteration 5 , Error = 0.035000   False positive = (6 / 300 )
Iteration 6 , Error = 0.035000   False positive = (6 / 300 )
Iteration 7 , Error = 0.035000   False positive = (6 / 300 )
Iteration 8 , Error = 0.255000   False positive = (9 / 300 )
Iteration 9 , Error = 0.075000   False positive = (7 / 300 )
 Iterative Rebustness Algorithm Reduces Error !!!
 ```

 
