# Diffusion-based-Gaussian-Process
This is the code of the Graph Laplacian based Gaussian Process(GL-GP) algorithm and the datasets in the paper "GRAPH BASED GAUSSIAN PROCESSES ON RESTRICTED DOMAINS" by David B Dunson, Hau-Tieng Wu and Nan Wu

Explanation to files 

1.CovMatrix.m

The Matlab code of the GL-GP for predictions.

2. Nys.m

The code of Nystrom type extension. Suppose the GL-GP covariance parameters are determined by a small number of labeled points, the corresponding response variables, and a small number of unlabelled points. Nys.m can extend the prediction to a large number of unlabel points without determing the GL-GP paramters again. 

2. hand.mat

The dataset in Raynauld disease example. 

3. hand.txt

The explanation of the variables in the dataset "hand.mat".

4. spiral example.txt

The Matlab code to generate a random dataset and a regression function on a spiral in R^2.

5. spiral.mat

The dataset on the spiral that we use in the paper.

6. two spheres example.txt

The Matlab code to generate a random dataset and a regression function on the spheres example in R^3.

7. two spheres.mat

The dataset on the two spheres example that we use in the paper.
