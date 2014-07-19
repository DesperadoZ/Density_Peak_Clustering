"Density Peak Clustering" a Matlab implementation
==================

This is a practical MATLAB implementation of the density-based clustering algorithm described in [1].

The current implementation includes a sample dataset: the sEMG data sampled from both healthy and neurologically-impaired subjects(stroke patients). 10 Features without PCA. The last column of the feature matrix is the Brunnstrom Classification labels which can be used for performance evaluation (data with label 6 are sampled from healthy subject). 

Class 'nodez' was created to store structured sample/observations in various applications. In this project, it has been modified to suit density-based algorithm. The input dataset should be a M by N feature matrix which consists of M observations/samples and N features. It is recommended that the input dataset to contain only 1 feature matrix to avoid confusion.

[1] Rodriguez, A. and A. Laio (2014). "Clustering by fast search and find of density peaks." Science 344(6191): 1492-1496.

