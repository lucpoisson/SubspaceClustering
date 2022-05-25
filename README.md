# SubspaceClustering

Repository for the paper [*Subspace clustering in high-dimensions: Phase transitions \& Statistical-to-Computational gap*](https://arxiv.org/abs/XXXX.XXXX). 

<p float="left">
  <img src="https://github.com/IdePHICS/RandomLabelsUniversality/blob/main/Figures/zero_reg.jpeg" height="270" />
</p>

*An illustration of the phase diagram for the clustering of two-classes sparse homogeneous Gaussian Mixture. We colour each region according to the associated *reconstruction phase*. The different phases are separated by different thresholds. 

## Structure

In this repository we provide the code and some guided example to help the reader to reproduce the figures of the paper [1]. The repository is structured as follows.

| File                          | Description                                                                                                                                                    |
|-------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------|
| ```/K2mix``` | Solver for the fixed point equations of the order parameters in the case of classification tasks on K=2 Gaussian clusters with nonhomogeneous diagonal covariances. The notebook [```how_to.ipynb```](https://github.com/gsicuro/GaussMixtureProject/blob/main/K2mix/how_to.ipynb) provides a step-by-step explanation on how to use the package. This implementation has been used to produce the results in Section 3.1 of the paper.           |
| ```/multiK``` | Solver for the fixed point equations of the order parameters in the case of classification tasks on K Gaussian clusters with diagonal covariances. The notebook [```how_to.ipynb```](https://github.com/gsicuro/GaussMixtureProject/blob/main/multiK/how_to.ipynb) provides a step-by-step explanation on how to use the package. This implementation has been used to produce the results in Section 3.2 of the paper.                                     |

The notebooks are self-explanatory.

## Reference

[1] *
