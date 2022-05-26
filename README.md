# SubspaceClustering  (:warning: under construction! :construction: :construction_worker:)

Repository for the paper [*Subspace clustering in high-dimensions: Phase transitions \& Statistical-to-Computational gap*](https://arxiv.org/abs/XXXX.XXXX). 

<p float="left">
  <img src="https://github.com/lucpoisson/SubspaceClustering/blob/main/Figures/subspace2.png" height="470" />
</p>
An illustration of the phase diagram for the clustering of two-classes sparse homogeneous Gaussian Mixture. We colour each region according to the associated *reconstruction phase* . The different phases are separated by the associated *thresholds* . 

## Structure

In this repository we provide the code and some guided example to help the reader to reproduce the figures of the paper [1]. The repository is structured as follows.

| File                          | Description                                                                                                                                                    |
|-------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------|
|/AMP_SE| Solver for AMP and its State Evolution (SE) equations. The notebook [```/AMP_SE/how_to_AMP/SE.ipynb```] provides a step-by-step explanation on how to use the package. This implementation has been used to produce the results in Sections 4/5/6 of the paper.           |
| /Miscellaneous | Solver for the general-purpose algorithms PCA, sparse PCA and Diagonal Thresholding. The notebook [```/Miscellaneous/how_to_mix.ipynb```] provides a step-by-step explanation on how to use the package. This implementation has been used to produce the results in Sections 4/5/6 of the paper.                  |


The notebooks are self-explanatory.

## Reference

[1]*Subspace clustering in high-dimensions: Phase transitions \& Statistical-to-Computational gap*,
Luca Pesce, Bruno Loureiro, Florent Krzakala, Lenka Zdeborová [arXiv: XXXX.XXXX](https://arxiv.org/abs/XXXX.XXXX)[stat.ML]

