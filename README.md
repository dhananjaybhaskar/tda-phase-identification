#### Summary

This public repository contains code to reproduce results published in [Bhaskar et al., Soft Matter, 2021](https://doi.org/10.1039/D1SM00072A "DOI: 10.1039/D1SM00072A"). The paper demonstrates that topology-based machine learning classifies complex spatial patterns exhibited by epithelial cells into distinct phases. The presence and stability of closed loops, obtained from Vietoris-Rips filtration applied to cell positions, can effectively distinguish between individual and collective phases, even when population size varies significantly due to proliferation.

---

#### Organization

| Source Code | Description |
|--|--|
| `model_param_sweep.m` | Simulate self-propelled particle model. |
| `density_map.m` | Plot particle density order parameter. |
| `Representative Barcode Statistics.ipynb` | Statistics for number and size of persistent loops in individual, branching and clustered phases. |
| `Pairwise Wasserstein Distances.ipynb` | Compute pairwise Wasserstein distances using persistence diagrams obtained from particle positions. |
| `Phase Diagram Classification.ipynb` | Permute pairwise Wasserstein distance matrix using hierarchical clustering and map classification results to phase diagram. |
| `Robustness (Julia Computation).ipynb` and `Robustness (Python Visualization).ipynb` | Quantify robustness of classification w.r.t. sparse sampling, spatial scaling and over time. |
| `Experimental Data Analysis.ipynb` | Compute persistent homology of cell nuclei positions and compare with simulations. |
