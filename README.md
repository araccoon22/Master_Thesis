# 🧠 Master Thesis – Functional Summaries of Persistent Homology for Neuronal Morphologies

This repository contains the code, figures, and analysis conducted as part of my Master’s thesis in Applied Mathematics at EPFL, under the supervision of Kathryn Hess Bellwald and Lida Kanari.

The project investigates the use of **persistent homology** and **functional summaries** to detect and quantify topological differences in neuronal morphologies, with a focus on Layer 2/3 pyramidal neurons from human and mouse cortex.

## 📂 Repository Structure

├── Images/                             # Figures used in the report
├── Training_Scripts/                  # Helper scripts for functions and automation
├── persistence_landscape_basal_axon/  # Full analysis for basal dendrites and axons
├── 0_check_data.ipynb                 # Data loading and consistency checks
├── 1_entropy.ipynb                    # Entropy summary analysis
├── 2_persistence_landscape.ipynb      # Persistence landscape summary analysis
├── 3_persistence_silhouette.ipynb     # Persistence silhouette summary analysis
├── 4_permutation_test.ipynb           # Permutation tests for statistical validation
├── functions.py                       # Shared utility functions
├── PDM_Report.pdf                     # Final Master’s thesis report
├── .gitignore                         # Ignored files and directories
└── README.md                          # Project documentation


## 🧪 Project Goals

This project answers four key questions:
1. Do human and mouse pyramidal neurons differ topologically in a statistically meaningful way?
2. Can functional summaries (entropy, landscape, silhouette) capture these differences?
3. Are the differences robust with respect to intra-group variability?
4. Which descriptors are most effective for statistical inference?

## 🔍 Methods

We use the **Topological Morphology Descriptor (TMD)** to compute persistence diagrams from neuronal reconstructions and analyze them using the following functional summaries:
- **Persistence Entropy**: Scalar measure of topological disorder.
- **Persistence Landscape**: Functional summary enabling rich statistical comparison.
- **Persistence Silhouette**: Weighted average function focusing on high-persistence features.

We apply **non-parametric permutation tests** and construct **bootstrap confidence bands** to validate statistical significance.

## 📊 Key Findings

| Descriptor             | Discriminative Power | Stability | Interpretability |
|------------------------|----------------------|-----------|------------------|
| Persistence Entropy    | Low                  | High      | High             |
| Persistence Landscape  | High                 | High      | Medium           |
| Persistence Silhouette | Medium               | Moderate  | Medium           |

- The **persistence landscape** consistently provided the strongest performance.
- Significant inter-species differences were detected in apical and basal dendrites.
- Axonal differences were not statistically significant, likely due to reconstruction noise and sample size limitations.

## 🚀 Future Work

One promising direction is to explore **generalized persistence landscapes** (Berry et al., 2020). These allow greater flexibility in kernel design and may further improve the discriminative power of functional summaries.


## 🔗 Links
- 📦 [TMD Documentation](https://tmd.readthedocs.io/en/latest/)
- 🧮 [Giotto-TDA Library](https://github.com/giotto-ai/giotto-tda)




