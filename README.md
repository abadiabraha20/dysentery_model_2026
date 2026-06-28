# dysentery_model_2026

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

**Python implementation of an SEIRS-B dysentery transmission model with environmental reservoir, stratified susceptibility, optimal control, and serotype-specific intervention strategies. Includes code for simulation, bifurcation analysis, sensitivity analysis (PRCC/LHS), optimal control (forward-backward sweep), and figure generation. Reproduces all results in the manuscript.**

---

## 📚 Overview

This repository contains the complete code, data, and supplementary materials for the paper:

> **"Global Dynamics, Optimal Control, and Serotype-Specific Intervention Strategies for Dysentery Transmission with Environmental Reservoir and Stratified Heterogeneity"**
>
> *Abadi Abraha Asgedom, Yohannes Yirga Kefela, Hailu Tkue Welu*
>
> *Department of Mathematics, Mekelle University, Mekelle, Ethiopia*

The model incorporates an SEIRS-B framework with an exposed compartment, logistic pathogen growth in an environmental reservoir, and stratified susceptibility to capture social heterogeneity in exposure risk. We derive the basic reproduction number $\mathcal{R}_0$ and prove global stability of disease-free and endemic equilibria using Lyapunov functions and geometric methods. The model is calibrated to Ethiopian outbreak data (2015-2019) with $R^2=0.94$.

---

## 🎯 Key Features

| Feature | Description |
|---------|-------------|
| **SEIRS-B Model** | Complete implementation with environmental reservoir and logistic pathogen growth |
| **Stratified Susceptibility** | High-risk and low-risk population modeling with 22% R₀ effect |
| **Optimal Control** | Forward-Backward Sweep algorithm with 43% cost savings |
| **Serotype Analysis** | *S. flexneri* (R₀=2.34) and *S. sonnei* (R₀=0.91) parameter sets |
| **Data Calibration** | Nonlinear least squares with Latin Hypercube sampling (R²=0.94) |
| **Sensitivity Analysis** | PRCC with 95% confidence intervals (10,000 LHS samples) |
| **Bifurcation Analysis** | Forward bifurcation at R₀=1 (Theorem 6) |
| **Global Stability** | Lyapunov and geometric methods (Theorems 4-5) |

---

## 📁 Repository Structure

```
dysentery_model_2026/
│
├── README.md                               # This file
├── LICENSE                                 # MIT License
├── requirements.txt                        # Python dependencies
├── Figures.py                              # Complete python code for all figures
├── data/
│   ├── ethiopia_outbreak.csv              # Ethiopian outbreak data (2015-2019)
│   └── parameters.yaml                    # Default parameter sets
│
├── figures                                # Main figures (Fig1-Fig12)
│
├── supplementary/
│   ├── Supplementary_Materials.tex        # LaTeX source
│   └── Supplementary_Materials.pdf        # Compiled PDF
```

---

## 🚀 Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/abadiabraha20/dysentery_model_2026.git
cd dysentery_model_2026

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Install package in development mode
pip install -e .
```

### Requirements

```txt
numpy>=1.21.0
scipy>=1.7.0
matplotlib>=3.4.0
pandas>=1.3.0
seaborn>=0.11.0
pyyaml>=5.4.0
tqdm>=4.62.0
SALib>=1.4.0
```

### Basic Usage

```python
import numpy as np
from src.model import seirsb_model, compute_R0, get_default_params
from scipy.integrate import odeint

# Load default parameters
params = get_default_params()

# Compute basic reproduction number
R0 = compute_R0(params)
print(f"Basic Reproduction Number: {R0:.4f}")

# Run simulation
t_span = np.linspace(0, 200, 2000)
y0 = [1000, 50, 10, 5, 100]
sol = odeint(seirsb_model, y0, t_span, args=(params,))

# Extract infected population
I_t = sol[:, 2]
```

---

## 📊 Generating Figures

### Main Figures (Fig1-Fig12)

```bash
# Generate all main figures
python scripts/generate_figures.py

# Generate specific figure
python scripts/generate_figures.py --figure 2  # Bifurcation diagram
python scripts/generate_figures.py --figure 3  # Optimal control
python scripts/generate_figures.py --figure 9  # R0 < 1 time series
```

### Supplementary Figures (S1-S5)

```bash
# Generate all supplementary figures
python scripts/generate_supplementary.py

# Generate specific supplementary figure
python scripts/generate_supplementary.py --figure S1  # Sigma bifurcation
python scripts/generate_supplementary.py --figure S4  # Joint sensitivity
```

---

## 🔬 Model Description

### SEIRS-B Model Equations

The model system is given by:

```
dS/dt = Λ + γ(1-ρ₁)I + ε₁u₁(1-ρ₂)I + αR - (β₁B/(K+B) + β₂I + μ)S
dE/dt = (β₁B/(K+B) + β₂I)S - (μ + κ)E
dI/dt = κE - (μ + d + γ + ε₁u₁)I
dR/dt = (γρ₁ + ε₁u₁ρ₂)I - (μ + α)R
dB/dt = εI + rB(1 - B/Bmax) - (σ + ε₂u₂)B
```

### Reproduction Numbers

**Basic Reproduction Number (no controls):**

```
R₀ = β₂κΛ/[μ(κ+μ)(μ+d+γ)] + β₁εκΛ/[μK(κ+μ)(μ+d+γ)(σ-r)]
```

**Control Reproduction Number (with interventions):**

```
Rc = β₂κΛ/[μ(κ+μ)(μ+d+γ+ε₁u₁)] + 
     β₁εκΛ/[μK(κ+μ)(μ+d+γ+ε₁u₁)(σ+ε₂u₂-r)]
```

### Key Parameters

| Parameter | Description | Value/Range |
|-----------|-------------|-------------|
| Λ | Recruitment rate | 0.1--100 |
| μ | Natural mortality | 0.00045662 day⁻¹ |
| β₁ | Environmental transmission | 10⁻⁵--10⁻³ |
| β₂ | Direct transmission | 10⁻⁶--10⁻⁴ |
| κ | Progression rate | 0.2--0.5 day⁻¹ |
| ε | Pathogen shedding | 0.1--1.0 |
| σ | Net decay | 0.033333 day⁻¹ |
| r | Pathogen growth | 0--0.02 day⁻¹ |
| α | Waning immunity | 0.002--0.01 day⁻¹ |

---

## 🧪 Running Simulations

### Model Simulation

```bash
# Run default simulation
python scripts/run_simulation.py

# Run with custom parameters
python scripts/run_simulation.py --lambda 0.5 --beta1 0.000075 --beta2 0.000011

# Run with controls
python scripts/run_simulation.py --u1 0.05 --u2 0.1
```

### Optimal Control

```bash
# Run optimal control
python scripts/run_optimal_control.py

# Run with custom weights
python scripts/run_optimal_control.py --A1 100 --A2 50 --B1 10 --B2 5
```

### Data Calibration

```bash
# Calibrate model to Ethiopian outbreak data
python scripts/run_calibration.py --data data/ethiopia_outbreak.csv

# Run with custom iterations
python scripts/run_calibration.py --iterations 10000
```

### Sensitivity Analysis

```bash
# Run PRCC sensitivity analysis
python scripts/run_sensitivity.py --samples 10000

# Run with specific parameters
python scripts/run_sensitivity.py --parameters beta1 beta2 epsilon
```

---

## 📈 Figure Descriptions

### Main Figures

| Figure | Description | Reference |
|--------|-------------|-----------|
| **Fig 1** | Transmission diagram of SEIRS-B model | Section 2 |
| **Fig 2** | Forward bifurcation at R₀=1 (Theorem 6) | Section 5 |
| **Fig 3** | Optimal control strategies (Theorem 7) | Section 6 |
| **Fig 4** | Calibration to Ethiopian data (R²=0.94) | Section 7 |
| **Fig 5** | Global sensitivity analysis (PRCC) | Section 8 |
| **Fig 6** | Serotype-specific comparison | Section 8 |
| **Fig 7** | Stratified susceptibility effect (22% increase) | Section 8 |
| **Fig 8** | Controllable region heatmap | Section 8 |
| **Fig 9** | Disease elimination (R₀<1, Theorem 4) | Section 9 |
| **Fig 10** | Disease persistence (R₀>1, Theorem 5) | Section 9 |
| **Fig 11** | Intervention efficacy | Section 9 |
| **Fig 12** | With vs without controls | Section 9 |

### Supplementary Figures

| Figure | Description | Reference |
|--------|-------------|-----------|
| **S1** | Bifurcation for varying σ-r | S3.1 |
| **S2** | Bifurcation for varying β₂ | S3.2 |
| **S3** | Bifurcation for varying K | S3.3 |
| **S4** | Joint sensitivity β₁, β₂ | S5.3.1 |
| **S5** | Joint sensitivity ε, σ | S5.3.2 |

---

## 📝 Key Results

### Serotype-Specific Thresholds

| Serotype | R₀ | Required Intervention |
|----------|-----|----------------------|
| *S. flexneri* | 2.34 | u₁ ≥ 0.25 + u₂ = 0.30 |
| *S. sonnei* | 0.91 | u₂ = 0.15 (sanitation alone) |

### Optimal Control Savings

| Metric | Value |
|--------|-------|
| Infection reduction vs no control | 78% |
| Infection reduction vs constant controls | 43% |
| Cost savings | 43% |

### Calibration Results

| Metric | Value |
|--------|-------|
| R² | 0.94 |
| RMSE | 127.3 |
| Calibrated R₀ | 2.18 (95% CI: 1.92-2.44) |

---

## 🧪 Testing

```bash
# Run all tests
pytest tests/

# Run specific test
pytest tests/test_model.py

# Run with coverage
pytest --cov=src tests/
```

---

## 🤝 Contributing

We welcome contributions! Please follow these steps:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

Please ensure:
- All tests pass
- Code is properly documented
- New features include tests

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 📧 Contact

**Abadi Abraha Asgedom**  
*Department of Mathematics, Mekelle University, Mekelle, Ethiopia*  
Email: abadi.abraha@mu.edu.et

**Yohannes Yirga Kefela**  
Email: yohannes.yirga@mu.edu.et

**Hailu Tkue Welu**  
Email: hailu.tkue@mu.edu.et

---

## 📚 Citation

If you use this code in your research, please cite:

```bibtex
@article{asgedom2026dysentery,
  title={Global Dynamics, Optimal Control, and Serotype-Specific Intervention Strategies for Dysentery Transmission with Environmental Reservoir and Stratified Heterogeneity},
  author={Asgedom, Abadi Abraha and Kefela, Yohannes Yirga and Welu, Hailu Tkue},
  journal={Journal of Biological Dynamics},
  year={2026},
  volume={XX},
  number={X},
  pages={XX-XX},
  doi={10.XXXX/XXXXXX}
}
```

---

## 📖 References

1. Berhe, H.W., Makinde, O.D., Theuri, D.M. (2021). Stability analysis and optimal control of dysentery diarrhea epidemic model. *Journal of Biological Dynamics*, 15(1), 194-216.

2. Kotloff, K.L., Riddle, M.S., Platts-Mills, J.A., et al. (2018). Shigellosis. *The Lancet*, 391(10122), 801-812.

3. Li, M.Y., Muldowney, J.S. (1996). A geometric approach to global-stability problems. *SIAM Journal on Mathematical Analysis*, 27(4), 1070-1083.

4. Lenhart, S., Workman, J.T. (2007). *Optimal Control Applied to Biological Models*. Chapman & Hall/CRC.

5. Marino, S., Hogue, I.B., Ray, C.J., Kirschner, D.E. (2008). A methodology for performing global uncertainty and sensitivity analysis in systems biology. *Journal of Theoretical Biology*, 254(1), 178-196.

6. Castillo-Chavez, C., Song, B. (2002). Dynamical models of tuberculosis and their applications. *Mathematical Biosciences and Engineering*, 1(2), 361-404.

---

## 🙏 Acknowledgements

The authors acknowledge institutional support from Mekelle University and thank the anonymous reviewers for their valuable comments. No specific funding was received for this research.

---

**Last Updated:** June 2026
