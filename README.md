# Wake Model in ABL Flows

This repository provides a MATLAB implementation for computing the steady-state wake structure behind yawed or unyawed wind turbines operating in Atmospheric Boundary Layer (ABL) flows. The model incorporates the effects of wind veer and thermal stratification present in ABL environments.

The extended wake model integrates a Gaussian wake formulation for yawed wind turbine wakes with an ABL model that self-consistently predicts the three-dimensional Ekman wind profile under both neutral and stable stratification. This coupling enables accurate prediction of wake behavior across a broad range of atmospheric conditions.

The model's predictions have been validated against Large Eddy Simulation (LES) data, showing good agreement and confirming its capability to capture key wake characteristics under various atmospheric conditions.

The implementation is based on the following peer-reviewed publication:

> **Ghanesh Narasimhan**, Dennice F. Gayme, Charles Meneveau (2025).  
> *An extended analytical wake model and applications to yawed wind turbines in atmospheric boundary layers with different levels of stratification and veer*.  
> *Journal of Renewable and Sustainable Energy*, **17**(3), 033302.  
> [https://doi.org/10.1063/5.0251305](https://doi.org/10.1063/5.0251305)

Full-text also available on ResearchGate:  
[https://www.researchgate.net/publication/392165255_An_extended_analytical_wake_model_and_applications_to_yawed_wind_turbines_in_atmospheric_boundary_layers_with_different_levels_of_stratification_and_veer](https://www.researchgate.net/publication/392165255_An_extended_analytical_wake_model_and_applications_to_yawed_wind_turbines_in_atmospheric_boundary_layers_with_different_levels_of_stratification_and_veer)

---

## 📄 Files

- `extended_wake_model.mlx` – Main MATLAB Live Script containing the full wake model implementation. Run this script to compute the wake shape for yawed or unyawed wind turbines placed within neutral or stable ABL flows.

---

## ✅ Required Inputs

Define the following input parameters in the main script:

| Parameter     | Description                                                  | Units     |
|---------------|--------------------------------------------------------------|-----------|
| `G`           | Geostrophic wind speed                                       | m/s       |
| `phi`         | Latitude                                                     | degrees   |
| `z0`          | Surface roughness length                                     | m         |
| `lapse_rate`  | Free-stream potential temperature lapse rate                 | K/m       |
| `Cr`          | Cooling rate of surface potential temperature                | K/s       |
| `Theta0`      | Initial surface potential temperature                        | K         |
| `D`           | Rotor diameter of the turbine                                | m         |
| `yaw_angle`   | Yaw angle of the wind turbine                                | degrees   |
| `Ct`          | Local thrust coefficient                                     | —         |
| `zt`          | Hub height of the wind turbine                               | m         |

---

## 📌 Notes

- For questions or collaborations, please contact Dr. Ghanesh Narasimhan at naras062@umn.edu.




