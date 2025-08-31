# Four-VSCMG Spacecraft Simulator (Pyramid Configuration)

This repository provides a MATLAB-based simulation of a spacecraft equipped with four Variable Speed Control Moment Gyroscopes (VSCMGs) in a pyramid configuration. The simulator models nonlinear spacecraft attitude dynamics, gimbal angles, and wheel speeds, and includes a 3D visualization of spacecraft motion.

---
<img width="512" height="344" alt="f2" src="https://github.com/user-attachments/assets/d506201f-a04d-4c41-b7cb-e49df836d754" />

## 🚀 Features
- Nonlinear spacecraft attitude dynamics  
- Explicit VSCMG actuator modeling (gimbal angles & wheel speeds)  
- Simulation outputs: attitude, wheel speeds, gimbal angles  
- 3D visualization of spacecraft orientation and actuator motion  

---
<img width="550" height="422" alt="f1" src="https://github.com/user-attachments/assets/c6e87177-398a-4646-af79-068dbcdf422f" />

## 📂 How to Use
1. Run **`VscmgSim.m`**  
   - This generates simulation plots for:  
     - Spacecraft attitude  
     - Wheel speeds  
     - Gimbal angles  

2. Run **`Visualizations.m`**  
   - ⚠️ Must run **`VscmgSim.m`** first to generate simulation data  
   - Provides a 3D visualization of spacecraft motion and VSCMG actuator dynamics  

---

## 📖 Reference
For more background on VSCMG theory and modeling:  
[Hanspeter Schaub, *Variable Speed Control Moment Gyroscopes*](https://hanspeterschaub.info/PapersPrivate/vscmg.pdf)  

---

## 🔑 Keywords
VSCMG, Spacecraft Dynamics, Nonlinear Control, MATLAB Simulation, Attitude Control, 3D Visualization
