# **A Proximal Gradient Method with an Explicit Linesearch for Multiobjective Optimization**

This repository contains the source codes for the numerical experiments described in the article:

> **Y. Bello-Cruz, J.G. Melo, L.F. Prudente, and R.V.G. Serra,**  
> *A Proximal Gradient Method with an Explicit Linesearch for Multiobjective Optimization*, technical report, 2024.

**Date:** November 2024

---

## **Contents**

The repository includes the following routines:

### **Main Routine**
- **`main.m`**  
  Controls the execution of the experiments, allowing users to:
  - Select the problem to be solved.
  - Choose the optimization algorithm to apply.

### **Algorithms**
- **`ProxGrad.m`**  
  Contains implementations of:
  - **Explicit-MPG**: Explicit Multiobjective Proximal Gradient algorithm.
  - **Armijo-MPG**: Multiobjective Proximal Gradient with Armijo line search.
  - **Implicit-MPG**: Implicit Multiobjective Proximal Gradient algorithm.
  
- **`ProxGradAcc.m`**  
  Contains implementations of:
  - **Accelerated-MPG**: Accelerated Multiobjective Proximal Gradient algorithm.
  - **Normal-MPG**: Normal Proximal Gradient algorithm.

### **Auxiliary Routines**
- **`inip.m`**  
  Defines the initial parameters for the optimization problems.
  
- **`evalg.m`**  
  Implements the functions G_j.

- **`evalgradg.m`**  
  Computes the gradient of the functions G_j.

- **`evalh.m`**  
  Implements the functions H_j.

- **`datas.m`**  
  Generates the data for the functions H_j.

- **`subproblemProxGrad.m`**  
  Solves the subproblem required in the proximal gradient methods.

### **Line Search Routines**
- **`explicitLS.m`**  
  Implements the explicit line search used in Explicit-MPG.
  
- **`explicitLSone.m`**  
  Handles the line search for a single objective function in Explicit-MPG.
- **`armijo.m`**  
  Implements the Armijo line search used in Armijo-MPG.

- **`implicitLS.m`**  
  Implements the implicit line search used in Implicit-MPG.

---

## **Instructions**

To run the numerical experiments:
1. Open the file `main.m`.
2. Select the desired problem and algorithm in the code.
3. Run the script and view the results on the screen.

---

## **License**

This project is licensed under the **GNU General Public License (GPL) version 3**. For more information, see the LICENSE file.
