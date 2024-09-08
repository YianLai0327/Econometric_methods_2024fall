# HW1 0904: Linear Algebra

In HW1, we implement some basic linear algebra operation by R.

---

## Q1. **Data Visualization**

![Histogram of y](plots/histogram_y.png)
![Time series of y](plots/time_series_y.png)
![Histogram of x_dfy](plots/histogram_x_dfy.png)
![Time series of x_dfy](plots/time_series_x_dfy.png)
![Histogram of x_infl](plots/histogram_x_infl.png)
![Time series of x_infl](plots/time_series_x_infl.png)
![Histogram of x_svar](plots/histogram_x_svar.png)
![Time series of x_svar](plots/time_series_x_svar.png)
![Histogram of x_tms](plots/histogram_x_tms.png)
![Time series of x_tms](plots/time_series_x_tms.png)
![Histogram of x_tbl](plots/histogram_x_tbl.png)
![Time series of x_tbl](plots/time_series_x_tbl.png)
![Histogram of x_dfr](plots/histogram_x_dfr.png)
![Time series of x_dfr](plots/time_series_x_dfr.png)
![Histogram of x_dp](plots/histogram_x_dp.png)
![Time series of x_dp](plots/time_series_x_dp.png)
![Histogram of x_ltr](plots/histogram_x_ltr.png)
![Time series of x_ltr](plots/time_series_x_ltr.png)
![Histogram of x_ep](plots/histogram_x_ep.png)
![Time series of x_ep](plots/time_series_x_ep.png)
![Histogram of x_bmr](plots/histogram_x_bmr.png)
![Time series of x_bmr](plots/time_series_x_bmr.png)

---

## Q2. **Matrix Trace Calculations**

Two trace calculations were performed using R:

- **Q2.1:** Trace of $X (X'X)^{-1} X'$ is calculated as follow:
  
  $\text{trace}(X (X'X)^{-1} X') = 11$

- **Q2.2:** Trace of $I - X (X'X)^{-1} X'$ is computed as follow:
  
  $\text{trace}(I_n - X (X'X)^{-1} X') = 493$

---

## Q3. **Standardized Eigenvalues Calculation**

![Scree plots of Q3](scree_plots/scree_plot_q3.png)

---

## Q4. **Further Eigenvalue Analysis**

![Scree plots of Q4](scree_plots/scree_plot_q4.png)

---

## Q5. **Spectral Decomposition**

Using R, we apply spectral decomposition for the matrix $\widetilde{X}$, where:

$H = \widetilde{X}' \cdot \widetilde{X}$

The matrix $A$ is calculated as:

$A = H \Lambda^{-1} H'$

This result confirms the equation provided in the problem statement.

---

## Q6. **Best Fitted b Calculation**

The best fitted $b$ is calculated as the projection of $y$ in the row space of $\widetilde{X}$. The formula used is:

$b = (\widetilde{X}'\widetilde{X})^{-1} \widetilde{X}'y$

Using R, the resulting vector $b$ is computed successfully.

---
