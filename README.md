# Three-Segment Optimization Class
The optimization class uses the readout duration $T$, the resonator linewidth $\kappa$, the dispersive shift $\chi$, the maximum mean photon number $\rm n_{\mathrm{max}}$ and the  initial Three-segment pulse parameters ($\Delta t_{i} \, a_{i}$) for initialization and execution of the class.
Inside the class, from the initial only the free pulse parameters are passed to the python module \mintinline{python}{scipy.optimize.minimize()}.
Using the integrated sate separation ratio \autoref{eq:cost_three_seg} as a cost function, the optimal pulse parameters are returned \autoref{flow-chart:Three_seg_class}.
The robustness of the optimized pulse parameters to the values of $\kappa$  and $\chi$ on real setups is estimated by varying the values for $\kappa$ and $\chi$ while the optimized readout parameters are fixed and the instant mean photon number $n$ after the end of the Three-segment pulse is compared to the residual photon number $\rm n_0$ because of thermal occupation \autoref{tab:opt_three_sens_n_0} and more important to the instant mean photon number for a square pulse $\rm n /  \rm n_{sq} >1$ \autoref{tab:opt_three_sens_n_sq}. 


![](images/Three_Segment_Class-workflow.png)
