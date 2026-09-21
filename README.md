\section{Three-Segment Optimization Class}\label{app:three_seg_class}
The optimization class uses the readout duration $T$, the resonator linewidth $\kappa$, the dispersive shift $\chi$, the maximum mean photon number $\mathrm{n}_{\mathrm{max}}$ and the  initial Three-segment pulse parameters $(\Delta t_{i},\, a_{i})$ for initialization and execution of the class.
Inside the class, from the initial only the free pulse parameters are passed to the python module \mintinline{python}{scipy.optimize.minimize()}.
Using the integrated sate separation ratio \autoref{eq:cost_three_seg} as a cost function, the optimal pulse parameters are returned \autoref{flow-chart:Three_seg_class}.
The robustness of the optimized pulse parameters to the values of $\kappa$  and $\chi$ on real setups is estimated by varying the values for $\kappa$ and $\chi$ while the optimized readout parameters are fixed and the instant mean photon number $\mathrm{n}$ after the end of the Three-segment pulse is compared to the residual photon number $\mathrm{n}_{0}$ because of thermal occupation \autoref{tab:opt_three_sens_n_0} and more important to the instant mean photon number for a square pulse $\mathrm{n}/\mathrm{n}_{\mathrm{sq}}>1$ \autoref{tab:opt_three_sens_n_sq}. 
\begin{figure}[h]
\centering
\begin{tikzpicture}
\tikzset{box/.style={draw, rounded corners = 3pt, rectangle, minimum width= .8cm, minimum height=1cm,align=center, fill=brown},
arrow/.style= {-stealth,thick}
}
    \node[box] (raw) {$T,\, \kappa\,, \chi,\, \mathrm{n}_{\mathrm{max}}$, \\ $\text{Pulse pars}(\Delta t_{i},\, a_{i})$};
    \node[box, below=2cm] (init) {
    $\text{Three-segment}(T, \kappa, \chi, \mathrm{n}_{\mathrm{max}})$};
    \draw[arrow] (raw) -- node[midway, right]{Initialization} (init);
    \node[box, below =4cm] (model) {Cost function: $
    \mathrm{min}_{\Delta t_{i},\ a_{i}}\left( \dfrac{\int_{T}^{\infty} dt\, | \alpha_{\ket{g}} - \alpha_{\ket{e}}|^{2}}{\int_{0}^{T} dt \, | \alpha_{\ket{g}} - \alpha_{\ket{e}} |^{2}}. \right)$};
    \draw[arrow] (init) -- node[midway, right]{Computing} (model);
        \node[box, below =6cm] (res) {Returns optimized pulse parameters: \\ segment durations $\Delta t_{i}$ and segment amplitudes $a_{i}$};
    \draw[arrow] (model) -- node[midway, right]{Minimizing} (res);
\end{tikzpicture}
  \caption[Workflow behind the Three-segment pulse optimization class]{\textbf{Workflow behind the Three-segment pulse optimization class} Initializing the class with the setup parameters runs the optimization using the cost function \autoref{eq:cost_three_seg} and returns the optimized pulse parameters $\Delta t_{i}$ and $a_{i})$ with $i=1,\,2,\,3$.}
    \label{flow-chart:Three_seg_class}
\end{figure} 
