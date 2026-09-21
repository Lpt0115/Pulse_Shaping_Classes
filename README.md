
# Working principles of the classes



## Square pulse analysis class
This class uses the input output relations to analyse the readout signal, if a square pulse is used. Thereby, the relevant experimental parameters can be estimated.
The intra-cavity parameters such as the linewidth $\kappa$ and the dispersive shift $\chi$ can be experimentally obtained by the state-dependent resonator phase response. Similarly, these parameters can be found by analyzing the readout signal of a square pulse.
Parameter extraction is performed using a square pulse analysis class. 
The square pulse analysis class  is initialized by passing the time array, the I data, Q data and the photon calibration from an ac stark measurement:

IQ_traj(x_data,i_data,q_data,n_cal) 

The class normalizes the amplitude of both quadratures by the maximum of the absolute amplitude of the signal.
With the help of a function method, makes a rough estimate of readout parameters:



<p align="center">
  <img src="images/Square_Pulse_Readout_Class_workflow.png" width="300">
</p>

## Three-Segment Optimization Class
The optimization class uses the readout duration $T$, the resonator linewidth $\kappa$, the dispersive shift $\chi$, the maximum mean photon number $\rm n_{\mathrm{max}}$ and the  initial Three-segment pulse parameters ($\Delta t_{i} \, a_{i}$) for initialization and execution of the class.

Inside the class, from the initial only the free pulse parameters are passed to the python module <bf>scipy.optimize.minimize()</bf>

Using the integrated sate separation ratio as a cost function, the optimal pulse parameters are returned:

<p align="center">
  <img src="images/Three_Segment_Class_workflow.png" width="300">
</p>


