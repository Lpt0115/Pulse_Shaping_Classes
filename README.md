# Introduction
These classes are written for the master's thesis project "Pulse Shaping for Optimal Superconducting Qubit Readout", Philipp Thamm at Karlsruhe Institute of Technology (KIT), 2026.

As the name implies, this project is embedded in the research field of quantum computing and deals with a physics motivated engineering task.
Without going further into the details for readers with less physical background, a brief motivation will be outlined:

Measuring a qubit state is not as simple as for a classical measurement, since it is a quantum state.
A way to measure the qubit state is realized by coupling the system to an quantum-mechanical LC oscillator.
The resonance fequency of the resonator becomes dependent on the qubit state.
This allows to distinguish the qubit state by sending a readout microwave signal
with a moderate amplitude at frequency $\omega_{rf}$ between the two qubit state-dependent
resonance frequencies to the resonator
and measuring the reflected signal.


Before starting with pulse shaping, it is essential to know the sepcific parameters of the experimental setup for qubit readout.
A square pulse readout analysis class was written for parameter extraction of time domain readout signals that is as a prerequisite for pulse shaping.

Square pulse for superconducting qubit readout comes at cost of a rather slow decay of the readout signal after the information about the qubit state is obtained.
For a more optimized readout a method of adding square pulses was used [^1] [^2]. This method aims for a signal suppression after the state of the qubit has already been obtained.


# Square pulse analysis class
The class uses a linear intra-cavity field $\alpha(t)$ model to extract the linewidth of the resonator $\kappa/2\pi$ and the detuning $\delta/2\pi$ of the readout tone to the resonance frequency for a measured reflected output signal by fitting the model to the complex data. 
## Usage
----------
### Importing the class
-------------------------
```python
from Pulse_Shaping_Classes.src.Square_Pulse_Readout.square_readout_analysis import IQ_traj
```
### Initializing and executing the class
The following parameters are passed to the class for initializing:
|Parameter|Type|Description|
|---------|----|-----------|
|'x_data'|'np.ndarray'|time data (ns)|
|'i_data'|'np.ndarray'|I-quadrature data (V)|
|'q_data'|'np.ndarray'|Q-quadrature data (V)|
|'n_cal'|'float'|AC-Stark calibration factor ($\braket{\rm n} / \rm a^{2}$)|
```python
ro_square = IQ_traj(x_data,i_data,q_data,n_cal)
```

    

### Methods
#### These methods are executed internally by the class after initialization:
-------  

```python
   ro_square.plot()
```
1. Displays the drive amplitude $\epsilon /2\pi$, the linewidth of the resonator $\kappa / 2\pi$, and the detuning $\delta / 2\pi$ in MHz
2. Plots the real and imaginary normalized measured readout data and the fitted output field to 
    the intra-cavity field model over time
3. Parametric plot of the normalized data and model



```python
ro_square.plot_t_domain()
```
1. Plots the normalized drive field over time
2. Plots the real and imaginary part of the intra-cavity field over time  
3. Plots the real and imaginary normalized measured readout data and the fitted output field to 
    the intra-cavity field model over time
   
### Parameter storage
The extracted parameters can be accessed by the parameter dictionary attribute of the class 
```python
ro_square.ro_dict
```

### Documentation reference
For detailed information on specific modules of the class the help() module can be used.
```python
help(ro_square)
```
### Mitigation of issues    
Providing experimental readout data, that has a rather smooth and pronounced square pulse response is helpful for parameter estimation. This implies a readout duration that is sufficiently long for reaching a steady-state response. If the fit fails, consider to apply a moving average with a window size that is not too coarse, this might help the fitting process.
## Working principle
This class uses the input output relations to analyse the readout signal, if a square pulse is used. Thereby, the relevant experimental parameters can be estimated.
The intra-cavity parameters such as the linewidth $\kappa/2\pi$ and the dispersive shift $\chi$ can be experimentally obtained by the state-dependent resonator phase response. Similarly, these parameters can be found by analyzing the readout signal of a square pulse.
The class normalizes the amplitude of both quadratures by the maximum of the absolute amplitude of the signal.
With the help of a function method, makes a rough estimate of readout parameters:
<p align="center">
  <img src="images/Square_Pulse_Readout_Class_workflow.png" width="500">
</p>
Applying this analysis class to both qubit state-dependent readout signal, the difference in detuning gives the frequency shift between the two responses of the resonator $2 \chi / 2\pi$ which is linked to the coupling to the qubit.

[^1]: McClure et al., Phys. Rev. Applied 5,2016.
[^2]: Hazra et al., Phys. Rev. Lett. 134, 2025.
# Three-Segment Optimization Class
This class estimates optimal pulse parameters for a qubit readout using a three-segment pulse with total pulse duration T for a specific readout setup, characterized by the linewidth of the readout resonator $\kappa/2\pi$ and the dispersive shift $\chi/2\pi$ of the resonance frequency of the resonator dependent on the qubit state either in $\ket{\rm g}$ or in $\ket{\rm e}$. The drive amplitude is capped to a specific maximum mean photon number $\braket{\rm n}_{\rm max}$ induced by the three-segment pulse.
## Usage
----------
### Importing the class
-------------------------
```python
from Pulse_Shaping_Classes.src.Three_Segment_Optimization.opt_three_segment_final import three_seg
```
### Initializing and executing the class
For initializing the class these parameters are required:
|Parameter|Type|Description|
|---------|----|-----------|
|'T'|'float'|total pulse duration (µs)|
|'kappa'|'float'| decay rate $\kappa$ of the resonator ($10^{6}$ rad/s)|
|'chi'|'float'|dispersive shift $\chi$ ($10^{6}$ rad/s)|
|'n_max'|'float'|Maximum mean photon number $\braket{\rm n}_{\rm max}$|
|'three_seg_pulse_pars'|'np.ndarray'|Initial three-segment pulse parameters ($\Delta t_{i}$, $a_{i}$)|

The durations of the segments $\Delta t_{i}$ are passed in µs and the amplitudes $a_{i}$ of the segments are relative to the drive amplitude $\epsilon$, which leads to the maximum mean photon number $\rm n_{\rm max}$ number inside the resonator $\epsilon = \sqrt{\rm n_{\rm max}}\sqrt{(\kappa/2)^{2}+(\chi)^{2}}$

```python
three_opt = three_seg(T,kappa,chi,n_max,three_seg_pulse_pars)
````
### Methods
#### Internal methods called by the class after initialization:
-------  
```python
three_opt.opt_out(three_seg_pulse_pars)
````
Executes the optimization for the three-segment pulse by minimizing the cost function
##### Returns
##### 'res': Optimized three-segment pulse parameters
```python
three_opt.plot_pulse()
```
Plotting of the optimized three-segment pulse (normalized by $\sqrt{\kappa}}) over time 
```python
three_opt.plot(three_opt.res)
```
Plotting of the intra-cavity field trajectory for the qubit in $\ket{\rm{g}}$ and in $\ket{\rm{e}}$.
```python
three_opt.get_SNR_ratio(three_opt.res,three_opt.opt_square)
```
Plots the squared instantaneous state separation over time for the three-segment pulse compared to the instantaneous state separation for the square pulse
#### Additional methods the class provides:

```python
three_opt.get_n_t(three_opt.res,three_opt.opt_square)
```
Plots the instantaneous mean photon number over time for the three-segment pulse and the square pulse.
```python
three_opt.plot_out(three_opt.res)
```
1. Plots the output field trajectory on the complex plane.
2. Plots the output field over time
### Parameter storage
The optimized three-segment pulse paramteres are accessible through the class attribute
```python
   three_opt.res
```
### Mititgation of issues
If the optimization fails or if the "optimized" paramters are not as desired, the initial three-segment pulse parameters can be adapted or the total pulse duration can be adjusted.
Optimization will benefit, if experimental setup parameters such as the linewidth of the resonator $\kappa/2\pi$ and the dispersive shift $\chi/2\pi$ are extracted from measurements as precise as possible
### Documentation reference
```python
help(three_opt)
```

## Working principle
The optimization class uses the readout duration $T$, the resonator linewidth $\kappa$, the dispersive shift $\chi$, the maximum mean photon number $\rm n_{\mathrm{max}}$ and the  initial 
three-segment pulse parameters ($\Delta t_{i} \, a_{i}$) for initialization and execution of the class.

Inside the class, from the initial only the free pulse parameters are passed to the python module <bf>scipy.optimize.minimize()</bf>

Using the integrated sate separation ratio as a cost function, the optimal pulse parameters are returned:

<p align="center">
  <img src="images/Three_Segment_Class_workflow.png" width="500">
</p>

The four-segment class uses the same working principle and is used analogous to the three-segment class.
