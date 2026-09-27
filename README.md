# Introduction
These classes are written for the master's thesis project "Pulse Shaping for optimal superconducting qubit readout", Philipp Thamm at Karlsruhe Institut for Thehnology (KIT), 2026.

As the name implies, this project is embedded in the research field of quantum computing and deals with a physics motivated engineering task.
Without going further into the details for readers with less physical background, a brief motivation will be outlined:

Measuring a qubit state is not as simple as for a classical measurement, since it is a quantum state.
A way to measure the qubit state is realized by coupling the system to an quantum-mechanical LC oscillator.
The resonance fequency of the resonator becomes dependent on the qubit state.
This allows to distinguish the qubit state by sending a readout microwave signal
with a moderate amplitude at frequency between the two qubit state-dependent
resonance frequencies $\omega_{rf}$ and to the resonator
and measuring the reflected signal.


Before starting with pulse shaping, it is essential to know the sepcific parameters of the experimental setup for qubit readout.
A square pulse readout analysis class was written for parameter extraction of time domain readout signals that is as a prerequisite for pulse shaping.

Square pulse for superconducting qubit readout comes at cost of a rather slow decay of the readout signal after the information about the qubit state is obtained.
For a more optimized readout a method of adding square pulses was used [^1] [^2]. This method aims for a signal suppression after the state of the qubit has already been obtained.


# Square pulse analysis class
The class uses a linear intra-cavity field $\alpha(t)$ model to extract the linwidth of the resonator $\kappa$ and the detuning $\delta$ of the readout tone to the resonance frequnecy for a measured reflected output signal by fitting the model to the complex data. 
## Usage
----------
### Importing the class
-------------------------
```python
#Assuming ANAlysis_class_ro is the working directory
from ANAlysis_class_ro import RO_class_5
from ANAlysis_class_ro.RO_class_5 import IQ_traj
```
### Intitializing and executing the class
The folowing parameters are passed to the class for intitializing:
|Parameter|Type|Description|
|---------|----|-----------|
|'x_data'|'np.ndarray'|time data (ns)|
|'i_data'|'np.ndarray'|I-quadrature data (V)|
|'q_data'|'np.ndarray'|Q-quadrature data (V)|
|'n_cal'|'float'|AC-Stark calibration factor ($\braket{\rm n} / \rm a^{2}$)|
```python
IQ_traj(x_data,i_data,q_data,n_cal)
```

    

### Methods
#### These methods are executed internaly by the class after initialization:
-------  

```python
    plot(self)
```
1. Displays the drive amplitude, the decay rate kappa, and the detuning in MHz
2. Plots the real and imaginary normalized measured readout data and the fitted output field to 
    the intra-cavity field model over time
3. Parametric plot of the noramlized data and model



```python
plot_t_domain(self)
````
1. Plots the normalized drive field over time
2. Plots the real and imaginary part of the intra-cavity field over time  
3. Plots the real and imaginary normalized measured readout data and the fitted output field to 
    the intra-cavity field model over time
   
### Parameter storage
The extracted parameters can be accesed by the parameter dictionary attribute of the class 
```python
self.ro_dict
```

### Documentation reference
For detailed information on specific modules of the class the help() module can be used.
```python
#For documentation use help()
help(RO_class_5.IQ_traj)
```
### Mititgation of issues    
Providing experimental readout data, that has a rather smooth and pronunced square pulse response is helpful for parameter estimation. This implies a readout duration that is sufficiently long enough for reaching a steady-state response. If the fit fails, consider to apply a moving average with a window size that is not too coarse, this might help the fitting process.
## Working principle
This class uses the input output relations to analyse the readout signal, if a square pulse is used. Thereby, the relevant experimental parameters can be estimated.
The intra-cavity parameters such as the linewidth $\kappa$ and the dispersive shift $\chi$ can be experimentally obtained by the state-dependent resonator phase response. Similarly, these parameters can be found by analyzing the readout signal of a square pulse.
The class normalizes the amplitude of both quadratures by the maximum of the absolute amplitude of the signal.
With the help of a function method, makes a rough estimate of readout parameters:
<p align="center">
  <img src="images/Square_Pulse_Readout_Class_workflow.png" width="500">
</p>
Applying this analysis class to both qubit state-dependent readout signal, the difference in detuning gives the frequency shift between the two responses of the resonator which is linked to the coupling to the qubit.

[^1]: McClure et al.

[^2]: Hazra et al., 2024
# Three-Segment Optimization Class
This class estimates optimal pulse parameters for a qubit readout using a Three-segment pulse with total pulse duration T for a specific readout setup, characterized by the linewidth of the readout resonator $\kappa/2\pi$ and the dispersive shift $\chi/2\pi$ of the resonance frequency of the resonator dependent on the qubit state either in $\ket{\rm g}$ or in $\ket{\rm e}$. The drive amplitude is capped to a specific maximum mean photon number $\braket{\rm n}_{\rm max}$ induced by the Three-segment pulse.
## Usage
----------
### Importing the class
-------------------------
```python
from Pulse_Shaping_Classes.src.Three_Segment_Optimization import opt_three_segment_final
```
### Intitializing and executing the class
For initailizing the class these parameters are required:
|Parameter|Type|Description|
|---------|----|-----------|
|'T'|'float'|total pulse duration (µs)|
|'kappa'|'float'| decay rate $\kappa$ of the resonator (1e6 rad/s)|
|'chi'|'float'|dispersive shift $\chi$ (1e6 rad/s)|
|'n_max'|'float'|Maximum mean photon number $\braket{\rm n}_{\rm max}$|
|'three_seg_pulse_pars'|'np.ndarray'|Initial three-segment pulse parameters ($\sum_{i=1}^{3}{\Delta t_{i} \, a_{i}})$|
```python
three_seg(T,kappa,chi,n_max,three_seg_pulse_pars)
````

### Methods
#### Internal methods called by the class after initialization:
-------  
```python
opt_out(dt)
````
Executes the optimization for the Three-Segment Pulse by minimizing the cost function
##### Returns
##### 'res': Optimized Three-segment pulse parameters



```python
    plot_pulse()
```
Plotting of the optimized Three-Segment Pulse (normalized by sqrt(kappa)) over time 
```python
    plot(self.res)
```
Plotting of the intra-cavaity field trajectory for the qubit in $\ket{\rm{g}}$ and in $\ket{\rm{e}}$.
```python
    get_SNR_ratio(self.res,self.opt_square)
```
1. Plots the squared instantaneous state separation over time for the three-Segment Pulse compared to the instant state separation for the Square Pulse
2. plots the instantaneous mean photon number over time for the Three-segment Pulse and the SQuare pulse.
#### Additional methods the class provides:
```python
    plot_out(self.res)
```
1. Plots the output field trajectory on the complex plane.
2. Plots the output field over time
### Parameter storage
The optimized thrree-segment pulse paramteres are accessible through the class attribute
```python
    self.res
```
### Mititgation of issues
If the optimization fails or if the "otpimized" paramters are not as desired, the intitial three-segment pulse paramters can be adapted or the total pulse duration can be adjusted.
Optimization will benefit, if experimental setup parameters such as the linewidth of the resonator $\kappa/2\pi$ and the dispersive shift $\chi/2\pi$ are extracted from measurements as precise as possible
### Documentation reference
```python
help(opt_three_segment_final.three_seg)
```

## Working principle
The optimization class uses the readout duration $T$, the resonator linewidth $\kappa$, the dispersive shift $\chi$, the maximum mean photon number $\rm n_{\mathrm{max}}$ and the  initial Three-segment pulse parameters ($\Delta t_{i} \, a_{i}$) for initialization and execution of the class.

Inside the class, from the initial only the free pulse parameters are passed to the python module <bf>scipy.optimize.minimize()</bf>

Using the integrated sate separation ratio as a cost function, the optimal pulse parameters are returned:

<p align="center">
  <img src="images/Three_Segment_Class_workflow.png" width="500">
</p>


