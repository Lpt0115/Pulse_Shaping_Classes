 
# importing relevant pyhton modules
import h5py
import numpy as np
import scipy as sc
import sympy as sp
import matplotlib.pyplot as plt
from matplotlib.pyplot import colormaps
import matplotlib.animation as animation
from matplotlib.patches import Circle
import h5py
from scipy.optimize import curve_fit
from scipy.integrate import solve_ivp,odeint
from scipy.fft import fft,fftfreq,rfftfreq
from scipy.interpolate import interp1d
from scipy import signal
from IPython.display import HTML, display, Math

















class IQ_traj:
    """
    Workflow (in notebook file)
    --------
    >>> from ANAlysis_class_ro import RO_class_5
    >>> from ANAlysis_class_ro.RO_class_5 import IQ_traj
    >>> IQ_traj(x_data,i_data,q_data,n_cal)
    ...
    #For documentation use help()
    >>> help(IQ_traj)
    Parameters       
    ----------
    x: np.ndarray
        time data (ns)

    i: np.ndarray
        I-data

    q: np.ndarray
        Q-data

    n_cal: float
        photon number calibration factor 
 
    Attributes
    ----------
    a_out: np.ndarray, dtype=complex 
        normalized complex output field
        
    a_in: np.ndarraydtype=complex     
        normalized complex input field
        
    alpha: np.ndarray, dtype=complex   
        complex intra-cavity field
        
    pars: np.ndarray, dtype=float  
        (amplitde_drive,duration,time delay,phase_drive,_,kappa,detuning_resonator,_,_,_,phase_reference)
        
    ro_pars: np.ndarray, dtype = float 
        fitted pars   
        (amplitde_drive,duration,time delay,phase_drive,_,kappa,detuning_resonator,_,_,_,phase_reference)
        
    ro_dict: (parameter name and fitted paramtere value 
    

    Methods       
    -------  
    fid_out(self)
        Fitting the output field to the intra-cavity field model
        
        Returns
        -------
        popt: np.ndarray()
            optimized parameters
            
        pcov: (np.ndarray, np.ndarray)
            covariance matrix of the optimized parameters



    plot(self)
        1. Displays the drive amplitude, the decay rate kappa, and the detuning in MHz
         2. Plots the real and imaginary normalized measured readout data and the fitted output field to 
            the intra-cavity field model over time
         3. Parametric plot of the noramlized data and model
        
        Returns
        -------
        popt, pcov: np.ndarray(), (np.ndarray, np.ndarray)
                    fit result parameters and the covariance matrix




    plot_t_domain(self)
        1. Plots the normalized drive field over time
        2. Plots the real and imaginary part of the intra-cavity field over time  
        3. Plots the real and imaginary normalized measured readout data and the fitted output field to 
            the intra-cavity field model over time  
    """
   
            
  
            
    
    def phase_saving(self,y):
        """
        Ensures that phase is not NaN
        Parameters
        ---------
        y: np.ndarray, dtype =complex
            readout complex signal array 
            
        Returns 
        -------
        angle: float
            phase estimate that is not NaN
        """
        ph_default = 0.0
        angle = np.mean(np.angle(y))
        if np.isnan(angle):
            angle = ph_default
        else:
            pass
        return angle
        
    def d_est(self,x,y):
        """
        Estimates the detuning of the complex readout signal
        Paramters
        ---------
        x       : np.ndarray, dtype=float
                  measurment time array
        y       : np.ndarray, dtype =complex
                  complex signal array              
        Returns
        -------
        det_f: float
           estimated detuning of the signals FFT
        """
        t = x
        y_data = y
        start2 = np.argmax(y_data)
        y = y_data-np.mean(y_data)   # remove offset
        ring_dat = y[start2:]
        spec = np.abs(np.fft.rfft(np.real(ring_dat)))
        freqs = np.fft.rfftfreq(len(ring_dat),d = np.mean(np.diff(t[start2:])))
        dc_cut = min(100,max(1,len(spec)//10))
        det_f = freqs[dc_cut+np.argmax(np.abs(spec[dc_cut:]))]
        return det_f
            
    # defining initial ro parameters
    def initial_pars(self,x_data,y):
        """
        Estimation of readout paramters, extracted from measured readout signal
    
        Parameters
        ----------
        x_data: np.ndarray, dtype=float
            readout measurment time array
    
    
        y: np.ndarray, dtype = complex
            readout complex signal     
                               
        Returns 
        --------        
        pars: np.ndarray, dtype=float 
            estimated readout parameters
            [amp,dur,shift,ph,BW,kappa,deltar,K,chi,chi,ph_ro]
        """

        x_data=x_data[~np.isnan(y)]
        t = x_data
        y_data =np.asarray(y)
        y_data = y_data[~np.isnan(y)]
        y_max = np.nanmax(np.abs(y))
        if y_max == 0.0:
            y_max =1.0
        y_min = np.nanmin(y/y_max)
        amp = np.abs(y_min)
        
        
        # estimatng kappa via ring up and ring down
        start = np.nanargmin(y)
        maxamp = np.nanargmax(np.abs(y))
        if(maxamp<start):                          
            shift = t[int(maxamp)]
            dur = -(t[maxamp]-t[int(start)])
        else:
            shift = t[int(start)]
            dur = (t[maxamp]-t[int(start)])
        
        start2 = maxamp
        idx2 = np.nanargmin(np.abs(((y[start2:]/y_max)-(1./np.e))))
        t_down_s = t[start2]
        t_down_f = t[start2+idx2]
        
        t_dec = (t_down_f-t_down_s)
        kappa_est_down = 2./(t_dec)
        def kappa_down_est(x,y):
            y = np.asarray(y)
            y = y[~np.isnan(y)]
            x = x[~np.isnan(y)]
            y_max = np.nanmax(np.abs(y))
            y = y/y_max
            start =np.argmax(y)
            min_window = 50
            end = max(int(start*1.4),start+min_window)
            end = min(end,len(x))
            x = x[start:end]
            y= y[start:end]
         
            def at2(t,kappa,a):
                return a*np.exp(-2.*(t*kappa/2.))
            popt,pcov = curve_fit(at2,x,y,p0=[20.,1.0],maxfev=20000)
            return popt[0]
            
            
        def kappa_up_est(x,y):
            y = np.asarray(y)
            y = y[~np.isnan(y)]
            x = x[~np.isnan(y)]
            y=np.real(y)
            y_max = np.nanmax(np.abs(y))
            y = y/y_max
            y = np.real(y)
            start =np.argmin(y)
            end = int(np.argmax(y)*0.9)
            min_window = 50
            if (end<=start):
                end = start+min_window
            end = min(end,len(x))
            x = x[start:end]
            y= y[start:end]
            def at(t,kappa,a):
                return a*(1.0-np.exp(-(t*(kappa/2.))))-.1
            popt, pcov = curve_fit(at,x,y,p0=[9,1.0],maxfev=20000)
            return popt[0]
            
            
        ku = kappa_up_est(t,y_data)
        if np.isnan(ku):
            ku = 2*np.pi*1.
  
        
   
        kd = kappa_down_est(t,np.abs(y_data))
        kappa =np.mean([kd,ku])
        
   
      
        #estimation of detuning via fft
        d_est = self.d_est(t,y_data)/len(t)
        deltar =2*np.pi*d_est
     
        if(kappa<(0.8*2.*np.pi)):
            deltar = deltar
            
        else:
           deltar = kappa/4.
        
        
        # amplitude estimaiton
        amp = amp*kappa
        
        #phases of the drive field and the output field
        ph= ((np.pi/2)-self.phase_saving(y_data[int(maxamp*0.8):int(maxamp)]))
        ph_ro = (np.pi-self.phase_saving(y_data[:int(start*0.8)]))
        if (np.isnan(ph).any() ):
            ph = (np.pi/2.)
        else:
            pass
        if(np.isnan(ph_ro).any()):
            ph_ro = np.pi
        else:
            pass
        ph = np.angle(np.exp(1j*self.phase_saving(y_data[int(maxamp*0.8):int(maxamp)])))
        ph_ro= np.angle(np.exp(1j*self.phase_saving(y_data[:int(start*0.8)])))
        BW=2.0 # upper cutoff frequency
           
        pars = np.array([amp,dur,shift,ph,BW,kappa,deltar,ph_ro],dtype='float')
        return pars
        
        

    # Functions to calculate the model for the output field




    def square_pulse(self,t,pars):
        """
        Generates a square envelope pulse
        Parameters
        ---------
        t: np.ndarray
            time array
    
        pars: np.ndarray
            amp,dur, shift, ph,BW
                       
        Returns
        -------
        np.ndarray, dtype = complex
            complex square envelope pulse
        """
        t=t
        arr = pars
        amp,dur, shift, ph,BW = arr[:5]
        section = (np.heaviside((t-shift),1)-np.heaviside((t-(dur+shift)),1)) 
        dre= (amp*np.cos(ph))
        dimag= (amp*np.sin(ph))
        drive = (dre+1j*dimag)*section
        return drive
    
   
    def diffae(self,t,b,pars):             
        """
        Inta-cavity field differential equation Blais et al. Eq. (101)
        
        Paramters
        ---------
        t: np.ndarray, dtype = float
            time array
        
        b: [np.ndarray,np.ndarray]
            real part and imaginary part of the funciton
            
        pars: np.ndarray
            kappa, deltar, K, chi,chi2
        
        Returns
        -------
        [da_dt.real,da_dt.imag]: [np.ndarray,np.ndarray]
              real part and imaginary part of the differential equation
        """        
        pars_d = pars[:5]
        kappa, deltar = pars[5:7]
        drive = self.square_pulse(t,pars)
        a = b[0]+1j*b[1]
        da_dt = (drive.real+1j*drive.imag)-(1.0j*(deltar)+(kappa/2))*a 
        return [da_dt.real,da_dt.imag]

   
    def cav_sol(self,dt,t):     # solving of the ode with solve_ivp() method = LSODA
        """
        Solves the Drive frame intra-cavity field differential equation"
                
        Paramters
        ---------
        dt: np.ndarray, dtype = float
            kappa, deltar, K, chi,chi2
               
        t: np.ndarray, dtype = float
            time array

        Returns
        -------
        cav_g.t, alpha_g, alpha_e: np.ndarray, dtype=float,
            time array
            
        alpha_g: np.ndarray,dtype =complex
            complex Intra-cavity field solution for |g>    
            
        alpha_e: np.ndarray,dtype =complex
            complex Intra-cavity field solution for |e> 
        """
        ts = t[0]
        te=t[-1]
        neval = len(t)
        t_span = (ts,te)
        t_eval = np.linspace(*t_span,int(neval))
        t_eval = t
        mstep= np.mean(np.diff(t_eval))
        if mstep <=0:
            mstep=1e-3
           
        a0 = [0.0,0.0]
    
        cav_g = solve_ivp(self.diffae,t_span,a0,t_eval =t_eval,args=(dt,),method="LSODA",rtol=1e-8,atol=1e-10, max_step=mstep)
        alpha = cav_g.y[0]+1j*cav_g.y[1]
        return alpha
       


 
    
    def model(self,x,*p):
        """
        Calculates the complex otuput field
                
        Paramters
        ---------         
        x: [np.ndarray,np.ndarray], dtype = float
            staced time array for complexing fitting of curve_fit,
            cannot deal with complex function

        dt: np.ndarray, dtype = float
            kappa, deltar, K, chi,chi2
            
        Returns
        -------
        (np.real(a_out_nd),np.imag(a_out_nd): [np.ndarray,np.ndarray],dtype =float
            real and imaginary part of the output field staced
            for fitting to complex signal with curve_fit
        """
        t = x[:len(x)//2]
        alpha = self.cav_sol(p,t)
        pulse = self.square_pulse(t,p)
        a_in = pulse/np.sqrt(p[5])

          
        a_out = ((-a_in*np.exp(1j*p[-1])) +(alpha*np.sqrt(p[5])))
      # filtering of the sharp jumps of the output field model
        if (p[4]>0):
            sig =a_out
            sos = signal.butter(1, p[4], 'lp', fs=1.0/np.mean(np.diff(t)), output='sos')
            a_out = signal.sosfilt(sos, sig)
        a_out_nd = a_out/np.sqrt(p[5])
        return np.concatenate((np.real(a_out_nd),np.imag(a_out_nd)))

    def fid_out(self):
        """
        Fitting the output field to the intra-cavity field model
        
        Returns
        -------
        popt: np.ndarray()
            optimized parameters
            
        pcov: (np.ndarray, np.ndarray)
            covariance matrix of the optimized parameters
        """
        self.t_staced=np.concatenate((self.x,self.x))
        self.y_staced=np.concatenate(((self.i_data/self.y_max),(self.q_data/self.y_max)))
        
        
        
        # preparing and running of opt with curve fit
        
        # defining bounds for the guessed parameters
        dt = self.pars
        p0 = self.pars
        
        
     
            
        tol =.3
        min_tol = 1e-10
        lb =np.abs(np.array(p0))*(1-tol)
        ub =np.abs(np.array(p0))*(1+tol)+min_tol
        
        
        pars_names = [
            'amplitude_drive',
            'duration',
            'shift',
            'ph_drive',
            'BW',
            'kappa',
            'detuning',
            'ph_ref',
        ]
        
        pars_bounds =   {
        
            'amplitude_drive':
                {   'lb': p0[0]*0.5,
                    'ub': p0[0]*10.,  
                }
            ,
            'duration': 
               {  
                    'lb': 0.999*p0[1],
                    'ub':   1.001*ub[1],   
                }
            ,
            'shift': 
               {   
                    'lb': p0[2]*0.999,
                    'ub':    p0[2]*1.001,   
                }
            ,
            'ph_drive': 
               {  
                    'lb': -np.pi,
                    'ub': np.pi, 
                }
            ,
            'BW': 
               {   
                    'lb': 0.0,
                    'ub':   10.0  
                }
            ,
            'kappa':
               {   
                    'lb': 0.01*p0[5],
                    'ub': 4.*p0[5] ,   
                }
            ,
            'detuning':
               {   
                    'lb': -4.*p0[5],
                    'ub':    4.*p0[5],   
                }
            ,
            'ph_ref':
               {  
                    'lb': -np.pi,
                    'ub':  np.pi*.9,
                }
        
        }
        
        for i, name in enumerate(pars_names):
            if name in pars_bounds:
                if 'lb' in pars_bounds[name]:
                    lb[i] = float(pars_bounds[name]['lb'])
                if 'ub' in pars_bounds[name]:
                    ub[i] = float(pars_bounds[name]['ub'])
        
        
        
        for i in range (len(p0)):
            if((lb[i]>p0[i]) or (p0[i]>ub[i])):
                print(f"Traceback: Intital pars {pars_names[i][:]} : {p0[i]} not within boundaries'")
            elif(i==len(p0)):
                print("Initial pars within boundaries:...  \n")
        
        
        
        #Fitting to the model via curve_fit()
        popt, pcov = curve_fit(self.model, self.t_staced,self.y_staced.reshape(-1),p0=p0,bounds=(lb,ub),method='dogbox',check_finite=True )
        self.popt,self.pcov = popt, pcov
        self.ro_pars = self.popt
        self.ro_dict = dict(zip(self.ro_dict.keys(),self.ro_pars))
        return popt, pcov
    def plot(self):
        """
        1. Displays the drive amplitude, the decay rate kappa, and the detuning in MHz
        2. Plots the real and imaginary normalized measured readout data and the fitted output field to 
            the intra-cavity field model over time
        3. Parametric plot of the noramlized data and model
        """
 
        popt, pcov = self.popt, self.pcov
        t_staced = self.t_staced
        y_staced = self.y_staced
        display(Math(fr"$\epsilon / 2\pi = \left( {popt[0]/(2*np.pi):.3f} \pm {np.sqrt(np.diag(pcov))[0]/(2*np.pi):.3f}\right) \, \rm{{MHz}}$"))
        for i in range(2):
            name_vals = np.array([r"\kappa / 2\pi", r"\Delta / 2\pi"])
            display(Math(fr"{name_vals[i]} = \left( {popt[5+i]/(2*np.pi):.3f} \pm {np.sqrt(np.diag(pcov))[5+i]/(2*np.pi):.3f}\right) \, \rm{{MHz}}"))




        # Plotting of fit
        r = popt[0]/popt[5]
        circ_a_out = Circle((0.0,0.0),r,edgecolor=(24/255,126/255,73/255),linestyle='--',facecolor='none',label='steady state') 
        fig, ax = plt.subplots(2,1,figsize=(6,5))
        

        ax[0].plot(self.x,self.model(t_staced,*popt)[:len(self.x)],'--',label='model real fit')
        ax[0].plot(self.x,y_staced[:len(self.x)],label='data real',alpha=0.5)
        ax[1].plot(self.x,self.model(t_staced,*popt)[len(self.x):],'--',label='model imag fit')
        ax[1].plot(self.x,y_staced[len(self.x):],label='data imag',alpha=0.5)
       
        for i in range(2):
            ax[i].grid(True,'major')
            ax[i].legend()
            if (i==1):
                ax[i].set_xlabel("t (µs)")
        ax[0].set_ylabel(r"$\rm{Re}\left(\rm{a}_{\rm{out}}\right) / \sqrt{\kappa}$")
        ax[1].set_ylabel(r"$\rm{Im}\left(\rm{a}_{\rm{out}}\right) / \sqrt{\kappa}$")
        
        plt.tight_layout()
        fig,ax = plt.subplots(1,1)
        ax.plot(self.model(t_staced,*popt)[:len(self.x)],self.model(t_staced,*popt)[len(self.x):],'--',color='brown',label='model fit')
        bar = ax.scatter(y_staced[:len(self.x)],y_staced[len(self.x):],c=self.x,cmap='Purples')
        fig.colorbar(bar,label="t (µs)")
        ax.add_patch(circ_a_out)
        ax.set_aspect(1)
        ax.set_xlim()
        ax.set_ylim()
        ax.set_xlabel(r"$\rm{Re} \left(\rm{a}_{\rm{out}}\right)/ \sqrt{\kappa}$")
        ax.set_ylabel(r"$\rm{Im} \left(\rm{a}_{\rm{out}}\right) / \sqrt{\kappa}$")
        ax.legend()
        ax.grid(True,'major')            


    def plot_t_domain(self):
        """
        1. Plots the normalized drive field over time
        2. Plots the real and imaginary part of the intra-cavity field over time  
        3. Plots the real and imaginary normalized measured readout data and the fitted output field to 
            the intra-cavity field model over time
        """
        r = np.abs(self.ro_pars[0])/self.ro_pars[5]
        fi, ax = plt.subplots(4,1,sharex=True)
        ax[2].plot(self.x,self.y.real)
        ax[3].plot(self.x,self.y.imag)
        ax[1].plot(self.x,self.alpha.real,color='blue')
        ax[1].plot(self.x,self.alpha.imag,color='orange')
        ax[2].plot(self.x,self.a_out.real)
        ax[3].plot(self.x,self.a_out.imag)
        ax[0].plot(self.x,self.ain.real/np.sqrt(self.ro_pars[5]),color='blue')
        ax[0].plot(self.x,self.ain.imag/np.sqrt(self.ro_pars[5]),color='orange')
        ax[2].plot(self.x,-self.ain.real/np.sqrt(self.ro_pars[5]),color='blue')
        ax[3].plot(self.x,-self.ain.imag/np.sqrt(self.ro_pars[5]),color='orange')
        ax[3].set_xlabel("t (µs)")
        ax[0].set_ylabel(r"$\rm{Re}\left( \rm{a}_{\rm{in}} \right) / \sqrt{\kappa}$",color='blue')
        ax[1].set_ylabel(r"$\rm{Re}\left( \alpha \right)$",color='blue')
        ax[2].set_ylabel(r"$\rm{Re}\left( \rm{a}_{\rm{out}} \right) / \sqrt{\kappa}$")
        ax[3].set_ylabel(r"$\rm{Im} \left( \rm{a}_{\rm{out}} \right) / \sqrt{\kappa}$")
        ax[0].twinx().set_ylabel(r"$\mathrm{Im} \left( \mathrm{a}_{\mathrm{in}} \right) / \sqrt{\kappa}$",color='orange')
        ax[1].twinx().set_ylabel(r"$\mathrm{Im} \left( \alpha \right)$",color='orange')









   
   
   

    def __init__(self,x,i,q,n_cal):
    
        """
        Parameters       
        ----------
        x: np.ndarray
            time data (ns)
    
        i: np.ndarray
            I-data 
    
        q: np.ndarray
            Q-data
    
        n_cal: float
            photon number calibration factor 
        """

       

        # relevant paramters 
        self.ro_dict = {
            'amplitude_drive': 0.0,
            'duration': 0.0,
            'shift': 0.0,
            'ph_drive': 0.0,
            'BW': 0.0,
            'kappa': 0.0,
            'detuning': 0.0,
            'ph_ref':0.0,
        }
        
       
        #normalizing the amplitudes and scaling the time unit
        x_data = x
        self.x = (x_data/1e-6)                # conversion from ns to µs
                                              
        i_data = i
        q_data = q
        self.i_data = i_data
        self.q_data = q_data
        y_data_abs = np.sqrt(i_data**(2)+q_data**(2))
        y_max = np.max(y_data_abs)
        
        self.y_max = y_max
        self.sqrt_n = np.sqrt(n_cal)*self.y_max # output field scaling from max Voltage to max_sqrt(n)  
        y_data = np.array(i_data+1j*q_data)
        self.y = (y_data/self.y_max)
        
        self.pars = self.initial_pars(self.x,self.y)

        
            
            # Printing the fit pars drive amp., kappa and Detunug
        self.popt,self.pcov = self.fid_out()
        self.ro_pars = self.popt
        self.ro_dict = dict(zip(self.ro_dict.keys(),self.ro_pars))

        self.square = self.square_pulse(self.x,self.ro_pars)
        
        self.ain = self.square/np.sqrt(self.ro_pars[5])
        
        self.alpha = self.cav_sol(self.ro_pars,self.x)
        
        a_out = ( -(self.ain*np.exp(1j*self.ro_pars[-1])) + np.sqrt(self.ro_pars[5])*self.alpha )
        
        self.a_out = a_out / np.sqrt(self.ro_pars[5])
        # converting amp in MHz into sqrt(n) sqrt(n) = amp/(2*np.pi)/sqrt((kappa/2)**(2)+deltar**(2)) , steady state realtion
        self.amp_sqrt_n = (self.ro_dict['amplitude_drive'])/np.sqrt((self.ro_dict['kappa']/2.)**(2)+(self.ro_dict['detuning'])**(2))
        
        # Calculating the photon number 
        n_int=np.sum(np.abs(self.alpha)**(2))*(self.x[-1]/len(self.x))*self.ro_dict['kappa']
        self.amp =np.sqrt((n_int/2200))
        
       
        #self.fid_out()
        self.plot()
        self.plot_t_domain()
