import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import scipy as sc
import sympy as sp
from sympy import diff,symbols,sin,pi,I,cos
from sympy.utilities.lambdify import lambdastr
from sympy.printing.numpy import NumPyPrinter
from scipy.integrate import odeint,solve_ivp
from numpy import heaviside
from scipy import optimize
from scipy.optimize import minimize, NonlinearConstraint








class four_seg:
    """
    Workflow (in notebook file)
    --------
    >>>import opt_four_segment_final
    >>>from opt_four_segment_final import four_seg
    >>>four_seg(T,kappa,chi,n_max,four_seg_pulse_pars)
    ...
    ### for documentation use help()
    >>>help(opt_four_segment_final.four_seg)

    Parameters 
    ----------
    T_p: float 
        total pulse duration (µs)

    kappa: float 
        linewidth of the resonator (MHz)

    chi: float
        effective dispersive shift (MHz)

    n_p: float
        Maximum instant mean photon number

    initial_pars:
        Parameter vector for otimization,
        with duartions d1 to d4 and amplitudes a1 to a4 of the pulse,
        [d1,d2,d3,d4,a1,a2,a3,a4]
    
    

    
    Attributes
    ----------
    T, kappa, deltar, chi, kerr, n_p: float
        physical paramters (kerr: default value= 0.0)
        
    chi2: float
        dispersive shift for |e> state
       
    t_span: tuple
        time intervall (0,6*T)
      
    t_evaluate: ndarray
        number of evaluation points =100000 for the time intervall t_span
    
    cutoff_des: float
        upper integration limit for the instant state separation for t>= T*1
        
    guessed: ndarray
        Initial paramters (durations/T, ampliutdes/amp_max)
     
    res: ndarray
        Optimized Four-Segment Pulse parameters
        
    opt_square: ndarray
        Corresponding parameters of the Square Pulse (duration, amplitude)
    
    
    
       
    Methods
    -------
    optout(dt)
        Executes the optimization for the Four-Segment Pulse by minimizing the 'snr_ratio'
        the optimized parameters are stored in 'res' and the parameters of the Square Pulse in 'opt_square' 

    snr_ratio(dt)
        Computes the ratio between the integrated state separation after the total duration (T<= t<=t_cutoff)
        and the integrated state separation untill the end of the pulse (0<=t<=T)

    get_SNR_ratio(self,args,args2)
        Plotting of the instant state separation over time for the Four-Segment Pulse (args) 
        compared to the instant state separation for the Square Pulse(args2)
        and plotting of the instant mean photon number over time for the Four-Segment Pulse 
        and for the Square Pulse 
   
    plot_pulse()
        Plotting of the optimized Four-Segment Pulse (normalized by sqrt(kappa)) over time 

    plot()
        Plotting of the intra-cavity field trajectories for |g> and |e> state on complex plane
        for the optimized Four-Segment Pulse and as reference for the Square Pulse
    """  
    
    
    
    
    
    
    def square_pulse(self,t,arr):
        """
        Generates a square envelope pulse
        
        Parameters
        ---------
        t: np.ndarray
            time array
                  
                  
        pars: np.ndarray
            amplitude, duration
                           
        Returns
        -------
        drive: np.ndarray, dtype = float
            square envelope pulse
        """
        amp = arr[1]
        duration = arr[0]
        section = (np.heaviside(t,1)-np.heaviside((t-duration),1))
        drive = amp*section
        
        return drive
        
        
    def four_segments(self,t,params):
        """
        Generates a Four-Segment envelope pulse
        
        Parameters
        ---------
        t: np.ndarray
            time array
                  
                  
        pars: np.ndarray
            [t1,t2,t3,t4,amp1,amp2,amp3,amp4]
                           
        Returns
        -------
        pulse4: np.ndarray, dtype = float
            Four-Segment envelope pulse
        """
        first_segment = params[4]*(np.heaviside(t,1)-np.heaviside(t-params[0],1))
        second_segment = params[5]*(np.heaviside(t-params[0],1)-np.heaviside(t-(params[0]+params[1]),1))
        third_segment = params[6]*(np.heaviside(t-(params[0]+params[1]),1)-np.heaviside((t-(params[0]+params[1]+params[2])),1))
        fourth_segment = params[7]*(np.heaviside(t-(params[0]+params[1]+params[2]),1)-np.heaviside((t-(params[0]+params[1]+params[2]+params[3])),1))
        pulse4  = (first_segment+second_segment+third_segment+fourth_segment) 
        
        return pulse4
        
        
        
    # Inta-cavity field differential equation for |g> and |e> state Blais et al. Eq. (101)

    def diffae(self,t,b,params):
        """
        Inta-cavity field differential equation for |g> state Blais et al. Eq. (101)
        
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

        if (len(params)>2):
            drive = self.four_segments(t,params)
        else:
            drive = self.square_pulse(t,params)
        deltar =self.deltar   
        a = b[0]+1j*b[1]
        da_dt = drive-(1.0j*(deltar+self.chi)+(self.kappa/2))*a
        
        return [da_dt.real,da_dt.imag]
        

    def diffaem(self,t,b,params):
        """
        Inta-cavity field differential equation for |e> state Blais et al. Eq. (101)
        
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
        if (len(params)>2):
            drive = self.four_segments(t,params)
        else:
            drive =self.square_pulse(t,params)
        deltar=self.deltar   
            
        a = b[0]+1j*b[1]
        da_dt =drive-(1.0j*(deltar+self.chi2)+(self.kappa/2))*a
        
        return [da_dt.real,da_dt.imag]
        
        
        
           
        
        
    def cav_sol(self,dt):
        """
        Solves the Drive frame intra-cavity field differential equation 
                
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
        
        
        a0=[0.0,0.0] # initial condition of alpha : alpha_(i,j)(t=0) =0+0j
        cav_g = solve_ivp(self.diffae,self.t_span,a0,t_eval =self.t_evaluate ,args=(dt,))
        cav_e = solve_ivp(self.diffaem,self.t_span,a0,t_eval =self.t_evaluate ,args=(dt,))
        alpha_g = cav_g.y[0]+1j*cav_g.y[1]
        alpha_e = cav_e.y[0]+1j*cav_e.y[1]
        return cav_g.t, alpha_g, alpha_e
        
        
    # fnding the index for t=t_stamp (T is just an instance in time)
    def get_idx_T(self,t,t_stamp):
        return np.searchsorted(t, t_stamp, side='left')
    
    def snr_ratio(self,dt): 
        """
        Computes the ratio between the integrated state separation after the total duration (T<= t<=t_cutoff)
        and the integrated state separation untill the end of the pulse (0<=t<=T)
        Paramters
        ---------
        dt: np.ndarray, dtype = float
            kappa, deltar, K, chi,chi2
              
        Returns
        -------
        R: float
            Integrated state separation ratio ([Int. state sep. for t in (T,t_cutoff)] / [Int. state sep. for t in (0,T)])
        """
        T = self.T
        kappa = self.kappa
        cutoff_des = self.cutoff_des
        t, a_g, a_e = self.cav_sol(dt) # intra-cavitiy solution for both states
        discr =self.get_idx_T(t,T) # index where t=T
        cutoff = self.get_idx_T(t,cutoff_des)# index where t=cutoff_desired
        det = (t[1]-t[0])
        det2 = det 
        SNR_rungup = (np.sum(np.abs(a_g[:discr]-a_e[:discr])**(2)*det))
        
        SNR_afterT = (np.sum(np.abs(a_g[discr:cutoff]-a_e[discr:cutoff])**(2)*det2))  
        R =(SNR_afterT/(SNR_rungup))
        return R  
  
    def get_n_int(self,dt,tf):
        """
        Computes the integrated mean photon number over time t in (0,t_final)
       
        Paramters
        ---------
        dt: np.ndarray, dtype = float
            kappa, deltar, K, chi,chi2
               
        tf: float
            final time of integration

        Returns
        -------
        n_int: float
            integrated mean photon number
        """
        T = tf
        t, ag, ae = self.cav_sol(dt)
        discr =self.get_idx_T(t,T)
        Dt = (t[1]-t[0])
        n_int = self.kappa*np.sum((np.abs(ag[:discr])**(2)+np.abs(ae[:discr])**(2))/2)*Dt
        return n_int
        
       
        
    def n_eq(self,params,tf,n_tot):
        """
        Computes the total pulse duration to reach a certain
        integrated mean photon number 'n_tot' for given pulse parameters 'params'
        
        Parameters
        ----------
        params: 
        
        """        
        def residual(params):
            n_int = self.get_n_int(params,tf)
            return (n_tot-n_int)**(2)
            
        r = sc.optimize.minimize(residual,params)
        
        return r.x
        
    
    
    
    
        

    def get_SNR_ratio(self,args,args2):
        """
        Plotting of the instant state separation over time for the Four-Segment Pulse (args) 
        compared to the instant state separation for the Square Pulse(args2)
        and plotting of the instant mean photon number over time for the Four-Segment Pulse 
        and for the Square Pulse 
        """
        kappa = self.kappa
        T = self.T
        t, a_g, a_e = self.cav_sol(args)
        t, a_g2, a_e2 = self.cav_sol(args2)
        discr =self.get_idx_T(t,T)
        discr2 = self.get_idx_T(t,(args[0]+args[1]))
        cutoff = self.get_idx_T(t,self.cutoff_des)
        det = (t[1]-t[0])
        
        
        def _inst_state_sep(g,e):
            sep = np.abs(g-e)**(2)
            return sep
        def _inst_n_mean(g,e):
            n_m = ((np.abs(g)**(2)+np.abs(e)**(2))/2)
            return n_m
            
            
        inst_state_sep_four_seg = _inst_state_sep(a_g,a_e)
        inst_state_sep_square = _inst_state_sep(a_g2,a_e2)
        inst_n_mean_four_seg = _inst_n_mean(a_g,a_e)
        inst_n_mean_square = _inst_n_mean(a_g2,a_e2)
        max_inst_n_mean_four_seg = np.max((np.abs(a_g)**(2)+np.abs(a_e)**(2))/2)
        
        
        

        plt.plot(t,inst_state_sep_four_seg,label='Four-Segment',color=(254/255,178/255,76/255))
        plt.plot(t,inst_state_sep_square,label='Square',color='black')
        #plt.plot(t[discr2:],np.abs(a_g2[discr2]-a_e2[discr2])**(2)*np.exp(-kappa*(t[discr2:]-self.T)),linestyle='dashed',label='exp. decay',color='black')
        plt.vlines(self.T,0,1e6,linestyle='dashed',color='green',label=fr'$\rm{{T}} = {self.T}\,\rm{{µs}}$')
        plt.xlabel("t (µs)",fontsize=18)
        plt.xlim(0.*T,T*1.5)
        plt.ylabel(r" $| \alpha_{e}(t) - \alpha_{\rm{g}}(t) |^{2} \ \left(\bar{\rm{n}}\right)$",fontsize =18)
        plt.yscale('log')
        plt.ylim(2e-3,1e1)
        plt.title(r"Instantaneous State Separation",fontsize=18)
        plt.legend()
        plt.grid(True,'major')
        plt.tick_params(axis='both',which='major',labelsize=18)
        #plt.savefig("Plot_graphics/four_seg_state_sep_opt_IQ_HQCuing",dpi=300)
        plt.show()



        plt.plot(t,inst_n_mean_four_seg,label='Four-Segement')
        plt.plot(t,inst_n_mean_square,label='Square')
        plt.vlines(self.T,0,1e6,linestyle='dashed',color='green',label=fr'$\rm{{T}} = {self.T}\,\rm{{µs}}$')
        plt.hlines(0.01523,0.1*args[0],T+.379,linestyle='dashed',color='brown',label=fr'Thermal population $\mathrm{{n}}_{{0}} = {0.0015}$')
        plt.xlabel(r"t (µs)",fontsize=18)
        plt.xlim(0.1*args[0],T*1.3)
        plt.ylabel(r" $\left(| \alpha_{e}(t)|^{2} + |\alpha_{\rm{g}}(t) |^{2}\right)/2 \ \left(\bar{\rm{n}}\right)$",fontsize=18)
        plt.yscale('log')
        plt.ylim(5e-3,max_inst_n_mean_four_seg*1.2)
        plt.hlines(max_inst_n_mean_four_seg,0.,1.6*T,label=f'n_max <= {max_inst_n_mean_four_seg}'
        ,linestyle='dashed')
        plt.title(r"Instantaneous mean photon number",fontsize=18)
        plt.legend(loc='best',bbox_to_anchor =(1.35,0.5),fontsize=18)
        plt.grid(True,'major')
        plt.tick_params(axis='both',which='major',labelsize=18)
        plt.show()




        n_t = self.get_n_int(args,T)
        n_t2 = self.get_n_int(args2,T)
        n_res2 = np.sum((np.abs(a_g2[discr2:])**(2)+np.abs(a_e2[discr2:])**(2))/2.)*kappa*det
        n_res = np.sum((np.abs(a_g[discr:])**(2)+np.abs(a_e[discr:])**(2))/2.)*kappa*det
        print(f"n_int pulse_4_seg = {n_t}")
        print(f"n_int pulse_square = {n_t2}")
        print(f"n residual_4_seg = {n_res}")
        print(f"n_int residual_square = {n_res2}")
        SNR_rungup = np.sqrt(2*kappa*np.sum(np.abs(a_g[:discr]-a_e[:discr])**(2)*det/len(t[:discr])))
        SNR_afterT = np.sqrt(2*kappa*np.sum(np.abs(a_g[discr:cutoff]-a_e[discr:cutoff])**(2)*det))
        SNR_rungup2 = np.sqrt(2*kappa*np.sum(np.abs(a_g2[:discr2]-a_e2[:discr2])**(2)*det/len(t[discr2:])))
        SNR_afterT2 = np.sqrt(2*kappa*np.sum(np.abs(a_g2[discr2:cutoff]-a_e2[discr2:cutoff])**(2)*det))
        SNR_ratio = (SNR_afterT/SNR_rungup)
        SNR_ratio2 = (SNR_afterT2/SNR_rungup2)
        imp = SNR_ratio2/SNR_ratio #improvement factor

        print("Int. SNR from 0 to T four segmenself:",SNR_rungup)
        print("Int. SNR from 0 to T square:",SNR_rungup2)
        print("Int. SNR/n_res from 0 to T four segment:",SNR_rungup/n_res)
        print("Int. SNR/n_res from 0 to T square:",SNR_rungup2/n_res2)
        #print("Int. SNR from T to cutoff three segmenself:",SNR_afterT)
        #print("Int. SNR from T to cutoff square:",SNR_afterT2)
        #print("SNR ratio three segmenself:",SNR_ratio)
        #print("SNR ratio square:",SNR_ratio2)
        #print("Improvementby factor of:",imp)
        #print(f"Improvement by : {20*np.log10(imp)} dB")
        print(f" Residual photon number  {(n_res/n_res2)*100:.3f} % of residual photons of a square pulse ")
        



    def optout(self,dt):
        """
        Executes the optimization for the Four-Segment Pulse by minimizing the 'snr_ratio'
        the optimized parameters are stored in 'res' and the parameters of the Square Pulse in 'opt_square'
        The scipy module optimize.minimize is used with the method 'trust-constr'
        """
        T =self.T

        
        def cost(dt):

            dt_new = self._pass_initials(dt)

            return self.snr_ratio(dt_new)        
    
        

        def cons2(dt): #max instant mean photon number bounded to n_max (see above)

            dt_new = self._pass_initials(dt)
            t, ag, ae = self.cav_sol(dt_new)
            discr =self.get_idx_T(t,T)
            Dt = (t[1]-t[0])
            n_c = np.max((np.abs(ag[:discr])**(2)+np.abs(ae[:discr])**(2))/2.)
       
            return (n_c**(20))**(1/20) 
    
        def cons3(dt):

            dt_new = self.pass_initials(dt)
            T = self.T
            kappa = self.kappa
            cutoff_des = self.cutoff_des
            t, a_g, a_e = self.cav_sol(dt_new) # intra-cavitiy solution for both states
            discr =self.get_idx_T(t,T) # index where t=T
            cutoff = self.get_idx_T(t,cutoff_des)# index where t=cutoff_desired
            det = (t[1]-t[0])
            det2 = det 
            SNR_rungup = (np.sum(np.abs(a_g[:discr]-a_e[:discr])**(2)*det)) #Int. SNR from 0 to T    
            
            return SNR_rungup
            
            
            
            
       
 
        constraints = [NonlinearConstraint(cons2,self.n_p*0.99,self.n_p)]#,NonlinearConstraint(cons3,.5,np.inf)]
        dopt = np.append(np.array(dt[:2]),np.array(dt[2:]))
        print(dopt)
 
        lim = [((0.004/self.T) ,1.),((0.004/self.T) ,1.),((0.004/self.T) ,1.),(None,None),(None,None),(None,None),(None,None)] 
        result = sc.optimize.minimize(cost,dopt,method='trust-constr',constraints = constraints,options={ 'initial_tr_radius':1.0},callback=None)
        result.success
        print(result)
                
        

      
        result.x = self._pass_initials(result.x)
        n_tot = self.get_n_int(result.x,T)
        opt_square = np.concatenate([[result.x[0]+result.x[1]],[result.x[5]]])
        #self.get_SNR_ratio(result.x,opt_square)
        self.opt_square = opt_square
        self.result_opt = result.x
        print(f"Time segmenself: {(result.x[:4]/T)}")
        self.time_seg = (result.x[:4]/T)
        #print(f"Time segment: {self.time_seg}")
        maxamp = np.max(result.x[4:])
        #discr =self.get_idx_T(self.t,T)
        #n_t = np.sum((np.abs(a_g)**(2)+np.abs(a_e)**(2))/2.)*(t[1]-t[0])
        self.amp_n_est = np.sqrt(n_tot/2200)
        self.amp_seg_max =np.max(np.abs(result.x[4:]))
        self.amp_seg = (result.x[4:]/np.max(np.abs(result.x[4:])))
        print(f"amp segment: {self.amp_seg}")
        return result.x


       







 






   

   
    def plot_pulse(self):
        plt.figure(figsize=(4.5,2.75))
        #plt.plot(self.t_evaluate,self.square_pulse(self.t_evaluate,self.opt_square)/self.kappa,linestyle='dashed',label="Square",color='black')
        plt.plot(self.t_evaluate,self.four_segments(self.t_evaluate,self.res)/self.kappa,color=(254/255,178/255,76/255))
        #plt.vlines(self.T,1.4*np.min(self.res[4:]/self.kappa),1.1*self.res[4]/self.kappa,linestyle='dashed',color= 'maroon',label = fr"$T = {self.T}$ µs ")
        plt.xlim(0.0,self.T*1.1)
        plt.grid(True,'major')
        #plt.legend(fontsize=18)
        plt.xlabel("t (µs)",fontsize=18)
        plt.ylabel(r"$\rm a_{\rm in}/ \sqrt{\kappa}$",fontsize=18)
        #plt.title(r"Four segment pulse and square pulse",fontsize=18)
        plt.tick_params(axis='both',which='major',labelsize=18)
        plt.rcParams['lines.linewidth'] =3.0
        #plt.savefig("Plot_graphics/_four_seg_pulse_opt_IQ_HQCuing.svg",bbox_inches='tight')
        plt.show()
        """
        N= len(self.t_evaluate)
        spec_pulse = np.fft.fft(self.four_segments(self.t_evaluate,self.res))
        freq_pulse = np.fft.fftfreq(N,self.T)
        mask = freq_pulse >=0
        plt.plot(freq_pulse[mask],abs(spec_pulse[mask])*(2/N))
        plt.xlabel("f (MHz)")
        plt.ylabel(r"FFT($\mathrm{a}_{\mathrm{a_in}}$)")
        plt.xlim(0,0.1)
        plt.show()
        """
    def plot(self):
        plt.figure(figsize=(3.,3.))
        plt.plot(self.cav_sol(np.asarray(self.res))[1].real,self.cav_sol(np.asarray(self.res))[1].imag,label="|g>",color='b')
        plt.plot(self.cav_sol(np.asarray(self.res))[2].real,self.cav_sol(np.asarray(self.res))[2].imag,label="|e>",color='r')
        plt.plot(self.cav_sol(np.asarray(self.opt_square))[1].real,self.cav_sol(np.asarray(self.opt_square))[1].imag,linestyle='dashed',color='b')
        plt.plot(self.cav_sol(np.asarray(self.opt_square))[2].real,self.cav_sol(np.asarray(self.opt_square))[2].imag,linestyle='dashed',color='r')
        #plt.scatter(self.cav_sol(np.asarray(self.opt_square))[1][discr_n3].real,self.cav_sol(np.asarray(self.opt_square))[1].imag,label="square g")
        #plt.scatter(self.cav_sol(np.asarray(self.opt_square))[2].real,self.cav_sol(np.asarray(self.opt_square))[2].imag,label="square e")
        plt.xlabel(r"$\Re{\left(\alpha\right)}$",fontsize=18)
        plt.ylabel(r"$\Im{\left(\alpha\right)}$",fontsize=18)
        plt.legend(fontsize=18)
        plt.grid(True,'major')
        #plt.title(r"Intra-Cavity Field $\alpha(t)$" + "\n" + "on I-Q plane",fontsize=18)
        plt.rcParams['lines.linewidth'] =3.0
        #plt.savefig("Plot_graphics/_four_seg_alpha_opt_IQ_HQCuing.svg",bbox_inches='tight')
        plt.show()
     
    def __init__(self,T_p,kappa,chi,n_max,initial_guess):          #initializing,total pulse length,kappa,detuning,photon_number_callibration value 
        """
        Parameters
        ----------
        T_p: float 
            total pulse duration (µs)

        kappa: float 
            linewidth of the resonator (MHz)

        chi: float
            effective dispersive shift (MHz)

        n_p: float
            Maximum instant mean photon number

        initial_pars:
            Parameter vector for otimization,
            with duartions d1 to d4 and amplitudes a1 to a4 of the pulse,
            [d1,d2,d3,d4,a1,a2,a3,a4]
        """



        
        
        self.seg_dict = {
            'duration': T_p,
            'kappa': kappa,
            'detuning': 0.0,
            'chi': chi,
            'kerr':0.0,
            'n_p': n_max,
            'initial_pars':initial_guess,
        }
        self.T, self.kappa, self.deltar, self.chi, self.kerr,self.n_p,self.initial_guess = list(self.seg_dict.values())
        self.chi2 = -self.chi
        self.t_span = (0,(6*self.T)) # time range
        self.t_evaluate = np.linspace(*self.t_span,100000)
       
        T = self.T
        
      
        self.cutoff_des = self.t_evaluate[-1] # desired cutoff instance in time 
        self.guessed = self._initial_rel(self.initial_guess)

        self.res =self.optout(self.guessed)
        self.get_SNR_ratio(self.result_opt,self.opt_square)
        self.plot_pulse()
        self.plot()
        
        
        
    
        
        
        
        
        
        
        
    def int_snr_eq_square(self,arg,INT_SNR):    
        kappa = self.kappa
        def mod_snr_eq(arg):
            x =float(arg)
            args = np.array([x,self.res[3]],dtype='float')
            t, a_g, a_e = self.cav_sol(args)
            discr =self.get_idx_T(t,x)
            det = (t[1]-t[0])
            SNR_rungup_sq = np.sqrt(2*kappa*np.sum(np.abs(a_g[:discr]-a_e[:discr])**(2)*det/len(t[:discr])))
            return SNR_rungup_sq
            
        def res_snr_int_sq(arg):
            x= float(arg)
            int_snr = mod_snr_eq(x)
            return (INT_SNR-int_snr)
        r = sc.optimize.root_scalar(res_snr_int_sq,x0=arg,x1=arg*1.01)
        
        return r.root
        
        
    def desired_amp_add(self,dur,n_des,tf):
        def func_n(A):
                return (self.get_n_int((dur,A),tf)-n_des)
        r = sc.optimize.root_scalar(func_n,(0.0,self.n_p),method='brentq')
        return r.root    
        
        
        
        
        
        
        
        
        
    def _pass_initials(self,x):
        """
        Takes eight initial pulse parameters, 
        normalizes the time segments to the total pulse duration T 
        and the amplitudes to the maximal amplitude,
        and returns the free pulse paramters for the optimizer
            
        Parameters
        ----------
        x: ndarray, dtype=float
            inital pulse parameters
        
        Returns
        -------
        dopt_extracted: np.ndarray, dtype = float
            inital paramters for the optimizer,
            the last time segment is not a free parameter,
            since it is given by T and the durations of the previous
            time segments 
        """
        d1,d2,d3 = x[:3]*self.T
        d4 =(self.T - (d1+d2+d3))
        a1,a2,a3,a4 = x[3:]*np.max(self.initial_guess[5:])
        pass_pars = np.concatenate([[d1,d2,d3,d4],[a1,a2,a3,a4]])
        return pass_pars
        
    def _initial_rel(self,x):
        """
        Takes the free parameters, 
        scales the time segments with the total pulse duration T 
        and the amplitudes wtih the maximal amplitude,
        and returns all pulse paramters for the Four-Segment Pulse
            
        Parameters
        ----------
        x: ndarray, dtype=float
            free parameters
        
        Returns
        -------
        pass_pars: np.ndarray, dtype = float
            Four-Segment pulse parameters
        """
        t1,t2,t3 = (x[:3]/self.T)
        a1,a2,a3,a4 = (x[4:]/np.max(x[5:]))
        dopt_extracted = np.concatenate([[t1,t2,t3],[a1,a2,a3,a4]])
        return dopt_extracted
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    