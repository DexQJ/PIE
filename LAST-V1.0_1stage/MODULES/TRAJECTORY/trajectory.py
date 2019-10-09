# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 11:44:02 2019

@author: guila
"""


from openmdao.api import ImplicitComponent, ExplicitComponent, Problem, NewtonSolver, IndepVarComp, DirectSolver, ScipyKrylov, view_model, NonlinearBlockGS
import numpy as np
from MODULES.TRAJECTORY.EOM_calculations import EOM, nomass, hit_ground
import matplotlib.pyplot as plt
#from EOM_calculations import EOM, nomass, hit_ground
import scipy as sp
from scipy.integrate import solve_ivp
import constants as cte

class trajectory(ExplicitComponent):

    def setup(self):
        # declare inputs
        self.add_input('m_p_1', val=1.)                         # propellant mass first stage (kg)
        self.add_input('Isp1', val=1.)                          # specific impulse first stage (s)
        self.add_input('m_d_1', val=1)                          # dry mass first stage (kg)
        self.add_input('A_e_1', val=1)                          # nozzle exit surface (m²)
        self.add_input('TW1', val = 2)                          # initial Thrust to weigth ratio
        self.add_input('c_d', val = 0.5)                        # drag coefficient
        self.add_input('C1', val = 0.22869)                     # surface ratio between first stage and nozzle exit
        self.add_input('t_a', val=100.)                         # pitch over time (s) 
        self.add_input('t_vertical', val=100.)                         # pitch over time (s) 
        
        # declare outputs
        self.add_output('m_dot1', val=892.7921706374962)                # flow rate(kg/s)
        self.add_output('m_dot2', val=35.909266843036235)               # flow rate(kg/s)
        self.add_output('alt_final', val= 100)                          # final altitude (km)
        self.add_output('GLOW', val= 305.86586438)                      # Gross Lift Off Weight (t)
        self.add_output('v_end', val= 1000)                             # final speed
        self.add_output('A_s_1', val= 18.994020306876973)               # stage surface first stage (m²)
        self.add_output('A_s_2', val=2.7951690059015513)                # stage surface first stage (m²)
        self.add_output('Thrust1', val=  3554.6519381139733 * 1e3)      # Thrust first stage (kN)
        self.add_output('Thrust2', val= 142.9727423446221 * 1e3)        # Thrust second stage (kN)
        self.add_output('final_gamma', val= 1)                          # final flight path angle (deg)
        self.add_output('final_vl', val= 1000)                          # final tangencial speed (m/s)
        self.add_output('final_vr', val= 100)                           # final radial speed (m/s)
        self.add_output('nx_max',5.053814709229969)                     # maximal axial load factor
        self.add_output('Pdyn_max',2)                                   # maxQ (kPa)
        self.add_output('alpha_max',2)                                  # maximum angle of attack (deg)
        
        # Finite difference all partials.
        self.declare_partials('*', '*', method='fd')
        
        # Counter of module execution
        self.execution_count = 0

    def compute(self, inputs, outputs):
        
        # declare inputs
        m_p1 = inputs['m_p_1']
        m_d1 = inputs['m_d_1']
        m_p2 = cte.m_p_2
        m_d2 = cte.m_d_2
        Isp1 = inputs['Isp1']
        A_e1 = inputs['A_e_1']
        TW1 = inputs['TW1']
        c_d = inputs['c_d']
        C1 = inputs['C1']
        t_a = inputs['t_a']
        t_vertical = inputs['t_vertical']
        
        # module parameters
        g0 = cte.g0
        m_pl = cte.m_pl                                               # payload mass(kg)
        omega =2*sp.pi/ cte.earth_rotation_time                       # Rotation rate of the Earth (rad/s) 23h 56 min 4.0905 s
        R = cte.RT                                                    # earth radius (m)
        A_stage1 =   A_e1 / C1 * cte.n_engines
        A_stage2 =   cte.A_stage_2
        mu = cte.mu                                                   # standard gravitational parameter (m^3/s^2)
        
        # initial state conditions for the EOM
        r0 = R                                                        # Initial radial position = radius of Earth
        l0 = cte.longitude0 * sp.pi / 180                             # Initial angular position (Kourou)
        vr0 = 0.1                                                    # Initial radial velocity = 0.1
        vt0 =0.00 # omega*R                                           # Initial tangential velocity = 0 rotation rate of the Earth
        m0 = m_p1 + m_d1 + m_p2 + m_d2 + m_pl                         # Initial mass of the rocket (kg)
        outputs['GLOW'] = m0 / 1e3
        x0 = [r0, l0, vr0, vt0, m0]                                   # Initial state
        m02 = m_p2 + m_d2 + m_pl
        
        # set initial Thrust and integration parameters
        thrust1 = TW1 * m0 * g0 
        c1 = Isp1 * g0                                                # Effective exhaust velocity (m/s)
        m_dot1 = thrust1 / c1 / cte.n_engines                         # flow rate of one engine (kg/s)
        Thrust_time1 = m_p1 /( m_dot1 * cte.n_engines)                # Thrust until you consume all m_p (need to know m_p for this)
        int_time = Thrust_time1
        
        
    
        #Thrust_time2 = m_p2 /( m_dot2 * cte.n_engines)               # Thrust until you consume all m_p (need to know m_p for this)
        # define thrust and thrust angle
        #theta = lambda t: sp.pi/2 * (-t / t_a +1)                     # Thrust angle wrt. inertial frame
        Thrust = lambda t: thrust1*sp.heaviside(Thrust_time1-t,0)     # Thrust magnitude
        

        # Integration parameters 
        atol_integration =1e-12 # The solver keeps the local error estimates less than atol + rtol * abs(y)
        rtol_integration=1e-9 
        integration_method = 'BDF' # Implicit method based on backward-differentiation formulas (order 5): sum(s,k=0){a_k*y_{n+k}} = h*\beta*f(t_{n+s},y_{n+s})
        step=1.
        span_integration = (0.,int_time)

        # solving EOM
        integrate = True
        args = (omega, Thrust, c1, A_stage1,c_d, mu,R,integrate,t_a, t_vertical)
        sol = solve_ivp(lambda t, x: EOM(t, x, *args), span_integration, x0,atol=atol_integration,rtol=rtol_integration,method = integration_method,dense_output=True,events=(hit_ground,nomass))
        
        # used solution
        t = np.append(np.arange(sol.t[0],sol.t[-1],step),sol.t[-1])
        t_sol = sol.t
        x = sol.y
        r_sol = x[0]
        l_sol = x[1]
        vr_sol = x[2]
        vl_sol = x[3]
        m_sol = x[4]
        
         
        #-------------- get relevant data from the trajectory ----------------
        
        integrate = False # set to False to get only data
        args = (omega, Thrust, c1, A_stage1,c_d, mu,R,integrate,t_a, t_vertical)
       
        current_T = t_sol
        current_NX = np.zeros([len(current_T)])
        current_Mach = np.zeros([len(current_T)])
        current_Pdyn = np.zeros([len(current_T)])
        current_flux = np.zeros([len(current_T)])
        current_rho = np.zeros([len(current_T)])
        current_theta = np.zeros([len(current_T)])

        for i in range(len(t_sol)):
                (current_NX[i],
                 current_Mach[i],
                 current_Pdyn[i],
                 current_flux[i],
                 current_rho[i],current_theta[i]) = EOM(current_T[i],sol.sol(current_T[i]), *args)
        

        
        # modulus of speed and flight path angle
        v_sol = np.zeros([len(vr_sol)])
        gamma_sol = np.zeros([len(vr_sol)])
        alpha_sol = np.zeros([len(vr_sol)])
        
        for i in range(len(vr_sol)):
            v_sol[i] = np.sqrt(vr_sol[i] ** 2 + vl_sol[i] ** 2)
            if vr_sol[i] == 0 or vl_sol[i] == 0 :
                gamma_sol[i] = np.pi/2
                alpha_sol[i] = 0
            else: 
                gamma_sol[i] = np.arctan(vr_sol[i] / vl_sol[i])
                alpha_sol[i] = gamma_sol[i] - current_theta[i]
                
        # Outputs trajectory to CSV file
        size = 16
        csvData = [None]*size
        csvData[0] =current_T
        csvData[1] = (r_sol - cte.RT)  / 1e3
        csvData[2] =  v_sol
        csvData[3] =  m_sol /1e3
        csvData[4] = gamma_sol * 180 / np.pi
        csvData[5] =  current_theta * 180 / np.pi
        csvData[6] =  current_Pdyn / 1e3
        csvData[7] = l_sol * 180 / np.pi
        csvData[8] = l_sol * 180 / np.pi
        csvData[9] = current_Mach
        csvData[10] = current_flux / 1e3
        csvData[11] = current_rho
        csvData[12] = c_d * np.ones(np.size(t_sol))
        csvData[13] = Thrust(t_sol) / 1e3
        csvData[14] = current_NX
        csvData[15] = alpha_sol * 180 / np.pi
        np.savetxt("trajectory_png/trajectory_file.csv", csvData, delimiter=",") 
        
        csvData = [None]*2
        csvData[0] = A_stage1
        csvData[1] = A_stage2
        np.savetxt("OpenVSP/surface_file.csv", csvData, delimiter=",")


        # constraits for MDO
        nx_max = max(current_NX)
        Pdyn_max = max(current_Pdyn)
        alpha_max = max(alpha_sol)
        
        ######################################################################
# ------------------ OUTPUTS OF THE MODULE -----------------------------------
        ######################################################################
        
        v_end = v_sol[-1]
        #m_end = m_sol[-1]
        h_end = r_sol[-1] -R

        
        outputs['alt_final'] = h_end/1e3         
        outputs['v_end'] = v_end                    
        outputs['Thrust1'] = thrust1
        outputs['m_dot1'] = m_dot1                
        outputs['final_gamma'] = np.arctan(vr_sol[-1] / vl_sol[-1])
        outputs['final_vr'] = vr_sol[-1]
        outputs['final_vl'] = vl_sol[-1]
        outputs['A_s_1'] = A_stage1
        outputs['A_s_2'] = A_stage2
        outputs['nx_max'] = nx_max
        outputs['Pdyn_max'] = Pdyn_max / 1e3 
        outputs['alpha_max'] = alpha_max * 180 / np.pi

        