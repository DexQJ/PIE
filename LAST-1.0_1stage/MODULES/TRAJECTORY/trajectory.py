# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 11:44:02 2019

@author: guila

revised on Mer Aug 28 
@author: aurbano
"""


from openmdao.api import ImplicitComponent, ExplicitComponent, Problem, NewtonSolver, IndepVarComp, DirectSolver, ScipyKrylov, view_model
import numpy as np
from MODULES.TRAJECTORY.EOM_calculations import EOM, nomass, hit_ground
from MODULES.TRAJECTORY.command_law import command_law_1, command_law_2, command_law_3, command_law_4, command_law_5
import scipy as sp
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import constants as cte


class trajectory(ExplicitComponent):

    def setup(self):
    # declare inputs
        self.add_input('m_p_1', val=1.)                 # propellant mass (kg) 
        self.add_input('Isp1', val=1.)                  # specific impulse (m/s)
        self.add_input('m_d_1', val=1)                  # dry mass (kg)
        self.add_input('A_e1', val=1)                   # nozzle exit surface (m²)
        self.add_input('TW1', val = 2)                  # initial mass to weigth ratio
        self.add_input('c_d', val = 0.5)                # drag coefficient
        self.add_input('C1', val = 0.22869)             # surface ratio between first stage and nozzle exit
        self.add_input('dt_v', val=100.)                # duration of vertical phase (s) 
        self.add_input('dt_po', val=100.)               # duration of pich over (s) 
        self.add_input('dthe_po', val=90.)              # pitch over angle variation (°) 
        self.add_input('dt_d', val=100.)                # duration of angle variation towards gravity turn (s) 

        # declare outputs
        self.add_output('m_dot1', val=892.7921706374962)    # flow rate(kg/s)
        self.add_output('m_dot2', val=35.909266843036235)               # flow rate(kg/s)
        self.add_output('alt_f', val= 100)                  # final altitude (km)
        self.add_output('GLOW', val= 305.86586438)          # Gross Lift Off Weight (t)
        self.add_output('v_f', val= 1000)                   # final speed (m/s)
        self.add_output('A_s_1', val= 18.994020306876973)               # stage surface first stage (m²)
        self.add_output('A_s_2', val=2.7951690059015513)                # stage surface first stage (m²)
        self.add_output('ft1', val=  3554.65193 * 1e3)      # Thrust (N)
        self.add_output('ft2', val= 142.9727423446221 * 1e3)        # Thrust second stage (kN)
        self.add_output('fi_f', val= 1)                     # final flight path angle
        self.add_output('l_f', val= 1000)                   # final longitude
        self.add_output('alf_max', val= 5)                  # maximum incidence (deg)
        self.add_output('nx_max',5.053814)                  # maximal axial load factor 
        self.add_output('Pdyn_max',2)                                   # maxQ (kPa)       

        self.declare_partials('*', '*', method='fd')
        
        # Counter of module execution
        self.execution_count = 0

    def compute(self, inputs, outputs):
        m_p1 = inputs['m_p_1']
        m_d1 = inputs['m_d_1']
        m_p2 = cte.m_p_2
        m_d2 = cte.m_d_2
        Isp1 = inputs['Isp1']
        A_e1 = inputs['A_e1']
        TW1 = inputs['TW1']
        c_d = inputs['c_d']
        C1 = inputs['C1']
        dt_po = inputs['dt_po']
        dt_v = inputs['dt_v']
        dthe_po = inputs['dthe_po']*np.pi/180
        dt_d = inputs['dt_d']
        
        # module parameters
        g0 = cte.g0
        m_pl = cte.m_pl                                             # payload mass(kg)
        omega =2*sp.pi/ cte.earth_rotation_time                     # Rotation rate of the Earth (rad/s) 23h 56 min 4.0905 s
        R = cte.RT                                                  # earth radius (m)
#        rho0 = cte.rho0                                            # density at r=0m
        A_stage1 =   A_e1 / C1 * cte.n_engines                      # area of the stage
        A_stage2 =   cte.A_stage_2                                                                    # standard gravitational parameter (m^3/s^2)
        
        # initial state conditions for the EOM
        r0 = R  + 100                                                   # Initial radial position = radius of Earth
        l0 = cte.longitude0 * sp.pi / 180                           # Initial angular position (Kourou)
        v0 = 10.00                                                  # Initial velocity = 0
        fi0 =sp.pi/2                                                # Initial flight path angle
        m02 = m_p2 + m_d2 + m_pl
        m0 = m_p1 + m_d1 + m02
 #       m0 = m_p1 + m_d1 + m_p2 + m_d2 + m_pl                            # Initial mass of the rocket (kg)
        x0 = [r0, l0, v0, fi0, m0]                                  # Initial state
        #print(m_p2,m_d2)
        # set initial Thrust and integration parameters

        ft1 = TW1 * m0 * g0                                           # thurst first stage (N)
        c1 = Isp1 * g0                                                # Effective exhaust velocity (m/s)
        m_dot1 = ft1 / c1 / cte.n_engines                             # flow rate of one engine (kg/s)
        Thrust_time1 = m_p1 /( m_dot1 * cte.n_engines)                # Thrust until you consume all m_p (need to know m_p for this)
        int_time = Thrust_time1                                       # integration time
        
              
 #       theta = lambda t: sp.pi/2 * (-t / t_a +1)                     # pitch angle
        ft = lambda t: ft1*sp.heaviside(Thrust_time1-t,0)          # Thrust magnitude
        

        # Integration parameters 
        atol_integration =1e-12 # The solver keeps the local error estimates less than atol + rtol * abs(y)
        rtol_integration=1e-9 
        integration_method = 'BDF' # Implicit method based on backward-differentiation formulas (order 5): sum(s,k=0){a_k*y_{n+k}} = h*\beta*f(t_{n+s},y_{n+s})
#        step=1.
        span_integration = (0.,int_time)
 
       # solving EOM
        integrate = True

        args = (omega, ft, c1, A_stage1,c_d,integrate,dt_po,dt_v,dthe_po,dt_d)
        sol = solve_ivp(lambda t, x: EOM(t, x, *args), span_integration, x0,atol=atol_integration,rtol=rtol_integration,method = integration_method,dense_output=True,events=(hit_ground,nomass))
        
        # used solution
        t_sol = sol.t
        x = sol.y
        r_sol = x[0]
        l_sol = x[1]
        v_sol = x[2]
        fi_sol = x[3]
        m_sol = x[4]
        
         
        #-------------- get relevant data from the trajectory ----------------
        integrate = False # set to False to get only data
        args = (omega,ft, c1, A_stage1,c_d,integrate,dt_po,dt_v,dthe_po,dt_d)
       
        current_T = t_sol
        size_t=len(t_sol)
        current_NX = np.zeros(size_t)
        current_Mach = np.zeros(size_t)
        current_Pdyn = np.zeros(size_t)
        current_flux = np.zeros(size_t)
        current_rho = np.zeros(size_t)
        

        for i in range(len(t_sol)):
                (current_NX[i],
                 current_Mach[i],
                 current_Pdyn[i],
                 current_flux[i],
                 current_rho[i]) = EOM(current_T[i],sol.sol(current_T[i]), *args)

        
        nx_max = max(current_NX)

        theta = np.zeros([len(t_sol)])
        
        for i in range(len(t_sol)):
            theta[i] = command_law_5(dt_po,dthe_po,dt_v,dt_d,fi_sol[i],t_sol[i])
            
        # constraits for MDO
        alf = theta - fi_sol
        nx_max = max(current_NX)
        alf_max = max(alf)*180/sp.pi 
        Pdyn_max = max(current_Pdyn)
        h_sol = r_sol -R
        
        v_f = v_sol[-1]
        m_f = m_sol[-1]
        l_f = l_sol[-1]
        fi_f = fi_sol[-1]
        h_f = r_sol[-1] - R




        # Outputs trajectory to CSV file
        size = 15
        csvData = np.zeros((size_t,size))
        csvData[:,0] = t_sol
        csvData[:,1] = (r_sol - cte.RT)  / 1e3
        csvData[:,2] = v_sol
        csvData[:,3] = m_sol /1e3
        csvData[:,4] = l_sol * 180 / np.pi
        csvData[:,5] = fi_sol * 180 / np.pi
        csvData[:,6] = theta * 180 / np.pi
        csvData[:,7] = alf * 180 / np.pi
        csvData[:,8] = ft1 / 1e3 * np.ones(size_t)
        csvData[:,9] = current_Pdyn / 1e3 
        csvData[:,10] = current_Mach
        csvData[:,11] = current_flux / 1e3
        csvData[:,12] = current_rho
        csvData[:,13] = current_NX
        csvData[:,14] = c_d * np.ones(size_t)
        np.savetxt("trajectory_png/trajectory_file.csv", csvData, delimiter=",") 
        
        csvData = [None]*2
        csvData[0] = A_stage1
        csvData[1] = A_stage2
        np.savetxt("OpenVSP/surface_file.csv", csvData, delimiter=",")


        #output of the module

        outputs['GLOW'] = m0 / 1e3            # GLOW (t)
        outputs['ft1'] = ft1                  # Thrust (N)
        outputs['m_dot1'] = m_dot1            # flow rate single engine (kg/s)
        outputs['A_s_1'] = A_stage1           # stage surface
        outputs['nx_max'] = nx_max            # maximum axial load factor
        outputs['alf_max'] = alf_max          # maximum incidence
        outputs['alt_f'] = h_f/1e3            # final altitude (km) = h inertial ref
        outputs['v_f'] = v_f                  # final speed (m/s) non inertial
        outputs['l_f'] = l_f*180/sp.pi        # final longitude (deg)
        outputs['fi_f'] = fi_f*180/sp.pi      # final flight path angle (deg)
        outputs['nx_max'] = nx_max
        outputs['Pdyn_max'] = Pdyn_max / 1e3 
        outputs['A_s_1'] = A_stage1
 #       outputs['A_s_2'] = A_stage2
 
         
 
if __name__ == "__main__":


    probl = Problem()
    indeps = probl.model.add_subsystem('indeps', IndepVarComp(), promotes=['*'])
    indeps.add_output('m_p_1', 198758.2907090509) # inert fraction f = md / (md + mp)
    indeps.add_output('Isp1', 365.4569598869584) 
    indeps.add_output('m_d_1', 26591.24524225397) 
    indeps.add_output('A_e1', 4.271520730041623) 
    indeps.add_output('TW1', 1.224521781) 
    indeps.add_output('c_d', 0.5) 
    indeps.add_output('C1', 0.2289697) 
#    indeps.add_output('dt_po', 20) 
#    indeps.add_output('dt_v', 25) 
#    indeps.add_output('dthe_po', 1) 
#    indeps.add_output('dt_d', 20)
    indeps.add_output('dt_po', 20) 
    indeps.add_output('dt_v', 25) 
    indeps.add_output('dthe_po', 0.5) 
    indeps.add_output('dt_d', 20)


    probl.model.add_subsystem('trajectory',trajectory(),promotes_inputs=['*'], promotes_outputs=['*'])
    
    probl.model.nonlinear_solver = NewtonSolver()
    probl.model.linear_solver = DirectSolver()
    #probl.model.linear_solver = ScipyKrylov()
    probl.model.nonlinear_solver.options['maxiter'] = 10**3
    probl.model.nonlinear_solver.options['rtol'] = 10**-6    
    probl.set_solver_print(level=1)
    
    probl.setup()
    
    probl.run_model()
    view_model(probl, outfile="trajectoryn2.html", show_browser=False)    
    test = probl
    print('\n')
    print('final altitude (km) =',test['alt_f'][0])
    print('final speed (m/s)=',test['v_f'][0])
    print('final  fi= ',test['fi_f'][0])
    print('final long= ',test['l_f'][0])
    print("\n")
    print("m_dot1:", test['m_dot1'][0])
#    print("m_dot2:", test['m_dot2'][0])
    
    print("Thrust 1 :", test['ft1'][0])
    print('nx_max ', test['nx_max'][0])
    print('alf_max ', test['alf_max'][0])
    
#"""
