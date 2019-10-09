# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 16:00:32 2019

@author: guila
"""

from openmdao.api import ExplicitComponent
from MODULES.MASS_STRUCTURE.mass_calculations import mass_additional, mass_tank, mass_press_gaz, mass_turbopump, mass_c_chamber, mass_nozzle
import numpy as np
class mass_structure1(ExplicitComponent):
    
    def setup(self):
        self.add_input('pc1', val=0.)            # chamber pressure (MPa)
        self.add_input('OF1', val=0.)            # mixture ratio 
        self.add_input('epsilon1', val = 0)      # expansion ratio
        self.add_input('c_star1', val=0.)        # combustion characteristic speed (m/s)
        self.add_input('m_dot1', val=100.)       # mass flow rate (kg/s)
        self.add_input('m_p_1', val = 1)
        self.add_input('nx_max',val = 4)
        
        self.add_output('m_d_1', val=1e4)         # dry mass (kg)
        self.add_output('A_t_1', val=0.01)        # throat surface (m²)
        self.add_output('A_e_1', val=5.)          # exit surface (m²)

        
        
        # Finite difference all partials.
        self.declare_partials('*', '*', method='fd')

        
    def compute(self, inputs, outputs):
        # inputs of module
        pc = inputs['pc1']
        OF = inputs['OF1']
        m_p = inputs['m_p_1']
        c_star = inputs['c_star1']
        m_dot= inputs['m_dot1']
        eps = inputs['epsilon1']
        nx_max = inputs['nx_max']       
        
        # module calculations
        
        # Throat surface (m²)
        A_t = c_star * m_dot / pc
        
        # exit surface (m²)
        A_e = eps * A_t
        

        # mass and volumes tanks  (kg,m³)
        tank = mass_tank(m_p, OF, A_e)
        m_t_lox = tank[0]
        m_t_h2 = tank[1]
        v_lox =  tank[2]
        v_h2 =  tank[3]
        p_lox = tank[3]
        p_h2 = tank[4]

        
        # mass turbopump (kg)
        m_tp = mass_turbopump(m_p, OF, pc, m_dot)
        m_tp_lox = m_tp[0]
        m_tp_h2 = m_tp[1]

        # mass pressurised gaz (kg)
        m_pg = mass_press_gaz(OF, m_p)
        m_pg_lox = m_pg[0]
        m_pg_h2 = m_pg[1]
        
        # mass propellant systems (kg)
        m_sys_lox = m_t_lox + m_tp_lox + m_pg_lox
        m_sys_h2 = m_t_h2 + m_tp_h2 + m_pg_h2
        
        # mass, lengths and Area combustion chamber(kg,m,m²)
        chamber = mass_c_chamber(pc, A_t)
        m_cc = chamber[0]
        l_cc = chamber[1]
        h_cone = chamber[2]
        A_c = chamber[3]   
        e_cc = chamber[4]   
        
        # mass length nozzle (kg,m)
        nozzle = mass_nozzle(pc, A_e, A_t)
        m_n = nozzle[0]
        L_n = nozzle[1]
        e_n = e_cc
        
        # mass engine (kg)
        m_eng = m_cc + m_n
        
        # additional mass (kg)
        m_add = mass_additional(m_p)
        
       # connect trajectory to dry mass
        n_max = nx_max
        c = 0.02*n_max + 1
        
        # dry mass (kg)
        m_d = c * (m_add + m_eng + m_sys_lox + m_sys_h2)
        
         # Outputs geometry to CSV file
        size = 18
        csvData = [None]*size
        csvData[0] =m_d
        csvData[1] = A_e
        csvData[2] =  A_t
        csvData[3] =  m_add
        csvData[4] = m_t_lox
        csvData[5] =  m_t_h2
        csvData[6] =  v_lox
        csvData[7] = v_h2
        csvData[8] = p_lox
        csvData[9] = p_h2
        csvData[10] = m_n
        csvData[11] = A_c
        csvData[12] = L_n
        csvData[13] = m_cc
        csvData[14] = l_cc
        csvData[15] = h_cone
        csvData[16] = e_n
        csvData[17] = e_cc
        np.savetxt("OpenVSP/geom_stage1_file.csv", csvData, delimiter=",") 
        """
        csvData[1] =
        csvData[1] =
        csvData[20] =
        csvData[2] =
        csvData[2] =
        csvData[2] =
        csvData[2] =
        csvData[2] =
        csvData[2] =
        csvData[2] =
        csvData[2] =
        csvData[2] =
        csvData[2] =
        """
        
        
        
        # output of module 
        outputs['m_d_1'] = m_d
        outputs['A_e_1'] = A_e
        outputs['A_t_1'] = A_t
        
        

class mass_propellant1(ExplicitComponent):
    
    def setup(self):
        self.add_input('m_d_1', val=1.)                         # dry mass (kg)
        self.add_input('f1', val =0.)                           # inert fraction
        self.add_output('m_p_1', val = 232174.800041)           # propellant mass (kg)
        
        # Finite difference all partials.
        self.declare_partials('*', '*', method='fd')
        
    def compute(self, inputs, outputs):
        # inputs of the module
        m_d = inputs['m_d_1']
        f = inputs['f1']
        
        # outputs of module 
        m_p = (1 - f) / f * m_d
        outputs['m_p_1'] = m_p

        
        
        
        
        
        