# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 15:53:42 2019

@author: guila
"""

from openmdao.api import ExplicitComponent, Problem, IndepVarComp, LinearBlockGS
import numpy as np
from rocketcea.cea_obj import CEA_Obj

class propulsion1(ExplicitComponent):
    
    def setup(self):
        # declare inputs
        self.add_input('pc1', val=25. * 1e6)                            # chamber pressure (Pa)
        self.add_input('pe1', val=101325)                               # exit pressure (pa)
        self.add_input('OF1', val=3.)                                   # mixture ratio
        
        # declare outputs
        self.add_output('Isp1', val = 365.3321606193021)              # specific impulse (m/s)
        self.add_output('c_star1', val= 2290.5814135894693)           # combustion characteristic speed (m/s)
        self.add_output('epsilon1', val=21.50588729916232)            # expansion ratio
        
        # Finite difference all partials.
        self.declare_partials('*', '*', method='fd')
        
        # Counter of module execution
        self.execution_count = 0
        
    def compute(self, inputs, outputs):
        
        # inputs of the module
        pc = inputs['pc1']
        pe = inputs['pe1']
        OF = inputs['OF1']
        
        # Rocket_CEA (only need pc(psi) and OF, eps no influence for CEA_output)           
        CEA_output = CEA_Obj(fuelName='LH2',oxName='LOX')
        value = CEA_output.get_IvacCstrTc_ChmMwGam(Pc=pc * 0.000145038, MR=OF, eps=21.5)
        
        Tc =  value[2] * 5 / 9                      # chamber temperature (K rocket CEA gives it in deg Rankine)
        M = value[3]                                # molecular mass at combustion (kg/kmol)
        gam = value[4]                              # isentropic coefficient at throat (SI)
        #ispvac = value[0]                          # vacuum Isp (s) 
         
        # Module parameters
        R    = 8314/ M 
        g0   = 9.81       # gravitational acceleration (m/s²)
        pa   = 101325     # atmospheric pressure at z = 0m (Pa) 
        etac = 0.99       # combustion efficiency 
        etan = 0.938      # nozzle efficiency
        
        # Module calculations
        gamma = np.abs(np.sqrt(gam * (2/(gam + 1)) ** ((gam+1)/(gam-1))))
        
        epsilon = gamma / np.sqrt(2 * gam / (gam - 1) * 
                                  (pe / pc) ** (2/gam) * (1 - (pe/pc) ** ((gam - 1)/gam)))
        
        c_star = etac * np.sqrt(R * Tc) / gamma
        
        cf = gamma * np.sqrt(2 * gam / (gam - 1) * (1 - (pe / pc) **
                                        ((gam-1)/gam))) + epsilon / pc * (pe -pa)

        Isp =  etan*cf*c_star/g0
        
        # outputs of the module
        outputs['c_star1'] = c_star
        outputs['epsilon1'] = epsilon
        outputs['Isp1'] = Isp
        self.execution_count += 1
