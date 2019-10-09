# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 16:00:32 2019

@author: guila
"""

from openmdao.api import  ExplicitComponent


class aerodynamics(ExplicitComponent):
    
    def setup(self):
        # declare outputs
        self.add_output('c_d', val=0.5)                     # Drag coefficient
        
        # Define analytically all partials.
        self.declare_partials('*', '*', method ='fd')
        
        # Counter of module execution
        self.execution_count = 0
        
    def compute(self, inputs, outputs):
        #outputs of the module
        c_d = 0.5
        outputs['c_d'] = c_d
        self.execution_count += 1
    
