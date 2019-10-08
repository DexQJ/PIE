# -*- coding: utf-8 -*-
"""
Created on Mon May 13 11:29:52 2019

@author: guila
"""
from openmdao.api import Group, LinearBlockGS, NonlinearBlockGS, NewtonSolver
from MODULES.MASS_STRUCTURE.mass_structure_stage1 import mass_structure1, mass_propellant1 
class mass_group(Group):
        def setup(self):
            self.add_subsystem('mass_structure1', mass_structure1() , promotes_inputs=['*'], promotes_outputs=['*'])
            self.add_subsystem('mass_propellant1', mass_propellant1() , promotes_inputs=['*'], promotes_outputs=['*'])
            #self.nonlinear_solver = NonlinearBlockGS(maxiter =500)#  NewtonSolver(maxiter =1000) #
            #self.nonlinear_solver.options['rtol'] = 10**-4
            #self.linear_solver = LinearBlockGS(maxiter =1000)
