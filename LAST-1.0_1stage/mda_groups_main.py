# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 15:33:26 2019

@author: guila
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from openmdao.api import Problem, Group, IndepVarComp, DirectSolver
from openmdao.api import NewtonSolver, NonlinearBlockGS, ScipyKrylov, LinearBlockGS, view_model, ExecComp
import result_vizualization


from MODULES.PROPULSION.propulsion_group import propulsion_group 
from MODULES.MASS_STRUCTURE.mass_group import mass_group
from MODULES.AERODYNAMICS.aerodynamics import aerodynamics
from MODULES.TRAJECTORY.trajectory import trajectory



class LAST_MDA(Group):
    """
    Group containing the Last MDA.
    """
    def setup(self):
        
        indeps = self.add_subsystem('indeps', IndepVarComp(), promotes=['*'])
        
        # variables corresponding to a RS-68 engine of a Delta IV-M first stage
        indeps.add_output('pc1',  102.6 * 100000)               # chamber pressure (Pa)
        indeps.add_output('pe1', 57525)                         # exit pressure (Pa)
        indeps.add_output('OF1', 5.97)                          # mixture ratio
        indeps.add_output('f1', 0.118)                          # inert fraction f = md / (md + mp)
 #       indeps.add_output('TW1', 1.185073797)                   # initial Thrust to weight ratio
        indeps.add_output('TW1', 1.22521781)                   # initial Thrust to weight ratio

        indeps.add_output('C1',0.2289697)                       # surface ratio between first stage and nozzle exit
        
        indeps.add_output('dt_v', val=25.)                # duration of vertical phase (s) 
        indeps.add_output('dt_po', val=20.)               # duration of pich over (s) 
        indeps.add_output('dthe_po', val=1)              # pitch over angle variation (°) 
        indeps.add_output('dt_d', val=20.)                # duration of angle variation towards gravity turn (s) 


        # Disciplines of the problem 
        self.add_subsystem('propulsion_group', propulsion_group() ,  promotes_inputs=['*'], promotes_outputs=['*'])
        self.add_subsystem('aerodynamics', aerodynamics() , promotes_outputs=['*'])
        cycle = self.add_subsystem('cycle', Group(), promotes=['*'])
        cycle.add_subsystem('mass_group',mass_group() , promotes_inputs=['*'], promotes_outputs=['*'])
        cycle.add_subsystem('trajectory',trajectory() , promotes_inputs=['*'], promotes_outputs=['*'])

        # Linear and non-linear Solvers for MDA
        cycle.nonlinear_solver =    NonlinearBlockGS(maxiter =100) #NewtonSolver(maxiter =100)# 
        cycle.nonlinear_solver.options['rtol'] = 10**-6
        cycle.nonlinear_solver.options['atol'] = 10**-6
        cycle.linear_solver = DirectSolver()# LinearBlockGS(maxiter =100) #ScipyKrylov()
        self.linear_solver = DirectSolver()# LinearBlockGS(maxiter =100) #ScipyKrylov()
        
        
# test of the MDA
test = Problem()
test.model = LAST_MDA()
test.set_solver_print(level=2)
test.setup()
test.run_model()
view_model(test, outfile="mda2.html", show_browser=False)   # N² diagramm


# generate as file for OpenVSP
view_values =    True # False #
plot_trajectory = True # False #
MDA_value = True

# visualize output values
if plot_trajectory == True:
    result_vizualization.plot_trajectory(MDA_value)
if view_values == True:
    result_vizualization.view_values(test)
