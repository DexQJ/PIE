# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 15:33:26 2019

@author: guila
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from openmdao.api import Problem, Group, IndepVarComp, DirectSolver
from openmdao.api import NewtonSolver, NonlinearBlockGS, ScipyKrylov, LinearBlockGS, view_model, ExecComp, ScipyOptimizeDriver
import result_vizualization

from MODULES.PROPULSION.propulsion_group import propulsion_group 
from MODULES.MASS_STRUCTURE.mass_group import mass_group
from MODULES.AERODYNAMICS.aerodynamics import aerodynamics
from MODULES.TRAJECTORY.trajectory import trajectory


class LAST_MDO(Group):
    """
    Group containing the Last MDA.
    """

    def setup(self):
        
        indeps = self.add_subsystem('indeps', IndepVarComp(), promotes=['*'])   # these are the desgin variables
        
        # variables corresponding to the first stage of Delta IV-M rocket a RS-68 engine of 
        indeps.add_output('pc1',  102.6 * 100000, ref = 10000000)           # chamber pressure (Pa)
        indeps.add_output('pe1', 57525)                                     #  exit pressure (Pa)
        indeps.add_output('OF1', 5.97)                                      # mixture ratio
        indeps.add_output('f1', 0.118)                                      # inert fraction f = md / (md + mp)
        indeps.add_output('TW1', 1.185073797)                               # initial Thrust to weight ratio
        indeps.add_output('C1',0.2289697)                                   # surface ratio between first stage and nozzle exit
               
        indeps.add_output('t_a',250, ref = 100)                             # command law (time por pitch over)
        indeps.add_output('t_vertical',30, ref=10)                          # command law (time of vertical phase)
        
        # Disciplines of the problem 
        self.add_subsystem('propulsion_group', propulsion_group() ,  promotes_inputs=['*'], promotes_outputs=['*'])
        self.add_subsystem('aerodynamics', aerodynamics() , promotes_outputs=['*'])
        cycle = self.add_subsystem('cycle', Group(), promotes=['*'])
        cycle.add_subsystem('mass_group',mass_group() , promotes_inputs=['*'], promotes_outputs=['*'])
        cycle.add_subsystem('trajectory',trajectory() , promotes_inputs=['*'], promotes_outputs=['*'])

        # Linear and non-linear Solvers for MDA
        cycle.nonlinear_solver =NonlinearBlockGS(maxiter =50) #  NewtonSolver(maxiter =10)# 
        cycle.nonlinear_solver.options['rtol'] = 10**-4
        cycle.nonlinear_solver.options['atol'] = 10**-4
        cycle.linear_solver = DirectSolver()# LinearBlockGS(maxiter =100) #ScipyKrylov()
        
        
        # Objective function and constraints for the MDO
        self.add_subsystem('obj_cmp', ExecComp('obj = GLOW',
                           GLOW = 260),
                           promotes=['obj', 'GLOW'])
        
        self.add_subsystem('con_cmp1', ExecComp('con1 =  alt_final'), promotes=['con1', 'alt_final'])
        self.add_subsystem('con_cmp2', ExecComp('con2 = A_s_1  '), promotes=['con2', 'A_s_1'])
        self.add_subsystem('con_cmp3', ExecComp('con3 = nx_max  '), promotes=['con3', 'nx_max'])
        self.add_subsystem('con_cmp5', ExecComp('con5 = Pdyn_max  '), promotes=['con5', 'Pdyn_max'])
        
        self.add_design_var('pc1', lower=0.996, upper=1.126)
        self.add_design_var('TW1', lower=1.01, upper=2.0)
        self.add_design_var('f1', lower=0.08, upper=0.16)
        self.add_design_var('C1', lower=0.2, upper=0.3)
        self.add_design_var('OF1', lower=5.92, upper=6.05)
        self.add_design_var('t_a', lower=2.40, upper=3.10)
        self.add_design_var('t_vertical', lower=1.0, upper=4.0)        
        
        self.add_objective('obj')
        
        self.add_constraint('con1',lower = 80) # final altitude  
        #self.add_constraint('con2',lower = 16) # first stage area   
        self.add_constraint('con3',upper = 4.5 ) # maximum axial load factor
        self.add_constraint('con5',upper = 40 ) # maximum dynamic pressure allowed
        
        


prob = Problem()
prob.model = LAST_MDO()
prob.driver = ScipyOptimizeDriver()
prob.driver.options['optimizer'] = 'SLSQP'#'L-BFGS-B'#
prob.driver.options['maxiter'] = 200
prob.driver.options['tol'] = 1e-5
prob.driver.options['debug_print'] = ['desvars','ln_cons','nl_cons','objs','totals']
prob.driver.options['disp'] = True
prob.set_solver_print(level=1)
prob.setup()

# Ask OpenMDAO to finite-difference across the model to compute the gradients for the optimizer
prob.model.approx_totals()
prob.run_driver()

view_model(prob, outfile="mdao2.html", show_browser=False)


# generate as file for OpenVSP
generate_geometry =       False # True#

# view values calculations
view_values =         False #  True #

# plot the trajectory
plot_trajectory =            False # True #
MDA_value = False

# visualize output values
if plot_trajectory == True:
    #result_vizualization.view_traj_values(test)
    result_vizualization.plot_trajectory(MDA_value)
if view_values == True:
    result_vizualization.view_values(prob)

