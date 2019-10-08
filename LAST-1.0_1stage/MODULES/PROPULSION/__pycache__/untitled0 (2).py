""" Definition of the Paraboloid component, which evaluates the equation
(x-3)^2 + xy + (y+4)^2 = 3
"""
from __future__ import division, print_function
from openmdao.core.explicitcomponent import ExplicitComponent
import scipy as sp
from openmdao.api import Problem, Group, ImplicitComponent, IndepVarComp, NewtonSolver, ScipyKrylov, ExecComp,NonlinearBlockGS, LinearBlockGS, DirectSolver


class fuel(ExplicitComponent):
    """
    Evaluates the equation L1 = M1*(1-exp(-v1/w1))
    """

    def setup(self):
        self.add_input('M1', val=1.0)#Total mass of rocket
        self.add_input('v1', val=1.0)#Delta v
        self.add_input('w1', val=1.0)#Exhaust velocity

        self.add_output('L1', val=1.0)

        # Finite difference all partials.
        self.declare_partials('*', '*', method='fd')

    def compute(self, inputs, outputs):
        """
        L1 = M1*(1-exp(-v1/w1))

        """
        M1 = inputs['M1']
        v1 = inputs['v1']
        w1 = inputs['w1']

        outputs['L1'] = M1*(1-sp.exp(-v1/w1))
        
        
class struct(ExplicitComponent):
    """
    Evaluates the equation S1 = beta1*M1
    """

    def setup(self):
        self.add_input('M1', val=1.0)#Total mass of rocket
        self.add_input('beta1', val=0.2)#Structure ratio

        self.add_output('S1', val=1.0)

        # Finite difference all partials.
        self.declare_partials('*', '*', method='fd')

    def compute(self, inputs, outputs):
        """
        S1 = beta1*M1

        """
        M1 = inputs['M1']
        beta1 = inputs['beta1']

        outputs['S1'] = beta1*M1
              
             
class tank(ImplicitComponent):
    """
    Evaluates the equation T1 = alpha1*(T1+L1)
    """

    def setup(self):
        self.add_input('L1', val=1.0)#Amount of fuel
        self.add_input('alpha1', val=0.1)#Tank ratio

        self.add_output('T1', val=1.0)

        # Finite difference all partials.
        self.declare_partials(of='*', wrt='*')

    def apply_nonlinear(self, inputs, outputs,residuals):
        """
        T1 = alpha1*(T1+L1)

        """
        L1 = inputs['L1']
        alpha1 = inputs['alpha1']
        T1 = outputs['T1']
        
               
        residuals['T1'] = T1-alpha1*(T1+L1)
'''
class tank(ExplicitComponent):
    """
    Evaluates the equation S1 = beta1*M1
    """

    def setup(self):
        self.add_input('L1', val=1.0)#Total mass of rocket
        self.add_input('alpha1', val=0.05)#Structure ratio

        self.add_output('T1', val=1.0)

        # Finite difference all partials.
        self.declare_partials('*', '*', method='fd')

    def compute(self, inputs, outputs):
        """
        T1 = alpha1*L1/(1-alpha1)

        """
        L1 = inputs['L1']
        alpha1 = inputs['alpha1']

        outputs['T1'] = alpha1*L1/(1-alpha1)
'''




class totalmass(ExplicitComponent):
    """
    Evaluates the equation M1 = P+S1+T1+L1
    """

    def setup(self):
        self.add_input('P', val=1.0)
        self.add_input('S1', val=1.0)
        self.add_input('T1', val=1.0)
        self.add_input('L1', val=1.0)

        self.add_output('M1', val=1.0)

        # Finite difference all partials.
        self.declare_partials('*', '*', method='fd')

    def compute(self, inputs, outputs):
        """
        M1 = P+S1+T1+L1

        """
        P = inputs['P']
        S1 = inputs['S1']
        T1 = inputs['T1']
        L1 = inputs['L1']

        outputs['M1'] = P+S1+T1+L1


if __name__ == "__main__":



    probl = Problem()
    
    probl.model.add_subsystem('P',IndepVarComp('P', 1.0),promotes=['P'])
    probl.model.add_subsystem('alpha1',IndepVarComp('alpha1', 0.05),promotes=['alpha1'])
    probl.model.add_subsystem('beta1',IndepVarComp('beta1', 0.2),promotes=['beta1'])
    probl.model.add_subsystem('v1',IndepVarComp('v1', 1),promotes=['v1'])
    probl.model.add_subsystem('w1',IndepVarComp('w1', 1),promotes=['w1'])
    
    probl.model.add_subsystem('fuel',fuel(),promotes=['L1','M1','v1','w1'])
    probl.model.add_subsystem('struct',struct(),promotes=['S1','beta1','M1'])
    probl.model.add_subsystem('tank',tank(),promotes=['alpha1','T1','L1'])
    probl.model.add_subsystem('totalmass',totalmass(),promotes=['M1','P','S1','T1','L1'])
    
    #probl.model.linear_solver = LinearBlockGS()
    #probl.model.linear_solver.options['maxiter'] = 10
    probl.model.nonlinear_solver = NewtonSolver()
    #probl.model.linear_solver = DirectSolver()
    probl.model.nonlinear_solver.options['maxiter'] = 10**6
    probl.model.nonlinear_solver.options['rtol'] = 10**-3
    #probl.model.nonlinear_solver.options['max_sub_solves'] = 1000
    probl.set_solver_print(level=2)
    
    probl.setup()
    
    probl.run_model()
    
    print('P =',probl['P'])
    print('L1 =',probl['L1'])
    print('T1 =',probl['T1'])
    print('S1 =',probl['S1'])
    print('M1 =',probl['M1'])
    print('v1 =',probl['v1'])
    print('w1 =',probl['w1'])

#    model = Group()
#    ivc = IndepVarComp()
#    ivc.add_output('M1', 1.0)
#    ivc.add_output('v1', 1.0)
#    ivc.add_output('w1', 1.0)
#    model.add_subsystem('des_vars', ivc)
#    model.add_subsystem('fuel_comp', fuel())
#
#    model.connect('des_vars.M1', 'fuel_comp.M1')
#    model.connect('des_vars.v1', 'fuel_comp.v1')
#    model.connect('des_vars.w1', 'fuel_comp.w1')
#
#    prob = Problem(model)
#    prob.setup()
#    prob.run_model()
#    print(prob['fuel_comp.L1'])




#    prob['des_vars.x'] = 5.0
#    prob['des_vars.y'] = -2.0
#    prob.run_model()
#    print(prob['parab_comp.f_xy'])