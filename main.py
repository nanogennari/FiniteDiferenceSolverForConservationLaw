from solver import solver
import misc
import numpy as np

########################################################################
#             Begin of parameters
title = 'riemman u_l=1, u_d=0'
def u0(x):
    if x < -0.5:
        return 1.
    else:
        return 0.

#title = 'gaussian mu=0.5 sig=0.05'
#def u0(x):
#    return misc.gaussian(x,0.5,0.05)
    
step_x = 500
step_t = 1000

range_x = [-1.,1.]
#range_x = [0.,1.]
range_t = [0.,1.]

ulim = [-0.2,1.2]

a = 1.0

#finite_diference = "BackwardEuler"
#finite_diference = "BackwardEulerNoCicle"
#finite_diference = "LaxFriedrichs"
finite_diference = "LaxFriedrichsNoCicle"
#finite_diference = "BeamWarming"
#finite_diference = "BeamWarmingNoCicle"
#finite_diference = "LaxWendroff"

#             End of Parameters
########################################################################

# Load the selected function to do the Finite Diference
if (finite_diference == "BackwardEuler"):
    from solver import BackwardEuler as FiniteMethod
elif (finite_diference == "BackwardEulerNoCicle"):
    from solver import BackwardEulerNoCicle as FiniteMethod
elif (finite_diference == "LaxFriedrichs"):
    from solver import LaxFriedrichs as FiniteMethod
elif (finite_diference == "LaxFriedrichsNoCicle"):
    from solver import LaxFriedrichsNoCicle as FiniteMethod
elif (finite_diference == "BeamWarming"):
    from solver import BeamWarming as FiniteMethod
elif (finite_diference == "BeamWarmingNoCicle"):
    from solver import BeamWarmingNoCicle as FiniteMethod
elif (finite_diference == "LaxWendroff"):
    from solver import LaxWendroff as FiniteMethod

#Set-up the solver
s = solver(range_x,step_x,range_t,step_t,u0,FiniteMethod,a)
#Run the solver
s.solve()
np.save('data',s.U)
#s.U = np.load('data.npy')
#Animate solution
misc.animate(s,range_x,ulim,finite_diference+'.mp4',finite_diference+', '+title)
