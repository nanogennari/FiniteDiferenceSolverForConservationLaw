from solver import solver
import misc

########################################################################
#             Begin of parameters
def u0(x):
    return misc.gaussian(x,0.5,0.05)
    
step_x = 500
step_t = 1000

range_x = [0.,1.]
range_t = [0.,1.]

a = 1.

#finite_diference = "BackwardEuler"
finite_diference = "BeamWarming"
#finite_diference = "LaxWendroff"

#             End of Parameters
########################################################################

# Load the selected function to do the Finite Diference
if (finite_diference == "BackwardEuler"):
    from solver import BackwardEuler as FiniteMethod
elif (finite_diference == "BeamWarming"):
    from solver import BeamWarming as FiniteMethod
elif (finite_diference == "LaxWendroff"):
    from solver import LaxWendroff as FiniteMethod

#Set-up the solver
s = solver(range_x,step_x,range_t,step_t,u0,FiniteMethod,a)
#Run the solver
s.solve()
#Animate solution
misc.animate(s,finite_diference+'.mp4',finite_diference+', gaussian mu=0.5 sig=0.05')
