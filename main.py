from solver import solver
import misc
import numpy as np

########################################################################
#             Begin of parameters
title = 'discontinuity'
def u0(x):
    if x < 0:
        return 1.
    else:
        return 0.2

#title = 'Gaussian'
#def u0(x):
#    return misc.gaussian(x,-0.5,0.15)
    
tmax = 10.
h = 0.01
koverh = 0.8
range_x = [-1.,6.]

ulim = [-0.0,1.19]
xlim = [1.5,3.5]
#xlim = range_x
a = 1.0

#finite_diference = "BackwardEuler"
#finite_diference = "BackwardEulerNoCicle"
#finite_diference = "OneSidedUp"
#finite_diference = "OneSidedUpNoCicle"
#finite_diference = "OneSidedDown"
#finite_diference = "LaxFriedrichs"
#finite_diference = "LaxFriedrichsNoCicle"
#finite_diference = "BeamWarming"
#finite_diference = "BeamWarmingNoCicle"
#finite_diference = "LaxWendroff"
#finite_diference = "OneSidedUpBurgers"
finite_diference = "OneSidedUpBurgersNoCicle"

#             End of Parameters
########################################################################
k = h*koverh
step_t = int(tmax/k) + 1
step_x = int((range_x[1] - range_x[0])/h) + 1
range_t = [0.,tmax]
print h
print k
# Load the selected function to do the Finite Diference
if (finite_diference == "BackwardEuler"):
    from solver import BackwardEuler as FiniteMethod
elif (finite_diference == "BackwardEulerNoCicle"):
    from solver import BackwardEulerNoCicle as FiniteMethod
if (finite_diference == "OneSidedUp"):
    from solver import OneSidedUp as FiniteMethod
if (finite_diference == "OneSidedUpNoCicle"):
    from solver import OneSidedUpNoCicle as FiniteMethod
if (finite_diference == "OneSidedDown"):
    from solver import OneSidedDown as FiniteMethod
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
if (finite_diference == "OneSidedUpBurgers"):
    from solver import OneSidedUpBurgers as FiniteMethod
if (finite_diference == "OneSidedUpBurgersNoCicle"):
    from solver import OneSidedUpBurgersNoCicle as FiniteMethod

#Set-up the solver
s = solver(range_x,step_x,range_t,step_t,u0,FiniteMethod,a)
#Run the solver
s.solve()
np.save('/home/nano/DataDumpTCC/'+title+finite_diference,s.U)
#s.U = np.load('data.npy')
#s.U = np.load('/home/nano/DataDumpTCC/'+title+finite_diference+'.npy')

#Plot solution
misc.plot_ts(s,[5],xlim,ulim,finite_diference,'pdf','Burgers Finite , t=5 ,h='+str(h)+',k='+str(k),'','x','u')

#Animate solution
misc.animate(s,range_x,ulim,1,finite_diference+'.mp4',finite_diference+', '+title+', a='+str(a)+',h='+str(h)+',k='+str(k))
