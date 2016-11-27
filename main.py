from solver import solver
import misc
import numpy as np

########################################################################
#             Begin of parameters
#title = 'discontinuity'
#def u0(x):
#    if x < 0:
#        return 1.
#    else:
#        return 0.2

#title = 'Gaussian'
#def u0(x):
#    return misc.gaussian(x,0,0.15)+0.2

title = 'sin'
def u0(x):
    if x > 0 and x < np.pi:
        return 5*np.sin(x)
    else:
        return 0

#def f(u,a):
#    return a[0]*u
#def fl(u,a):
#    return a[0]
def f(u,a):
    return (u*u)/2
def fl(u,a):
    return u

tmax = 5.
h = 0.04
koverh = 0.1
range_x = [0.,10.]
k = h*koverh
k = 0.004
step_t = int(tmax/k) + 1
step_x = int((range_x[1] - range_x[0])/h) + 1
range_t = [0.,tmax]

ulim = [0.,5.]
#xlim = [1.5,3.5]
a = 1.0
xlim = range_x
#fpar = a
bound = [0.,0.]
fpar = [f,fl,bound,a]

#method = "BackwardEuler"
#method = "BackwardEulerNoCicle"
#method = "OneSidedUp"
#method = "OneSidedUpNoCicle"
#method = "OneSidedDown"
#method = "LaxFriedrichs"
#method = "LaxFriedrichsNoCicle"
#method = "BeamWarming"
#method = "BeamWarmingNoCicle"
#method = "LaxWendroff"
#method = "OneSidedUpBurgers"
#method = "OneSidedUpBurgersNoCicle"
#method = "Gudunov"
method = "vanLeer"
#method = "vanLeer minmod"
#method = "vanLeer superbee"

#             End of Parameters
########################################################################

print h
print k
# Load the selected function to do the Finite Diference
if (method == "BackwardEuler"):
    from solver import BackwardEuler as Method
elif (method == "BackwardEulerNoCicle"):
    from solver import BackwardEulerNoCicle as Method
if (method == "OneSidedUp"):
    from solver import OneSidedUp as Method
if (method == "OneSidedUpNoCicle"):
    from solver import OneSidedUpNoCicle as Method
if (method == "OneSidedDown"):
    from solver import OneSidedDown as Method
elif (method == "LaxFriedrichs"):
    from solver import LaxFriedrichs as Method
elif (method == "LaxFriedrichsNoCicle"):
    from solver import LaxFriedrichsNoCicle as Method
elif (method == "BeamWarming"):
    from solver import BeamWarming as Method
elif (method == "BeamWarmingNoCicle"):
    from solver import BeamWarmingNoCicle as Method
elif (method == "LaxWendroff"):
    from solver import LaxWendroff as Method
elif (method == "OneSidedUpBurgers"):
    from solver import OneSidedUpBurgers as Method
elif (method == "OneSidedUpBurgersNoCicle"):
    from solver import OneSidedUpBurgersNoCicle as Method
elif (method == "Gudunov"):
    from solver import Gudunov as Method
elif (method == "vanLeer"):
    from solver import vanLeer as Method
elif (method == "vanLeer minmod"):
    from solver import vanLeerMinmod as Method
elif (method == "vanLeer superbee"):
    from solver import vanLeerSuperbee as Method


#Set-up the solver
s = solver(range_x,step_x,range_t,step_t,u0,Method,fpar)
#Run the solver
s.solve()
np.save('/home/nano/DataDumpTCC/'+title+method,s.U)
#s.U = np.load('data.npy')
#s.U = np.load('/home/nano/DataDumpTCC/'+title+method+'.npy')

#Plot solution
misc.plot_ts(s,[0,1,2,3,4,4.99],xlim,ulim,method,'pdf',title+', '+method+', h='+str(h)+',k='+str(k),'','x','u')

#Animate solution
misc.animate(s,range_x,ulim,10,method+'.mp4',title+', '+method+', a='+str(a)+',h='+str(h)+',k='+str(k))
