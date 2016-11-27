import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

class solver:
    
    def __init__(self,range_x,steps_x,range_t,steps_t,u0,fdisc,fargs):
        #Store everything we need
        self.n = 0
        self.h = (range_x[1] - range_x[0])/steps_x
        self.x0 = range_x[0]
        self.k = (range_t[1] - range_t[0])/steps_t
        self.u0 = u0
        self.fdisc = fdisc
        self.fargs = fargs
        # Create array U
        self.U = np.zeros([steps_t,steps_x],dtype=float)
        #Initialize U[0] with cell mean values
        for i in range(0,len(self.U[0])):
            x0 = self.x0 + i*self.h
            x = self.x0 + (i + 1.)*self.h
            self.U[0][i] = quad(u0,x0,x)[0]/self.h
    
    def step(self):
        n = self.n + 1
        # Call finite diference method
        self.U[n] = self.fdisc(self.U[n-1],self.k,self.h,self.fargs)
        # Increment step
        self.n = n
        
    def get_t_step(self,t):
        return int(t/self.k)
    
    def solve(self,rep=10):
        for i in range(1,len(self.U)):
            self.step()
            if (i%rep == 0):
                print "  t = "+str(i*self.k)
                
    # Method used for printing data of time step n, returns U[n], a list of Xs and t at step n
    def get_snap(self,t_step):   
        x = []
        for i in range(0,len(self.U[0])):
            x.append(self.x0 + i*self.h)     
        return self.U[t_step],x,t_step*self.k

# Finite diference using BackwardEuler
def BackwardEuler(Un,k,h,a):
    m = len(Un)
    A = np.zeros([m,m],dtype=float)
    for i in range(0,m):
        A[(i-1)%m][i] = (k*a)/(2*h)
        A[i][i] = 1
        A[(i+1)%m][i] = -(k*a)/(2*h)
    B = Un
    return np.linalg.solve(A,B)
    
# Finite diference using BackwardEuler without cicle boundaries
def BackwardEulerNoCicle(Un,k,h,a):
    m = len(Un)
    A = np.zeros([m,m],dtype=float)
    for i in range(0,m):
        if i > 1:
            A[i-1][i] = (k*a)/(2*h)
        A[i][i] = 1
        if i < m-1:
            A[i+1][i] = -(k*a)/(2*h)
    B = Un
    return np.linalg.solve(A,B)

# Finite diference using OneSided Upwind
def OneSidedUp(Un,k,h,a):
    m = len(Un)
    Unp = np.zeros([m],dtype=float)
    C = (a*k)/(h)
    for i in range(0,m):
        Unp[i] =  Un[i] - C*(Un[i]-Un[(i-1)%m])
    return Unp

# Finite diference using OneSided Upwind
def OneSidedUpNoCicle(Un,k,h,a):
    m = len(Un)
    Unp = np.zeros([m],dtype=float)
    C = (a*k)/(h)
    Unp[0] = 1
    Unp[m-1] = 0
    for i in range(1,m-1):
        Unp[i] =  Un[i] - C*(Un[i]-Un[(i-1)%m])
    return Unp
    
# Finite diference using OneSided Downwind
def OneSidedDown(Un,k,h,a):
    m = len(Un)
    Unp = np.zeros([m],dtype=float)
    C = (a*k)/(h)
    for i in range(0,m):
        Unp[i] =  Un[i] - C*(Un[(i+1)%m]-Un[i])
    return Unp

# Finite diference using LaxFriedrichs
def LaxFriedrichs(Un,k,h,a):
    m = len(Un)
    Unp = np.zeros([m],dtype=float)
    C = (a*k)/(h)
    for i in range(0,m):
        Unp[i] =  ( (Un[(i-1)%m]+Un[(i+1)%m]) - C*(Un[(i+1)%m]-Un[(i-1)%m]))/2.
    return Unp
    
# Finite diference using LaxFriedrichs without cicle boundaries
def LaxFriedrichsNoCicle(Un,k,h,a):
    m = len(Un)
    Unp = np.zeros([m],dtype=float)
    C = (a*k)/(h)
    for i in range(0,m):
        if i == 0:
            Unp[i] =  ( (1+Un[i+1]) - C*(Un[i+1]-1))/2.
        elif i == m-1:
            Unp[i] =  ( Un[i-1] - C*(-Un[i-1]))/2.
        else:
            Unp[i] =  ( (Un[i-1]+Un[i+1]) - C*(Un[i+1]-Un[i-1]))/2.
    return Unp

# Finite diference using BeamWarming
def BeamWarming(Un,k,h,a):
    m = len(Un)
    Unp = np.zeros([m],dtype=float)
    C = (a*k)/(h)
    for i in range(0,m):
        Unp = Un[i] + ( (-1.)*C*(3.*Un[i]-4.*Un[(i-1)%m]+Un[(i-2)%m]) + (C*C)*(Un[i]-2.*Un[(i-1)%m]+Un[(i-2)%m]))/2.
    return Unp

# Finite diference using BeamWarming without cicle boundaries
def BeamWarmingNoCicle(Un,k,h,a):
    m = len(Un)
    Unp = np.zeros([m],dtype=float)
    C = (a*k)/(h)
    for i in range(0,m):
        if i > 1:
            Unp[i] = Un[i] + ( (-1.)*C*(3.*Un[i]-4.*Un[(i-1)%m]+Un[(i-2)%m]) + (C*C)*(Un[i]-2.*Un[(i-1)%m]+Un[(i-2)%m]))/2.
        elif i == 1:
            Unp[i] = Un[i] + ( (-1.)*C*(3.*Un[i]-4.*Un[(i-1)%m]+1) + (C*C)*(Un[i]-2.*Un[(i-1)%m]+1))/2.
        else:
            Unp[i] = Un[i] + ( (-1.)*C*(3.*Un[i]-4.*1+1) + (C*C)*(Un[i]-2.*1+1))/2.
    return Unp

# Finite diference using LaxWendroff
def LaxWendroff(Un,k,h,a):
    m = len(Un)
    Unp = np.zeros([m],dtype=float)
    C = (a*k)/(h)
    for i in range(0,m):
        Unp = Un[i] + ( (-1.)*C*(Un[(i+1)%m]-Un[(i-1)%m]) + (C*C)*(Un[(i+1)%m]-2.*Un[i]+Un[(i-1)%m]))/2.
    return Unp

# Finite diference using OneSided Upwind for Burgers Equation
def OneSidedUpBurgers(Un,k,h,a):
    m = len(Un)
    Unp = np.zeros([m],dtype=float)
    C = k/h
    for i in range(0,m):
        Unp[i] =  Un[i] - Un[i]*C*(Un[i]-Un[(i-1)%m])
    return Unp
    
# Finite diference using OneSided Upwind for Burgers Equation without cicle boundaries
def OneSidedUpBurgersNoCicle(Un,k,h,a):
    m = len(Un)
    Unp = np.zeros([m],dtype=float)
    Unp[0] = 1
    Unp[m-1] = 0.2
    for i in range(1,m-1):
        Unp[i] =  Un[i] - Un[i]*(k/h)*(Un[i]-Un[i-1])
    return Unp

# An rimman solver for Gudunov and vanLeer
def riemman_solver(ul,ur,f,fl,fpar):
    flul = fl(ul,fpar)
    flur = fl(ur,fpar)
    if flul >= 0:
        if flur >= 0:
            u = ul
        else:      
            foveru = (f(ur,fpar) - f(ul,fpar))/(ur - ul)
            if foveru >= 0:
                u = ul
            else:
                u = ur        
    else:
        if flur > 0:
            if Un[(i-1)%m] <= Un[i]:
                res = minimize_scalar(f,bounds=(ul,ur),args=(fpar),method='bounded')
                u = res.x
            else:
                res = minimize_scalar(-f,bounds=(ul,ur),args=(fpar),method='bounded')
                u = res.x  
        else:
            u = ur
    return u

# Conservative Gudunov
def Gudunov(Un,k,h,par):
    m = len(Un)
    bounds = par[2]
    f = par[0]
    fl = par[1]
    fpar = par[3:]
    
    Unp = np.zeros([m],dtype=float)
    if bounds != "cicle":
        Unp[0] = bounds[0]
        Unp[m-1] = bounds[1]
        steps = range(1,m-1)
        ur = riemman_solver(Un[0],Un[1],f,fl,fpar)
    else:
        steps = range(0,m)
        ur = riemman_solver(Un[m-1],Un[0],f,fl,fpar)
    
    for i in steps:
        ul = ur
        ur = riemman_solver(Un[i],Un[(i+1)%m],f,fl,fpar)
        ful = f(ul,fpar)
        fur = f(ur,fpar)
        Unp[i] =  Un[i] - (k/h)*(fur - ful)
    return Unp

# S functions for vanLeer
def vanLeerS(Un,i):
    m = len(Un)
    if (Un[(i+1)%m]-Un[i]) != 0:
        theta = (Un[i]-Un[(i-1)%m])/(Un[(i+1)%m]-Un[i])
        print theta
        if (theta > 0.):
            phi = (abs(theta) - theta)/(1+abs(theta))
            return (Un[(i+1)%m]-Un[i])*phi
        else:
            return 0
    else:
        return 0

def minmod(Un,i):
    m = len(Un)
    if (Un[(i+1)%m]-Un[i]) != 0:
        theta = (Un[i]-Un[(i-1)%m])/(Un[(i+1)%m]-Un[i])
        if theta > 0:
            phi = min(theta,1.)
            return (Un[(i+1)%m]-Un[i])*phi
        else:
            return 0
    else:
        return 0

def superbee(Un,i):
    m = len(Un)
    if (Un[(i+1)%m]-Un[i]) != 0:
        theta = (Un[i]-Un[(i-1)%m])/(Un[(i+1)%m]-Un[i])
        if theta > 0:
            phi = max(min(2*theta,1),min(theta,2))
            return (Un[(i+1)%m]-Un[i])*phi
        else:
            return 0
    else:
        return 0

# Conservative vanLeer
def vanLeerFunc(Un,k,h,par,S):
    m = len(Un)
    bounds = par[2]
    f = par[0]
    fl = par[1]
    fpar = par[3:]
    
    Unp = np.zeros([m],dtype=float)
    if bounds != "cicle":
        Unp[0] = bounds[0]
        Unp[m-1] = bounds[1]
        steps = range(2,m-2)
        ulr = 0
        urr = 0
    else:
        steps = range(0,m)
        
    for i in steps:
        ull = ulr
        ulr = urr
        url = Un[i] + .5*S(Un,i)
        urr = Un[i+1] - .5*S(Un,i+1)
        ul = riemman_solver(ull,ulr,f,fl,fpar)
        ur = riemman_solver(url,urr,f,fl,fpar)
        ful = f(ul,fpar)
        fur = f(ur,fpar)
        Unp[i] =  Un[i] - (k/h)*(fur - ful)
    return Unp

def vanLeer(Un,k,h,par):
    return vanLeerFunc(Un,k,h,par,vanLeerS)
    
def vanLeerMinmod(Un,k,h,par):
    return vanLeerFunc(Un,k,h,par,minmod)

def vanLeerSuperbee(Un,k,h,par):
    return vanLeerFunc(Un,k,h,par,superbee)

