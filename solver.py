import numpy as np
from scipy.integrate import quad

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

# Conservative Gudunov
def GudunovRoe(Un,k,h,fl):
    m = len(Un)
    Unp = np.zeros([m],dtype=float)
    
    for i in range(0,m):
        Unp[i] =  Un[i] - (k/h)*(f(Un[i]) +)
    return Unp
