import numpy as np
from scipy.integrate import quad

class solver:
    
    def __init__(self,range_x,steps_x,range_t,steps_t,u0,fcoef,fargs):
        #Store everything we need
        self.n = 0
        self.h = (range_x[1] - range_x[0])/steps_x
        self.k = (range_t[1] - range_t[0])/steps_t
        self.u0 = u0
        self.fcoef = fcoef
        self.fargs = fargs
        # Create array U
        self.U = np.zeros([steps_t,steps_x],dtype=float)
        #Initialize U[0] with cell mean values
        for i in range(0,self.steps_x):
            x0 = i*self.h
            x = (i + 1.)*self.h
            self.U[0][i] = quad(u0,x0,x)[0]/self.h
    
    def step(self):
        n = self.n + 1
        # Call finite diference method to calculate arrays A and B for the linear solver
        A,B = self.fcoef(self.U[n-1],self.k,self.h,self.fargs)
        # Call the linear solver
        self.U[n] = np.linalg.solve(A,B)
        # Increment step
        self.n = n
    
    def solve(self,rep=10):
        for i in range(1,self.steps_t):
            self.step()
            if (i%rep == 0):
                print "  t = "+str(i*self.k)
                
    # Method used for printing data of time step n, returns U[n], a list of Xs and t at step n
    def get_snap(self,t_step):   
        x = []
        for i in range(0,self.steps_x):
            x.append(i*self.h)     
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
    return A,B

# Finite diference using BeamWarming
def BeamWarming(Un,k,h,a):
    m = len(Un)
    A = np.zeros([m,m],dtype=float)
    B = np.zeros([m],dtype=float)
    C = (a*k)/(h)
    for i in range(0,m):
        A[i][i] = 1
        B[i] = Un[i] + ( (-1.)*C*(3.*Un[i]-4.*Un[(i-1)%m]+Un[(i-2)%m]) + (C*C)*(Un[i]-2.*Un[(i-1)%m]+Un[(i-2)%m]))/2.
    return A,B

# Finite diference using LaxWendroff
def LaxWendroff(Un,k,h,a):
    m = len(Un)
    A = np.zeros([m,m],dtype=float)
    B = np.zeros([m],dtype=float)
    C = (a*k)/(h)
    for i in range(0,m):
        A[i][i] = 1
        B[i] = Un[i] + ( (-1.)*C*(Un[(i+1)%m]-Un[(i-1)%m]) + (C*C)*(Un[(i+1)%m]-2.*Un[i]+Un[(i-1)%m]))/2.
    return A,B
