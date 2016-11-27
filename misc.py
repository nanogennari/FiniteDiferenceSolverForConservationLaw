import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

#Plot a graph at U[t_step)
def plot_snap(s,t_step):
    u,x,t = s.get_snap(t_step)
    plt.plot(x,u)
    plt.xlabel('x')
    plt.ylabel('u')
    plt.suptitle('t='+str(t), fontsize=14, fontweight='bold')
    plt.show()
    plt.close()

#Create a video from a solved system
def animate(s,xlim,ylim,skip='1',filename='video.mp4',text=''):
    # Set up the figure, the axis, and the plot element for animation
    fig = plt.figure()
    ax = plt.axes(xlim=(xlim[0], xlim[1]), ylim=(ylim[0], ylim[1]))
    line, = ax.plot([], [], lw=2)
    # Define function to update graph
    def animate(i):
        plt.clf()
        i = i*skip
        U,x,t = s.get_snap(i)
        u = np.zeros([len(x)],dtype=float)
        for i in range(0,len(x)):
            if x[i] < (0.6*t):
                u[i] = 1.
            else:
                u[i] = 0.2
        plt.ylim(ylim)
        plt.xlim(xlim)
        plt.xlabel('x')
        plt.ylabel('u')
        plt.suptitle(text, fontsize=14, fontweight='bold')
        plt.plot(x,U,'kx', label='Ujn')
        plt.plot(x,u, label='u(x,t)')
        plt.title('t='+str(round(t,3)))
        if (i%10 == 0):
            print "t="+str(t)
        return line,
    # Animate and save file
    anim = animation.FuncAnimation(fig, animate,frames=int(len(s.U)/skip), interval=20, blit=True)
    anim.save(filename, fps=30, extra_args=['-vcodec', 'libx264'])
    plt.close()

# Gaussian Function to calculate u0
def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2. * np.power(sig, 2.)))

def plot_ts(s,ts,xlim,ylim,filename,form,title,subtitle,xlabel='x',ylabel='u'):
    plt.figure(figsize=(8,6))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.ylim(ylim)
    plt.xlim(xlim)
    plt.suptitle(subtitle)
    for t in ts:
        i = s.get_t_step(t)
        U,x,t = s.get_snap(i)
        u = np.zeros([len(x)],dtype=float)
        for i in range(0,len(x)):
            if x[i] < (0.6*t):
                u[i] = 1.
            else:
                u[i] = 0.2
        plt.plot(x,U,'kx',label="Unj")
        plt.plot(x,u,label="u(x,t)")
    plt.title(title)
    plt.legend(ncol=3)
    plt.savefig(filename+'.'+form,format=form,dpi=200)
    
    plt.close()
