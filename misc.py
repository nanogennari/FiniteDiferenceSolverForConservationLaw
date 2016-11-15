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

def animate(s,filename='video.mp4',text=''):
    # Set up the figure, the axis, and the plot element for animation
    fig = plt.figure()
    ax = plt.axes(xlim=(0, 1), ylim=(0, 1))
    line, = ax.plot([], [], lw=2)
    plt.xlabel('x')
    plt.ylabel('u')
    plt.suptitle(text, fontsize=14, fontweight='bold')
    # Define function to update graph
    def animate(i):
        u,x,t = s.get_snap(i)
        line.set_data(x,u)
        plt.title('t='+str(t))
        return line,
    # Animate and save file
    anim = animation.FuncAnimation(fig, animate,frames=s.steps_t, interval=20, blit=True)
    anim.save(filename, fps=30, extra_args=['-vcodec', 'libx264'])

# Gaussian Function to calculate u0
def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2. * np.power(sig, 2.)))
