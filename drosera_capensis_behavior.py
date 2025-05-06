#Python Script to model the behavior of the carnivorous Drosera Capensis Plant in the presence of food.
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os

'''

Behavior of the plant depends on the location of the food. Use the following parameters for different 
food locations and spreads:

1. Food spread out around ~0.75*L: gamma = 2, lambda_ = 10, mu = 0.75*L, sigma = 0.25*L
2. Food concentrated at ~0.75*L: gamma = 1, lambda_ = 35, mu = 0.75*L, sigma = 0.05*L
2. Food concentrated at the tip: gamma = 0.5, lambda_ = 120, mu = L, sigma = 0.05*L

'''

# Parameters
L = 1.0           # Plant length
N = 100           # Number of segments
dt = 0.01         # Time step
T = 3.0           # Total simulation time
steps = int(T / dt)

gamma = 2      # Proprioception strength 
lambda_ = 10     # Food sensitivity strength
s = np.linspace(0, L, N)
ds = s[1] - s[0]

# Initial curvature: small random perturbation
curvature = 0.01 * np.random.randn(N)

# Food distribution: smooth Gaussian centered at mu
mu = 0.75*L #Location of the center of the food; 
sigma = 0.25 * L #Standard deviation of the food; how spread out it is
f = np.exp(-((s - mu)**2) / (2 * sigma**2))

#gaussian distribution of food centered around the location of the food (mu)
def food(s, mu, sigma):
    return (1/(np.sqrt(2*np.pi)*sigma))*(np.exp(-((s - mu)**2) / (2 * sigma**2))) 

# Manually highlight the food region: 
highlight_mask = (s >= 0.97 * L) & (s <= 1.0 * L)

# History for animation
x_history = []
y_history = []
highlight_indices = []

def compute_theta(curvature, ds):
    return np.cumsum(curvature) * ds

#simple differential equation inspired from graviproprioception paper, using food as a gaussian
def compute_dCdt(curvature, f, gamma, lambda_):
    return lambda_ * f - gamma * curvature

#detect when the plant touches itself
def detect_repeat(arr, sensitivity):
    for index, point in enumerate(arr):
        for i in range(index + 4, len(arr)):
            if np.linalg.norm(point - arr[i]) < sensitivity:
                print(f"Match at {index}, {i}")  # optional printed statement to see what indices touch
                return True, index  # return the lower match index. 
    return False, None

#initializing constants needed for sensing if it touches itself
old_index = 99
index = 100        

# Time integration
for step in range(steps):
    dCdt = compute_dCdt(curvature, f, gamma, lambda_)
    curvature[:index] += dt * dCdt[:index]
    curvature[:N//3] = 0 #bottom third of the plant does not move, as recorded in experiment

    # Convert curvature to angle
    theta = compute_theta(curvature, ds)

    # Convert angle to (x, y)
    x = np.cumsum(np.sin(theta)) * ds
    y = np.cumsum(np.cos(theta)) * ds
    
    #detect when the plant touches itself
    points = np.array([x,y]).T
    intersect, index = detect_repeat(points,1.03e-02)

    #new, smaller curvature when the plant touches itself
    if intersect == True and index < old_index:
        print("new hit",index,", ",old_index)
        distance_from_food = np.linalg.norm(points[index]-points[int(mu*N-1)])
        if distance_from_food < 0.2:  #if the food is close to where it hits itself, it rolls 
            sigma += 0.05*L      #slightly increase the SD each time the plant hits itself 
            mu = float(index*L/N) #the effective new location of the food is where the plant hits itself
            f = food(s,mu, sigma) 
        else:                     #if the food isn't close to where it hits itself, it stops moving
            index = 10            #movement is now completely restricted  

    x_history.append(x.copy())
    y_history.append(y.copy())
    highlight_indices.append(np.where(highlight_mask)[0])

# Animation
fig, ax = plt.subplots(figsize=(5, 5))
line_full, = ax.plot([], [], lw=2, color='green')
line_highlight, = ax.plot([], [], lw=2, color='red')
ax.set_xlim(-L / 2, L / 2)
ax.set_ylim(0, L)
ax.set_title("Carnivorous Plant Under Food Influence")
ax.set_xlabel("x")
ax.set_ylabel("y")

def init():
    line_full.set_data([], [])
    line_highlight.set_data([], [])
    return line_full, line_highlight

def update(frame):
    x = x_history[frame]
    y = y_history[frame]
    line_full.set_data(x, y)
    highlight_idx = highlight_indices[frame]
    line_highlight.set_data(x[highlight_idx], y[highlight_idx])
    return line_full, line_highlight

ani = FuncAnimation(fig, update, frames=len(x_history), init_func=init, blit=True, interval=50)
#plt.show()

desktop_path = os.path.join(os.path.expanduser("~"), "Desktop", "plant_bending.mp4")
# Save animation (requires ffmpeg for mp4, or use Pillow for gif)
#ani.save("my_animation.gif", writer="pillow", fps=40)
#ani.save(desktop_path, writer='ffmpeg', fps=30)
plt.show()
