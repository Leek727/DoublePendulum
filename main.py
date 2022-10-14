import pygame
import sys
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import math

pygame.init()

screen = pygame.display.set_mode([500, 500])

# ----------------------- CONSTANTS AND CALCULATIONS -----------------------
# Pendulum rod lengths (m), bob masses (kg).
L1, L2 = 1, 1
m1, m2 = 1, 1
# The gravitational acceleration (m.s-2).
g = 9.81

def deriv(y, t, L1, L2, m1, m2):
    """Return the first derivatives of y = theta1, z1, theta2, z2."""
    theta1, z1, theta2, z2 = y

    c, s = np.cos(theta1-theta2), np.sin(theta1-theta2)

    theta1dot = z1
    z1dot = (m2*g*np.sin(theta2)*c - m2*s*(L1*z1**2*c + L2*z2**2) -
             (m1+m2)*g*np.sin(theta1)) / L1 / (m1 + m2*s**2)
    theta2dot = z2
    z2dot = ((m1+m2)*(L1*z1**2*s - g*np.sin(theta2) + g*np.sin(theta1)*c) + 
             m2*L2*z2**2*s*c) / L2 / (m1 + m2*s**2)
    return theta1dot, z1dot, theta2dot, z2dot

def calc_E(y):
    """Return the total energy of the system."""

    th1, th1d, th2, th2d = y.T
    V = -(m1+m2)*L1*g*np.cos(th1) - m2*L2*g*np.cos(th2)
    T = 0.5*m1*(L1*th1d)**2 + 0.5*m2*((L1*th1d)**2 + (L2*th2d)**2 +
            2*L1*L2*th1d*th2d*np.cos(th1-th2))
    return T + V

# ----------------------- PENDULUM CLASS DEFINITON -----------------------
class pendulum:
    def __init__(self, theta1, theta2):

        # Maximum time, time point spacings and the time grid (all in s).
        tmax, dt = 30, 0.001
        t = np.arange(0, tmax+dt, dt)
        # Initial conditions: theta1, dtheta1/dt, theta2, dtheta2/dt.
        y0 = np.array([math.radians(theta1), 0, math.radians(theta2), 0])

        # Do the numerical integration of the equations of motion
        y = odeint(deriv, y0, t, args=(L1, L2, m1, m2))

        # Unpack z and theta as a function of time
        theta1, theta2 = y[:,0], y[:,2]

        # Convert to Cartesian coordinates of the two bob positions.
        self.x1 = L1 * np.sin(theta1)
        self.y1 = -L1 * np.cos(theta1)
        self.x2 = self.x1 + L2 * np.sin(theta2)
        self.y2 = self.y1 - L2 * np.cos(theta2)


    def make_plot(self, screen, i):
        # rods
        pygame.draw.line(screen, (255,0,0), (0+250, 0+250), (self.x1[i]*100+250, -self.y1[i]*100+250), 3)
        pygame.draw.line(screen, (255,0,0), (self.x1[i]*100+250, -self.y1[i]*100+250), (self.x2[i]*100+250, -self.y2[i]*100+250), 3)
        #ax.plot([0, x1[i], x2[i]], [0, y1[i], y2[i]], lw=2, c='k')
        # bobs
        pygame.draw.circle(screen, (0,0,255), (0+250, 0+250), 5)
        pygame.draw.circle(screen, (0,0,255), (self.x1[i]*100+250, -self.y1[i]*100+250), 5)
        pygame.draw.circle(screen, (0,0,255), (self.x2[i]*100+250, -self.y2[i]*100+250), 5)

# Make an image every di time points, corresponding to a frame rate of fps
# frames per second.
# Frame rate, s-1
p = []
for x in range(1,100):
    p.append(pendulum(170+x/1000000,180+x/1000000))
    print(x)

i = 0
clock = pygame.time.Clock()
running = True
while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False
        
    screen.fill((0, 0, 0))
    #s = pygame.Surface((500,500), pygame.SRCALPHA)   # per-pixel alpha
    #s.fill((0,0,0,1))                         # notice the alpha value in the color
    #screen.blit(s, (0,0))
    for pend in p:
        pend.make_plot(screen, i)



    pygame.display.flip()
    clock.tick(1000)
    i += 1

pygame.quit()
