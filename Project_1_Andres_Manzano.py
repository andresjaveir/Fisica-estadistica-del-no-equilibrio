# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 23:50:10 2022

@author: Usuario
"""

import math
import numpy as np
import matplotlib.pyplot as plt


#%% 2)

## The goal is to generate 10^4 random numbers that follow a normal distribution.
## In order to do so, the Box-Muller algorithm will be used.

## The first step is to generate 2 random arrays, U1 and U2
## In order to do so, a random seed is set so that the numbers are the most random
#  possible.

np.random.seed()

U1 = np.random.random_sample((10000,))
U2 = np.random.random_sample((10000,))

## Then we set R and O.

R = np.sqrt(-2*np.log(U1))
O = np.pi*2*U2

## Now we can generate 2 arrays of 10^4 random numbers, X and Y, which
#  distribution follow a normal distribution.

X = R*np.cos(O)
Y = R*np.sin(O)

## We are only interested in one normal distribution, let's take X and plot it
#  to check that everything is correct.

plt.figure(figsize=(16,16))
plt.subplot(2,2,1)
plt.title(r"X", fontsize=22)
plt.hist(X, bins=50, color = 'red')
plt.grid('on', linestyle='--', linewidth=1, alpha=0.7)
plt.subplot(2,2,2)
plt.title(r"Y", fontsize=22)
plt.hist(Y, bins=50)
plt.grid('on', linestyle='--', linewidth=1, alpha=0.7)
## We now calculate the mean value and the variance

## Mean value (MV)

S = 0

for i in range(len(X)):
    
    S = S + X[i]

    if i == 10000:
        break
    
MV = S/len(X)

## Variance (V)

D = 0

for i in range(len(X)):
    
    D = D + (X[i]-MV)**2
    
    if i == 10000:
        break

V = D/len(X)


#%% 3)

## We want to plot the trajectory of a single particle that follows the  
#  Langevin equation in an open box. We start by defining the parameters
#  of the equation (the time and the time-jumps, since it's discrete).

dt = 0.1                  ## Time jump
T = 100*dt                  ## Total time
n = int(T/dt)                    ## Number of jumps
t = np.linspace(0, T, 100)  ## Time vector

x = np.zeros(n)             ## x coordinate
y = np.zeros(n)             ## y coordinate

g = 1                    ## Amplitude of the noise term
for i in range(n-1):
    
    x[i+1] = x[i] + np.sqrt(2*g*dt)*X[i]
    y[i+1] = y[i] + np.sqrt(2*g*dt)*Y[i]
    
    if i == n-1:
        break

    
plt.plot(x, y, '--ko')
plt.title(r"Brownian dynamics in an open box")
plt.ylabel(r"y", labelpad=1, fontsize=12)
plt.xlabel(r'x', labelpad=1, fontsize=12)
plt.grid('on', linestyle='--', linewidth=1, alpha=0.7)


#%% 4)

## Same as before, but with periodic boundary conditions and 1000 particles
#  which start at random positions

N = 1000                ## Number of particles


xo = np.zeros(N)        ## Initial x value for the particles
yo = np.zeros(N)        ## Initial y value for the particles

## We need to create a normal distribution of at least 100000 random
#  numbers, to simulate the movement of each particle, just like we did
#  in exercise number 2

np.random.seed()

Z1 = np.random.random_sample((100000,))
Z2 = np.random.random_sample((100000,))

Q = np.sqrt(-2*np.log(Z1))
Theta = np.pi*2*Z2

X4 = Q*np.cos(Theta)
Y4 = Q*np.sin(Theta)



## Now we split the distribution in 1000 different ones of 100 values, one 
#  for each particle, and 100 valuse because the particle is going to move 
#  100 times

def split_list(alist, wanted_parts=1):
    length = len(alist)
    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts] 
             for i in range(wanted_parts) ]

K = split_list(X4, 1000)
L = split_list(Y4, 1000)
G = 1     ## New value for the amplitude of the noise term




## Now we simulate the movement of each particle, just like in the previous
#  exercise, but now with 100000 particles

for k in range(1000):
    for i in range(100):
        
        xo[k] = xo[k] + np.sqrt(2*G*dt)*K[k][i]
        yo[k] = yo[k] + np.sqrt(2*G*dt)*L[k][i]
        
            
    if (k == 999) and (i == 99):
        break        
    
    
    

## The initial positions of the particles, chosen from a random uniform 
#  distribution that goes from -50 to 50

xin = (-1+ np.random.random_sample((1000))*2)*50
yin = (-1+ np.random.random_sample((1000))*2)*50




## The final positions

xdef = xo + xin
ydef = yo + yin



## Now we set the boundary conditions

for k in range(1000):

    if(xdef[k]) > 50:
            
        xdef[k] = xdef[k] - 100
            
    if(xdef[k]) < -50: 
            
        xdef[k] = xdef[k] + 100
            
    if(ydef[k]) > 50:
            
        ydef[k] = ydef[k] - 100
            
    if(ydef[k]) < -50: 
            
        ydef[k] = ydef[k] + 100
            
## We plot the results

plt.plot(xin, yin, 'ko', markersize = 2, label='Initial position')
plt.plot(xdef, ydef, 'ro', markersize = 2, label='Final position')
plt.legend(loc="best")
plt.grid('on', linestyle='--', linewidth=1, alpha=0.7)
plt.ylabel("Y", fontsize=10)
plt.xlabel("X", fontsize=10)
plt.legend(bbox_to_anchor = (0.75, 1.15), ncol = 2)
plt.show()

#%% 5)

## We create a loop to graph the difference between the final and the initial 
#  positions for different values of Gamma.

## We restart the positions

xin = (-1+ np.random.random_sample((1000))*2)*50
yin = (-1+ np.random.random_sample((1000))*2)*50

xo = np.zeros(N)        
yo = np.zeros(N)

G = [0.1,1,10,30,50,100]
Delta = np.zeros(len(G))

for s in range(len(G)):
    for k in range(1000):
        for i in range(100):
        
            xo[k] = xo[k] + np.sqrt(2*G[s]*dt)*K[k][i]
            yo[k] = yo[k] + np.sqrt(2*G[s]*dt)*L[k][i]

    xdef = xo + xin
    ydef = yo + yin
    
    Delta[s] = (1/N)*np.sum((xdef-xin)**2+(ydef-yin)**2)/(4*T)  
    xo = np.zeros(N)
    yo = np.zeros(N)        
    if (k == 999) and (i == 99) and (s==len(G)):
        break       
        
        
plt.plot(G, Delta, '--ro')
plt.xscale(value="log")
plt.yscale(value="log")
plt.grid('on', linestyle='--', linewidth=1, alpha=0.7)
plt.ylabel(r"D", labelpad=1, fontsize=12)
plt.xlabel(r'Γ', labelpad=1, fontsize=12)
plt.title(r"D vs Γ")

#%% 6)

## We start by setting 5 different sets of N=1000 particles at the origin

x6 = np.zeros(1000)
y6 = np.zeros(1000)
z6 = np.zeros(1000)
h6 = np.zeros(1000)
f6 = np.zeros(1000)

## We then create a normal distribution with all the random numbers we will
#  need

M1 = np.random.random_sample((10000000,))
M2 = np.random.random_sample((10000000,))

M = np.sqrt(-2*np.log(M1))
B = np.pi*2*M2

X6 = M*np.cos(B)
Y6 = M*np.sin(B)

A = split_list(X6, 1000)
Z = split_list(Y6, 1000)

dt = 0.1        ## Value of the time jump

G6 = 1          ## Value of the amplitude of the noise

## Now we simulate the movement of the 5 sets of particle, to each set
#  corresponds a different amount of iterations, which is equal to a different
#  number of iterations, since each iteration is equal to going forward in 
#  time t = 0.1


for k in range(1000):
    for i in range(10):
        
        x6[k] = x6[k] + np.sqrt(2*G6*dt)*A[k][i]
        
            
    if (k == 999) and (i == 9):
        break        
    
for k in range(1000):
    for i in range(100):
        
        y6[k] = y6[k] + np.sqrt(2*G6*dt)*A[k][i]
        
            
    if (k == 999) and (i == 99):
        break           
    
for k in range(1000):
    for i in range(300):
        
        z6[k] = z6[k] + np.sqrt(2*G6*dt)*A[k][i]
        
            
    if (k == 999) and (i == 299):
        break        

for k in range(1000):
    for i in range(1000):
        
        h6[k] = h6[k] + np.sqrt(2*G6*dt)*A[k][i]
        
            
    if (k == 999) and (i == 999):
        break  
    
    
for k in range(1000):
    for i in range(10000):
        
        f6[k] = f6[k] + np.sqrt(2*G6*dt)*A[k][i]
        
            
    if (k == 999) and (i == 9999):
        break              


## We plot the results

plt.figure(figsize=(8,6))
plt.hist(x6, bins=20, label = "t=1")
plt.hist(y6, bins=40, label = "t=10")
plt.hist(z6, bins=60, label = "t=30")
plt.hist(h6, bins=80, label = "t=100")
plt.hist(f6, bins=150, label = "t=1000")
plt.grid('on', linestyle='--', linewidth=1, alpha=0.7)
plt.legend(loc="upper right")
plt.title("Diffusion of particles for Γ=1")






















