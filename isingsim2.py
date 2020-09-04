from random import random,choice
from math import exp

n = 20                  # Sites per edge for n x n system
n2 = n*n                # Precalculate n*n
nlist = list(range(n))  # List of sites 
ntrials = 200000        # Number Trials 
nequil = 100000         # Equilibration steps
Tmax = 4.0              #max temperature
Temp = 0.0              #initial temperature
dT = 0.1                #change in temperature

while Temp < Tmax:
    Temp += dT

    # Initialize sums for averages
    E_sum = 0
    E2_sum = 0

    # Create initial matrix of spins all spin up
    spins=[[1 for i in range(n)] for j in range(n)]

    # Run simulation

    for trial in range(1, (ntrials + nequil + 1)):
        #randomly pick a spin state
        i = choice(nlist)
        j = choice(nlist)

        #Calculate the change in energy if we flip this spin
        deltaE=2.*(spins[i][j]*\
           (spins[i][(j+1)%n]+spins[i][(j-1+n)%n]+\
             spins[(i+1)%n][j]+spins[(i+-1+n)%n][j]))

        #Flip the spin using Metropolis MC
        if exp(-deltaE / Temp) > random():
            spins[i][j] = - spins[i][j]
        else:
            deltaE = 0.0

        # Calculate system energy ONCE 
        if trial == nequil: #if equilibration trials are complete
            energy = 0.0  
            for i in range(n):
                for j in range(n):
                    energy -= (spins[i][j]*
                               (spins[i][(j+1)%n]+spins[i][(j-1+n)%n]+\
                               spins[(i+1)%n][j]+spins[(i+-1+n)%n][j]))
            energy /= n2 #calculate average energy 
                    
        # Update energy based on deltaE per spin
        if trial > nequil: #determine if we are in production portion of runs
            energy += 2 * deltaE / n2
            E_sum += energy
            E2_sum += energy*energy

    #calculate averages and heat capacity
    E_ave = E_sum / ntrials
    E2_ave = E2_sum / ntrials
    Cv = 1 / (Temp * Temp) * (E2_ave - (E_ave * E_ave))

    print('%8.4f %10.6f %10.6f'%(Temp, E_ave, Cv))
