"""
Daisyworld
Markus Roland Ernst, April 29th 2020

Program to calculate the percentage of "black" and "white" daisies covering a 
planet that has a surface comprised of only black and white daisies and bare 
ground, as well as the temperatures of the daisies and ground. Originally 
written in Fortran by Kirsten Menking and Joel Dashnaw, Dept. of Geology and 
Geography, Vassar College, Poughkeepsie, NY  12604, July 21, 2003.

Rewritten for the Python Language by Markus Ernst, Institute of Atmospheric
And Environmental Sciences, Altenhoefer Allee 1, Frankfurt am Main.
Based on the following readings: Watson, A.J., and Lovelock, J.E., 1983,
Biological homeostasis of the global environment: The parable of Daisyworld,
Tellus, v. 35B, p. 284-289.
"""

# _____________________________________________________________________________


# -----------------
# import libraries
# -----------------

# standard libraries
# -----

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
%matplotlib inline
import time
from io import BytesIO


# helper functions
# -----

# helper functions
# -----

def plot_world_trajectory(tmax, dt, daisydict, X, planet_temp):
    palette = sns.color_palette("Blues",10)
    tr = np.arange(0, tmax+dt, dt)
    #print(tr.shape)
    plt.subplot(2, 1, 1)
    for i, daisyname in enumerate(daisydict):
        plt.plot(tr, daisydict[daisyname], '-', color=palette[i+3+i*3], label=daisyname)
    plt.plot(tr, X, 'g--', label='fertile')
    plt.title('Daisyworld Output')
    plt.ylabel('% of area\ncovered with daisies')
    plt.legend()


    plt.subplot(2, 1, 2)
    plt.plot(tr, planet_temp, 'r-')
    plt.xlabel('time(yr)')
    plt.ylabel('temperature\nof the planet')
    
    #plt.legend()
    plt.show()


# constants
# -----

SOLCONST = 917.0      # J/(m^2*s)
STEFBOLTZ = 5.67e-8   # J/(m^2*s*degK^4)

# -----------------
# world classes
# -----------------

# this is your basic class for your daisyworld
class World(object):
    def __init__(self, name, fert=1.0, gralbedo=0.5, luminosity=0.5, lumacc=0.02, solconst=SOLCONST):
        self.name = name
        self.fert = fert
        self.gralbedo = gralbedo
        self.luminosity = luminosity
        self.lumacc = lumacc
        self.solconst = solconst
        
        self.time = 0.
        self.daisies = []
    
    def add_daisy(self, Daisy):
        self.daisies.append(Daisy)        
        area = 0
        for daisy in self.daisies:
            area += daisy.area
        if area > self.fert:
            self.daisies.pop()
            raise ValueError("daisies cannot occupy more than 100% of fertile area")
        else:
            pass
    
    def reset(self):
        self.time = [self.time[0]]
        self.luminosity = [self.luminosity[0]]
        for daisy in self.daisies:
            daisy.reset()
        pass
        
    def plot(self):
        pass
    
    def simulate_timestep(self, dt=0.2):
        self.X = self.fert
        for daisy in self.daisies:
            self.X -= daisy.area
        
        # planetary temperature
        # -----
        
        alb_planet = 0.
        remaining_frac = 1.
        for daisy in self.daisies:
            alb_planet += daisy.area * daisy.albedo
            remaining_frac -= daisy.area
        alb_planet += remaining_frac * self.gralbedo
        
        
        #self.luminosity = np.sin(np.pi/75 * self.time - np.pi) * 0.4 + 1.25    
        self.luminosity += self.lumacc * dt                 # unitless
        solflux = self.solconst * self.luminosity           # J/(m^2*s)
        q = (0.2 * solflux) / STEFBOLTZ                     # unitless
        
        self.planet_temp = (((solflux * (1 - alb_planet)) / STEFBOLTZ)**0.25) - 273.15
        
        # daisy temperature and growth
        # -----
        
        for daisy in self.daisies:
            daisy.temperature = (((q * (alb_planet - daisy.albedo)) + ((self.planet_temp + 273.16)**4))**(0.25)) - 273.16
            
            if (daisy.temperature > daisy.UT or daisy.temperature < daisy.LT):
                beta = 0.0
            else:
                beta = (1.0 - 0.003265 * ((((daisy.UT-daisy.LT)/2. + daisy.LT) - daisy.temperature)**2))
                
            dgrowth = daisy.area * self.X * beta
            
            if (daisy.area > daisy.seeds):
                ddeath = daisy.area * daisy.death_rate
            else:
                ddeath = 0.
            
            if (daisy.mutation > 0):
                daisy.mutate()
            
            
            # Quantity in each reservoir = quantity at previous timestep + (inflows - outflows) * dt    
            daisy.area += (dgrowth - ddeath) * dt
            
        self.time += dt
        
        return self.time, self.X, self.planet_temp, solflux
    
    

class Daisy(object):
    def __init__(self, name, albedo=0.5, area=0.001, death_rate=0.3, LT=5.0, UT=40.0, seeds=0.001, mutation=0):
        self.name = name
        self.albedo = albedo
        self.death_rate = death_rate
        self.UT = UT
        self.LT = LT
        self.seeds = seeds
        self.mutation = mutation
        
        self.temperature = 0.0
        self.area = area

    
    def __str__(self):
        return 'Daisy Object, name:{}'.format(self.name)
    
    def reset():
        self.temperature = [self.temperature[0]]
        self.beta = [self.beta[0]]
        self.area = [self.area[0]]         
    
    def mutate(self, rate=1.00):
        self.albedo = self.albedo if np.random.random() < self.mutation else self.albedo + self.albedo * rate * self.area * (np.random.random()*2 - 1)
        self.albedo = np.clip(self.albedo, 0., 1.)
        
        self.death_rate = self.death_rate if np.random.random() < self.mutation else self.death_rate + self.death_rate * rate * self.area * (np.random.random()*2 - 1)
        self.death_rate = np.clip(self.death_rate, 0., 1.)

        self.UT = self.UT if np.random.random() < self.mutation else self.UT + self.UT * rate * self.area * (np.random.random()*2 - 1)
        self.UT = np.clip(self.UT, -273.15, None)

        self.LT = self.LT if np.random.random() < self.mutation else self.LT + self.LT * rate * self.area * (np.random.random()*2 - 1)
        self.LT = np.clip(self.LT, -273.15, self.UT)




# -----------------
# main program
# -----------------

if __name__ == "__main__":
    
    
    # world building
    # -----
    
    daisyworld = World('daisyworld')
    white_daisy = Daisy('white', albedo=0.75, LT=5., UT=40.)
    black_daisy = Daisy('black', albedo=0.25, LT=5., UT=40.)
    grey_daisy = Daisy('grey', albedo=0.5, LT=5., UT=40.)    
    hyper_daisy = Daisy('hyper', albedo=1.0, LT=80., UT=100.)

    
    daisyworld.add_daisy(white_daisy)
    daisyworld.add_daisy(black_daisy)
    #daisyworld.add_daisy(grey_daisy)
    #daisyworld.add_daisy(hyper_daisy)
    
    # simulation loop
    # -----
    
        
    dt = 0.2
    tmax = 100.0
    imax = int((tmax + dt) / dt)
    
    planet_temp = np.zeros(imax)
    fertile_area = np.zeros(imax)
    daisydict = {}
    for daisy in daisyworld.daisies:
        daisydict[daisy.name] = np.zeros(imax)
    
    for i in range(imax):
        _, _, _, _ = daisyworld.simulate_timestep(dt)
        fertile_area[i], planet_temp[i] = daisyworld.X, daisyworld.planet_temp
        for daisy in daisyworld.daisies:
            daisydict[daisy.name][i] = daisy.area

    
    plot_world_trajectory(tmax, dt, daisydict, fertile_area, planet_temp)




# -----------------
# description of the variables
# -----------------

"""
The variables used in this model are as follows:


# reservoirs
# -----

area = The area of the planet covered by daisies, measured as a fraction of the
   total planetary area (unitless).

# fluxes
# -----

dgrowth = The white daisy growth flux is an increase in area over time.
ddeath = The white daisy death flux is a decrease in area over time.


# converters
# -----

dbeta = D_beta is the growth rate of daisies, which is a function of D_temp.
   In this model, the daisies, like most Earth species, have very specific ecological tolerances.
   Daisies cannot grow at temperatures colder than 5 degrees C or warmer than 40 degrees C.
   Their growth is maximized at 22.5 degrees.
wdtemp = Temperature of the white daisies population in degrees Celsius.
bdbeta = BD_beta is the growth rate of black daisies, which is a function of BD_temp.
   In this model, the daisies, like most Earth species, have very specific ecological tolerances.
   Daisies cannot grow at temperatures colder than 5 degrees C or warmer than 40 degrees C.
   Their growth is maximized at 22.5 degrees.
bdtemp = Temperature of the black daisies population in degrees Celsius.
X = The area of available fertile ground, not covered by either species, measured as a
   fraction of the total planetary area (dimensionless).
death_rate = The death rate of daisies (dimensionless).
fert = The proportion of the planet's area which is fertile ground (dimensionless).
daisyUT = The upper temperature limit of the daisies (40 degrees C).
daisyLT = The lower temperature limit of the daisies (5 degrees C).
alb_planet = The percentage of incoming solar radiation that is reflected back into space
   (dimensionless).
stefboltz = The Stefan Boltzmann constant {5.67e-8 (J/m^2*s*K^4)}.  The S-B constant relates
   the energy contained in a black body to its temperature.
solconst = The solar constant is the amount of light reaching the surface of the earth.  For
   the solar constant, Watson and Lovelock (1983) use a value of 9.17*105 ergs/cm2s, 
   or 917 J/(m^2*s).
wdalbedo = The percentage of incoming solar radiation that is reflected back into space by the
   white daisies (dimensionless).
bdalbedo = The percentage of incoming solar radiation that is reflected back into space by the
   black daisies (dimensionless).
gralbedo = The percentage of incoming solar radiation that is reflected back into space by the
   bare ground (dimensionless).
luminosity = Luminosity takes into account an increase in the amount of solar energy given off
   by the sun over time (dimensionless).  This equation begins with the sun at half its current
   brightness and ramps up the insolation by 2% a year.
solflux = Solar flux is the change in the amount of insolation reaching the earth over time.
Q = Q is a positive constant that expresses the degree to which absorbed solar energy is
   redistributed to the three types of surfaces on the planet (dimensionless).
planet_temp = Temperature of the entire planet in degrees Celsius.


# other
# -----

i = Counting loop incremental
imax = 
t = Time
tmax = The maximum time that the program runs to
dt = Time step in years.
"""

# _____________________________________________________________________________