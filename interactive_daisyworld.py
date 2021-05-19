"""
Daisyworld for Pythonista
Markus Roland Ernst, June 27th 2016

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

import numpy as np
import matplotlib.pyplot as plt
import ui
import time
from io import BytesIO

# UI connection Run Button


def prepare_daisyworld(sender):
    run_daisyworld(dt=0.2, t=0.0, tmax=100.0, death_rate=0.3, fert=float(view['valuelabel1'].text), daisyLT=float(view['valuelabel2'].text), daisyUT=float(
        view['valuelabel3'].text), solconst=917.0, luminosity=0.5, lumacc=0.02, wdalbedo=float(view['valuelabel4'].text), bdalbedo=float(view['valuelabel5'].text), gralbedo=0.5)

# UI update for the sliders


def update_sliders(sender):
    view['valuelabel1'].text = str(
        round(convert_slider_values(view['slider1'], 0, 1), 2))
    view['valuelabel2'].text = str(
        round(convert_slider_values(view['slider2'], 0, 100), 1))
    view['valuelabel3'].text = str(
        round(convert_slider_values(view['slider3'], 0, 100), 1))
    view['valuelabel4'].text = str(
        round(convert_slider_values(view['slider4'], 0, 1), 2))
    view['valuelabel5'].text = str(
        round(convert_slider_values(view['slider5'], 0, 1), 2))

# map 0,1 to minimum, maximum


def convert_slider_values(slider, mini, maxi):
    return (1 - slider.value) * mini + (slider.value) * maxi

# reset slider values to default


def reset_slider_values(sender):
    def animation():
        view['slider1'].value = 1.0
        view['slider2'].value = 0.05
        view['slider3'].value = 0.4
        view['slider4'].value = 0.75
        view['slider5'].value = 0.25
    ui.animate(animation, duration=1.0)
    update_sliders(view['slider1'])


# auto update (experimental)
def autoupdate(sender):
    time.sleep(1)
    #run_daisyworld(dt=0.2, t=0.0, tmax=100.0, death_rate=0.3, fert=float(view['valuelabel1'].text), daisyLT=float(view['valuelabel2'].text), daisyUT=float(view['valuelabel3'].text), solconst=917.0, luminosity=0.5, lumacc=0.02, wdalbedo=float(view['valuelabel4'].text), bdalbedo=float(view['valuelabel5'].text), gralbedo=0.5)


# The variables used in this model are as follows:
# RESERVOIRS:
# wdaisy = The area of the planet covered by white daisies, measured as a fraction of the
#    total planetary area (unitless).
# bdaisy = The area of the planet covered by black daisies, measured as a fraction of the
#    total planetary area (unitless).

# FLUXES:
# wdgrowth = The white daisy growth flux is an increase in wdaisy over time.
# wddeath = The white daisy death flux is a decrease in wdaisy over time.
# bdgrowth = The black daisy growth flux is an increase in bdaisy over time.
# bddeath = The black daisy death flux is a decrease in bdaisy over time.

# CONVERTERS:
# wdbeta = WD_beta is the growth rate of white daisies, which is a function of WD_temp.
#    In this model, the daisies, like most Earth species, have very specific ecological tolerances.
#    Daisies cannot grow at temperatures colder than 5 degrees C or warmer than 40 degrees C.
#    Their growth is maximized at 22.5 degrees.
# wdtemp = Temperature of the white daisies population in degrees Celsius.
# bdbeta = BD_beta is the growth rate of black daisies, which is a function of BD_temp.
#    In this model, the daisies, like most Earth species, have very specific ecological tolerances.
#    Daisies cannot grow at temperatures colder than 5 degrees C or warmer than 40 degrees C.
#    Their growth is maximized at 22.5 degrees.
# bdtemp = Temperature of the black daisies population in degrees Celsius.
# X = The area of available fertile ground, not covered by either species, measured as a
#    fraction of the total planetary area (dimensionless).
# death_rate = The death rate of daisies (dimensionless).
# fert = The proportion of the planet's area which is fertile ground (dimensionless).
# daisyUT = The upper temperature limit of the daisies (40 degrees C).
# daisyLT = The lower temperature limit of the daisies (5 degrees C).
# alb_planet = The percentage of incoming solar radiation that is reflected back into space
#    (dimensionless).
# stefboltz = The Stefan Boltzmann constant {5.67e-8 (J/m^2*s*K^4)}.  The S-B constant relates
#    the energy contained in a black body to its temperature.
# solconst = The solar constant is the amount of light reaching the surface of the earth.  For
#    the solar constant, Watson and Lovelock (1983) use a value of 9.17*105 ergs/cm2s,
#    or 917 J/(m^2*s).
# wdalbedo = The percentage of incoming solar radiation that is reflected back into space by the
#    white daisies (dimensionless).
# bdalbedo = The percentage of incoming solar radiation that is reflected back into space by the
#    black daisies (dimensionless).
# gralbedo = The percentage of incoming solar radiation that is reflected back into space by the
#    bare ground (dimensionless).
# luminosity = Luminosity takes into account an increase in the amount of solar energy given off
#    by the sun over time (dimensionless).  This equation begins with the sun at half its current
#    brightness and ramps up the insolation by 2% a year.
# solflux = Solar flux is the change in the amount of insolation reaching the earth over time.
# Q = Q is a positive constant that expresses the degree to which absorbed solar energy is
#    redistributed to the three types of surfaces on the planet (dimensionless).
# planet_temp = Temperature of the entire planet in degrees Celsius.

# OTHER VARIABLES:
# i = Counting loop incremental
# imax =
# t = Time
# tmax = The maximum time that the program runs to
# dt = Time step in years.


def run_daisyworld(dt=0.2, t=0.0, tmax=100.0, death_rate=0.3, fert=1.0, daisyLT=5.0, daisyUT=40.0, solconst=917.0, luminosity=0.5, lumacc=0.02, wdalbedo=0.75, bdalbedo=0.25, gralbedo=0.5):
    # dt=0.2             #timestep
    # t=0.0
    # tmax=100.0
    imax = int((tmax + dt) / dt)

    wdaisy = np.zeros(imax)
    bdaisy = np.zeros(imax)
    X = np.zeros(imax)
    planet_temp = np.zeros(imax)

    wdseeds = 0.001  # Arbitrary value used to limit the amount of depletion of wdaisy reservoir
    bdseeds = 0.001  # Arbitrary value used to limit the amount of depletion of bdaisy reservoir
    wdaisy[0] = 0.001  # Initial value of wdaisy: unitless
    bdaisy[0] = 0.001  # Initial value of bdaisy: unitless
    # death_rate=0.3      # unitless, but essentially per time
    # fert=1            # unitless
    # daisyUT=40.0        #degC
    # daisyLT=5.0         #degC

    # solconst=917.0              #J/(m^2*s)
    stefboltz = 5.67e-8  # J/(m^2*s*degK^4)
    #luminosity = 0.5
    #lumacc = 0.02
    # wdalbedo=0.75               #unitless
    # bdalbedo=0.25               #unitless
    # gralbedo=0.5                #unitless

    for i in range(1, imax):
        bdaisy[i] = bdaisy[i - 1]
        wdaisy[i] = wdaisy[i - 1]
        planet_temp[i] = planet_temp[i - 1]

        X[i] = fert - bdaisy[i] - wdaisy[i]  # unitless

        # ********************* planetary temperature**********************
        alb_planet = (bdalbedo * bdaisy[i]) + (wdalbedo *
                                               wdaisy[i]) + ((1. - bdaisy[i] - wdaisy[i]) * gralbedo)

        luminosity = luminosity + (lumacc * dt)  # unitless
        solflux = solconst * luminosity  # J/(m^2*s)
        q = (0.2 * solflux) / stefboltz  # unitless

        planet_temp[i] = (
            ((solflux * (1 - alb_planet)) / stefboltz)**0.25) - 273.15

        # ********************white daisy temp*****************************
        wdtemp = (((q * (alb_planet - wdalbedo)) +
                   ((planet_temp[i] + 273.16)**4))**(0.25)) - 273.16

        # ********************black daisy temp*****************************
        bdtemp = (((q * (alb_planet - bdalbedo)) +
                   ((planet_temp[i] + 273.16)**4))**(0.25)) - 273.16

        # ********************white daisy growth***************************
        if (wdtemp > daisyUT or wdtemp < daisyLT):
            wdbeta = 0.0
        else:
            wdbeta = (1.0 - 0.003265 *
                      ((((daisyUT - daisyLT) / 2. + daisyLT) - wdtemp)**2))

        wdgrowth = wdaisy[i] * X[i] * wdbeta

        if (wdaisy[i] > wdseeds):
            wddeath = wdaisy[i] * death_rate
        else:
            wddeath = 0.

        # Quantity in each reservoir = quantity at previous timestep + (inflows - outflows)*dt
        wdaisy[i] = wdaisy[i] + ((wdgrowth - wddeath) * dt)

        # ********************black daisy growth***************************
        if (bdtemp > daisyUT or bdtemp < daisyLT):
            bdbeta = 0.0
        else:
            bdbeta = (1.0 - 0.003265 *
                      ((((daisyUT - daisyLT) / 2. + daisyLT) - bdtemp)**2))

        bdgrowth = bdaisy[i] * X[i] * bdbeta

        if (bdaisy[i] > bdseeds):
            bddeath = bdaisy[i] * death_rate
        else:
            bddeath = 0.

        # Quantity in each reservoir = quantity at previous timestep + (inflows - outflows)*dt
        bdaisy[i] = bdaisy[i] + ((bdgrowth - bddeath) * dt)
        # print(i,t,bdaisy[i],wdaisy[i],planet_temp[i])    #print statement
        t = t + dt

    # plot

    xr = list(range(imax))
    tr = np.arange(0, tmax + dt, dt)

    sio = BytesIO()
    plt.subplot(2, 1, 1)
    plt.plot(tr, wdaisy, 'b-')
    plt.plot(tr, bdaisy, 'r-')
    plt.plot(tr, X, 'g-')
    plt.title('Daisyworld Output')
    plt.ylabel('% of area covered with daisies')

    plt.subplot(2, 1, 2)
    plt.plot(tr, planet_temp, 'r-')
    plt.xlabel('time(yr)')
    plt.ylabel('temperature of the planet')
    plt.savefig(sio)

    # send the plot to the interface
    view["image_result"].image = ui.Image.from_data(sio.getvalue())
    sio.close()
    plt.clf()


# show UI
view = ui.load_view()
# disable autoupdate
update_sliders(view['slider1'])
run_daisyworld()
view.present()
