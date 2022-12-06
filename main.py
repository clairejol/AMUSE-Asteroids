'''
Main file containing experiment information and plotting utilities.

@author: nazli
'''

import matplotlib.pyplot as plt
import seaborn as sb

from amuse.units import units as u
import numpy as np

from system import System
from bridge import Evolve

def get_config(): #Very crusty function to hold variable information. Cannot be easily made into a .json because of AMUSE units.

    experiment_config = {
        'time step': 0.05,
        'end time' : 10
    }

    system_config = {
        "stars" : { 
            "Sun" : {
                "name"     : "Sun",
                "mass"     : 1 | u.MSun,
                "radius"   : 1 | u.RSun,
                "semimajor_axis" : 0 | u.au,
                "orbital_phase" : 0,
                "luminosity" : 1 | u.LSun,
                },
            },
        
        "planets" : { 
            "Earth" : {
                "name"     : "Earth",
                "mass"     : 1 | u.MEarth,
                "radius"   : 1 | u.REarth,
                "semimajor_axis" : 1 | u.au,
                "orbital_phase" : np.pi/4,
                },
            
            "Jupiter" : {
                "name"     : "Jupiter",
                "mass"     : 1 | u.MJupiter,
                "radius"   : 1 | u.RJupiter,
                "semimajor_axis" : 5 | u.au,
                "orbital_phase" : np.pi/6,
                },
            
            },
        
        "asteroids" : {
            "Bennu" : {
                "name"     : "Bennu",
                "mass"     : 73e9 | u.kg,
                "radius"   : 0.24 | u.km,
                "semimajor_axis" : 1.126 | u.au,
                "orbital_phase" : np.pi/4,
                },
            }
        }
    return system_config, experiment_config

def plots(system, experiment_config):
    """___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___"""
    #Start with the position plot.
    fig, (ax_1) = plt.subplots(1,1)
    fig.set_size_inches(20, 20, forward=True)

    for item in ([ax_1.title, ax_1.xaxis.label, ax_1.yaxis.label] +
             ax_1.get_xticklabels() + ax_1.get_yticklabels()):
        item.set_fontsize(20)
    
    ax_1.set_ylabel('Heliocentric Y [AU]')
    ax_1.set_xlabel('Heliocentric X [AU]')

    ax_1.set_ylim(-2, 2)
    ax_1.set_xlim(-2, 2)

    #Manage the position and acceleration data to be plotted.
    asteroid_position = np.array(system.position_hist["asteroids"])
    star_position = np.array(system.position_hist["stars"])
    planet_positions = np.array(system.position_hist["planets"])

    acceleration_hist = np.array(system.acceleration_hist)
    acc_x = np.hstack((acceleration_hist[:,0][0], acceleration_hist[:,0][2:-2][::4], acceleration_hist[:,0][-1]))
    acc_y = np.hstack((acceleration_hist[:,1][0], acceleration_hist[:,1][2:-2][::4], acceleration_hist[:,1][-1]))
    acc_z = np.hstack((acceleration_hist[:,2][0], acceleration_hist[:,2][2:-2][::4], acceleration_hist[:,2][-1]))

    #Plot the asteroid position.
    ax_1.plot(asteroid_position[:,0], asteroid_position[:,1], color = 'grey', linewidth = 3.0, label = system.asteroids[0].name)

    #The planets' positions.
    for i in range(planet_positions.shape[2]):
        ax_1.plot(planet_positions[:,0,i], planet_positions[:,1,i], linewidth = 1.0, label = system.planets[i].name)

    #And finally the star.
    ax_1.scatter(star_position[0,0], star_position[0,1], color='orange', label = system.stars[0].name)

    #Plot the acceleration quivers.
    ax_1.quiver(asteroid_position[:,0], asteroid_position[:,1], acc_x, acc_y, color = 'grey', alpha=.5, label = 'YORP/Yarkovsky Acc.',
                    width = 0.005)

    plt.legend(loc = 'lower left')
    plt.show()
    """___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___"""
    #Add the semi-major axis plot.
    fig, (ax_2) = plt.subplots(1,1)
    fig.set_size_inches(20, 10, forward=True)

    for item in ([ax_2.title, ax_2.xaxis.label, ax_2.yaxis.label] +
             ax_2.get_xticklabels() + ax_2.get_yticklabels()):
        item.set_fontsize(20)
    
    ax_2.set_ylabel('Semi-major Axis [AU]')
    ax_2.set_xlabel('Time [years]')

    #Manage the semi-major axis data.
    semi_major_axis, ind = np.unique(system.semimajor_hist, return_index = True)
    semi_major_axis = semi_major_axis[np.argsort(ind)]
    times = np.arange(0, experiment_config['end time'], experiment_config['time step'])

    #Plot.
    ax_2.plot(times, semi_major_axis, color = 'black', linewidth = 3.0, label = system.asteroids[0].name)
    
    plt.legend(loc = 'lower left')
    plt.show()
    """___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___"""
    #Add the flux graph.
    fig, (ax_3) = plt.subplots(1,1)
    fig.set_size_inches(20, 10, forward=True)

    for item in ([ax_3.title, ax_3.xaxis.label, ax_3.yaxis.label] +
             ax_3.get_xticklabels() + ax_3.get_yticklabels()):
        item.set_fontsize(20)
    
    ax_3.set_ylabel(r'Flux Received from Observer [$W/(m^2)$]')
    ax_3.set_xlabel('Time [years]')

    #Manage the flux values.
    flux = np.array(system.light_curve)
    flux_value = np.zeros(len(flux))
    for i, row in enumerate(flux):
        flux_value[i] = row.number

    #Plot.
    ax_3.plot(times, flux_value, color = 'black', linewidth = 3.0, label = system.asteroids[0].name)

    plt.legend(loc = 'lower left')
    plt.show()
    """___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___"""
    return

def main(system_config, experiment_config):
    system = System(system_info = system_config)
    system = Evolve(system, experiment_config['time step'], experiment_config['end time'])
    plots(system = system, experiment_config = experiment_config)
    return

if __name__ == '__main__':
    system_config, experiment_config = get_config()
    main(system_config, experiment_config)