'''
Main file containing experiment information and ANIMATED plotting utilities.

@author: nazli
'''
from amuse.units import units as u

import numpy as np
import os
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({'font.size': 20})
import matplotlib.animation as animation
from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt

from system import System
from bridge import Evolve

DATA_DIR = '/home/kutaynazli/Project/Movies'
DATA_DIR = '/net/vdesk/data2/meijer/SMA/AMUSE-Asteroids-figures'

def get_config(): 
    '''
    Very crusty function to hold information. Cannot be easily made into a .json because of AMUSE units.
    '''
    experiment_config = {
        'time step': 0.01,
        'end time' : 1.5
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
    '''
    Creates good-old matplotlib plots for the simulations.
    '''
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
    acc_x = np.hstack((acceleration_hist[:,0][0], acceleration_hist[:,0][1:-1][::2], acceleration_hist[:,0][-1]))
    acc_y = np.hstack((acceleration_hist[:,1][0], acceleration_hist[:,1][1:-1][::2], acceleration_hist[:,1][-1]))
    acc_z = np.hstack((acceleration_hist[:,2][0], acceleration_hist[:,2][1:-1][::2], acceleration_hist[:,2][-1]))

    #Plot the asteroid position.
    ax_1.plot(asteroid_position[:,0], asteroid_position[:,1], color = 'grey', linewidth = 3.0, label = system.asteroids[0].name)

    #The planets' positions.
    for i in range(planet_positions.shape[2]):
        ax_1.plot(planet_positions[:,0,i], planet_positions[:,1,i], linewidth = 1.0, label = system.planets[i].name)

    #And finally the star.
    ax_1.scatter(star_position[:,0], star_position[:,1], color='orange', label = system.stars[0].name)

    #Plot the acceleration quivers.
    ax_1.quiver(asteroid_position[:,0], asteroid_position[:,1], acc_x, acc_y, color = 'grey', alpha=.5, label = 'YORP/Yarkovsky Acc.',
                    width = 0.005)

    plt.legend(loc = 'lower left')
    plt.show()
    plt.savefig(os.path.join(DATA_DIR, 'position_plot.png'))

    """___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___"""

    #Add the semi-major axis plot.
    fig, (ax_2, ax_2b) = plt.subplots(2,1)
    fig.set_size_inches(20, 10, forward=True)

    for item in ([ax_2.title, ax_2.xaxis.label, ax_2.yaxis.label] +
             ax_2.get_xticklabels() + ax_2.get_yticklabels()):
        item.set_fontsize(20)
    
    ax_2.set_ylabel('Semi-major Axis [AU]')
    ax_2.set_xlabel('Time [years]')
    
    ax_2b.set_ylabel('Eccentricity')
    ax_2b.set_xlabel('Time [years]')

    #Manage the semi-major axis data.
    semi_major_axis, ind = np.unique(system.semimajor_hist, return_index = True)
    semi_major_axis = semi_major_axis[np.argsort(ind)]
    times = np.arange(0, experiment_config['end time'], experiment_config['time step'])

    position_norm, ind = np.unique(system.position_norm, return_index = True)
    position_norm = position_norm[np.argsort(ind)]

    eccentricity, ind = np.unique(system.eccentricity_hist, return_index = True)
    eccentricity = eccentricity[np.argsort(ind)]
    
    #Plot.
    ax_2.plot(times, semi_major_axis, color = 'black', linewidth = 3.0, label = "semimajor axis")#system.asteroids[0].name)
    #ax_2.plot(times, position_norm, color = 'blue', linewidth = 3.0, label = "distance from center")
    ax_2b.plot(times, eccentricity, color = 'red', linewidth = 3.0, label = "eccentricity")
    
    plt.legend(loc = 'lower left')
    plt.show()
    plt.savefig(os.path.join(DATA_DIR, 'semimajor_plot.png'))

    """___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___"""

    #Add the flux graph.
    fig, (ax_3) = plt.subplots(1,1)
    fig.set_size_inches(20, 10, forward=True)
    
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
    plt.savefig(os.path.join(DATA_DIR, 'flux_plot.png'))

    """___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___o___"""

    return



def anim_position_plot(system, experiment_config):
    '''
    Animated plot routine for the position.
    '''
    #Create a blank figure.
    fig, ax = plt.subplots()
    fig.set_size_inches(20, 20, forward=True)

    #Define the writer for the plots.
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=10, bitrate=3000)

    #Create the data to be plotted, starting with the positions.
    asteroid_position = np.array(system.position_hist["asteroids"])
    star_position = np.array(system.position_hist["stars"])
    planet_positions = np.array(system.position_hist["planets"])

    #Accelerations of the asteroid.
    acceleration_hist = np.array(system.acceleration_hist)
    acc_x = np.hstack((acceleration_hist[:,0][0], acceleration_hist[:,0][1:-1][::2], acceleration_hist[:,0][-1]))
    acc_y = np.hstack((acceleration_hist[:,1][0], acceleration_hist[:,1][1:-1][::2], acceleration_hist[:,1][-1]))
    acc_z = np.hstack((acceleration_hist[:,2][0], acceleration_hist[:,2][1:-1][::2], acceleration_hist[:,2][-1]))

    #And the time steps of the simulation.
    times = np.arange(0, experiment_config['end time'], experiment_config['time step'])

    #Make a vanishing line for the trail of positions.
    class vanishing_line(object):
        def __init__(self, x_data, y_data, tail_length, rgb_color):
            self.tail_length = tail_length
            self.rgb_color = rgb_color
            self.segments = self._set_segments(x_data, y_data)
            
        def _set_segments(self, x_data, y_data):          
            _xy_array = np.array([x_data, y_data]).T.reshape(-1, 1, 2)
            _segments = np.concatenate([_xy_array[:-1], _xy_array[1:]], axis = 1)
            
            return _segments
        
        def get_segments(self, frame):
            length = min(self.tail_length, frame)
            subsegments = self.segments[frame-length:frame]

            alphas = np.linspace(0., 1., length).reshape(length, 1)
            colors = np.tile(np.array(self.rgb_color), length).reshape(length, 3)

            colors_alphas = np.hstack((colors, alphas))

            subsegments_lc = LineCollection([], animated = True)
            subsegments_lc.set_segments(subsegments)
            subsegments_lc.set_color(colors_alphas)

            return subsegments_lc

    asteroid_line = vanishing_line(asteroid_position[:,0], asteroid_position[:,1], tail_length = 100, rgb_color = [.2, .2, .2])      
    
    def animate(frame):
        ax.clear()

        #Plot initial parameters.
        ax.set_title(f'Position of the Asteroid {system.asteroids[0].name} at time {times[frame]:.3f} years.')
        ax.set_ylabel('Heliocentric Y [AU]')
        ax.set_xlabel('Heliocentric X [AU]')

        ax.set_ylim(-2, 2)
        ax.set_xlim(-2, 2)
        
        #Plot the star's position, independent of frame.
        ax.scatter(star_position[frame,0], star_position[frame,1], color='red', label = system.stars[0].name)

        #Plot the planets's positions.
        colors = ['blue', 'orange']
        for i in range(planet_positions.shape[2]):
            ax.scatter(planet_positions[frame,0,i], planet_positions[frame,1,i], color= colors[i], label = system.planets[i].name)

        #Plot the asteroid position and acceleration.
        ax.scatter(asteroid_position[frame,0], asteroid_position[frame,1], color = 'grey', linewidth = 3.0, label = system.asteroids[0].name)
        ax.add_collection(asteroid_line.get_segments(frame))
        ax.quiver(asteroid_position[frame,0], asteroid_position[frame,1], acc_x[frame], acc_y[frame], 
                        color = 'grey', alpha=.5, width = 0.005, label = 'YORP/Yarkovsky Acc.')
        
        plt.legend(loc = 'lower left')
    
    #Animate the figure.
    animated_graph = animation.FuncAnimation(fig, animate, frames = len(times), interval = 2000, blit=False)
    animated_graph.save(os.path.join(DATA_DIR,'position_plot.mp4'), writer=writer)
    return



def anim_semimajor_plot(system, experiment_config):
    '''
    Animated plot routine for the semi major axis.
    '''
    #Create a blank figure.
    fig, ax = plt.subplots()
    fig.set_size_inches(20, 10, forward=True)

    #Define the writer for the plots.
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=10, bitrate=3000)

    #Create the data to be plotted, first the semimajor axis.
    semi_major_axis, ind = np.unique(system.semimajor_hist, return_index = True)
    semi_major_axis = semi_major_axis[np.argsort(ind)]

    #And the time steps of the simulation.
    times = np.arange(0, experiment_config['end time'], experiment_config['time step'])

    ax.set_xlim(times[0]-0.01, times[0]+0.01)
    ax.set_ylim(semi_major_axis[0]-0.01, semi_major_axis[0]+0.01)

    def animate(frame):
        ax.clear()

        ax.set_title(f'SMA of the Asteroid {system.asteroids[0].name}.')
        ax.set_ylabel('Semi-major Axis [AU]')
        ax.set_xlabel('Time [years]')

        #Plot the semi-major axis.
        ax.plot(times[:frame], semi_major_axis[:frame], color = 'black', linewidth = 3.0, label = system.asteroids[0].name)
    
        plt.legend(loc = 'lower left')

    #Animate the figure.
    animated_graph = animation.FuncAnimation(fig, animate, frames = len(times), interval = 2000, blit=False)
    animated_graph.save(os.path.join(DATA_DIR,'semimajor_plot.mp4'), writer=writer)
    return



def main(system_config, experiment_config):
    #Make and evolve the system.
    system = System(system_info = system_config)
    print('Successfully created the system. Starting time evolution.')
    system = Evolve(system, experiment_config['time step'], experiment_config['end time'])
    print('Time evolution complete. Plotting...')

    #Make the static plots.
    plots(system, experiment_config)
    print('Static plots completed.')

    #Make the position plot.
    #print('Creating position plot animation...')
    #anim_position_plot(system, experiment_config)
    #print('Complete.')

    #Make the semi-major axis plot.
    #print('Creating the semi-major axis plot animation...')
    #anim_semimajor_plot(system, experiment_config)
    #print('Complete.')



if __name__ == '__main__':
    system_config, experiment_config = get_config()
    main(system_config, experiment_config)