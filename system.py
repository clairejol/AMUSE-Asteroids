#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 13:16:36 2022

@author: meijer
"""

from amuse.units import units as u
from amuse.units import constants as c
from amuse.lab import Particles
import numpy as np

class System:
    '''
    A class to house all the bodies in a solar system.

    Inputs:
    - config: a configuration .json file that details the particular system to be created.

    Returns:

    Variables:
    - star: the AMUSE particles object belonging to the star in the solar system.
    - planets: the list of all AMUSE particles objects belonging to the planets in the solar system.
    - asteroids: the list of all Asteroid class objects in the system.
    - observer: the location of the observer in the system.
    - observables: the list of all objects in the asteroids that are to be observed.
    - light_curve: the list of all flux calculations stored by calculate_flux for all observables.
    '''

    def __init__(self):
        '''
        Initializes the solar system with all of its objects.

        Inputs:
        - config: as defined in the class input.

        Returns:
        - System: a system class object.
        '''
        # minimum is sun, earth/observer, asteroid
        # properties we would want:
        # name, mass, radius, sma, orbital phase
        # make system coplanar
        system_info = {
            "stars" : { 
                "Sun" : {
                    "name"     : "Sun",
                    "mass"     : 1 | u.MSun,
                    "radius"   : 1 | u.RSun,
                    "semimajor_axis" : 0 | u.au,
                    "orbital_phase" : 0
                    # luminosity?
                    }
                },
            
            "planets" : { 
                "Earth" : {
                    "name"     : "Earth",
                    "mass"     : 1 | u.MEarth,
                    "radius"   : 1 | u.REarth,
                    "semimajor_axis" : 1 | u.au,
                    "orbital_phase" : 0
                    },
                
                "Jupiter" : {
                    "name"     : "Jupiter",
                    "mass"     : 1 | u.MJupiter,
                    "radius"   : 1 | u.RJupiter,
                    "semimajor_axis" : 5.2 | u.au,
                    "orbital_phase" : 2
                    }
                }
            }
        
        def pol2cart(r, phi):
            x = r * np.cos(phi)
            y = r * np.sin(phi)
            return(x.in_(u.au), y.in_(u.au), 0|u.au)
        
        
        stars = Particles(len(system_info["stars"]))
        planets = Particles(len(system_info["planets"]))
        particle_sets = {"stars" : stars, 
                         "planets" : planets}
        
        # properties that all bodies have
        for idx in system_info:
            bodies, bodies_info = particle_sets[idx], system_info[idx]
            i = 0
            for body in bodies_info:
                bodies[i].name   = bodies_info[body]["name"]
                bodies[i].mass   = bodies_info[body]["mass"]
                bodies[i].radius = bodies_info[body]["radius"]
                bodies[i].position = pol2cart(bodies_info[body]["semimajor_axis"], 
                                              bodies_info[body]["orbital_phase"])
                i += 1
        
        #print("stars:\n",stars)
        #print("planets\n",planets)
        
        def relative_orbital_velocity(mass, distance):
            return (c.G*mass/distance).sqrt()
        
        for planet in [planets[0]]:#planets:
            print(planet.position)
            vorb = relative_orbital_velocity(stars.mass.sum()+planet.mass, 
                                             planet.position.sum())
            print(vorb)
            print(planet.position)
            planet.velocity = [-planet.y, planet.x, planet.z]
            print(planet.velocity)
            #vec = pol2cart(bodies_info[body]["semimajor_axis"], 
            #               bodies_info[body]["orbital_phase"]+np.pi/2)
            #print(vec)
            #planet.velocity = (0, 1, 0) * vorb
        #print(planets)
        '''
        planets = Particles(2)
        earth = planets[0]
        earth.name = "Earth"
        earth.mass = 1 | u.MEarth
        earth.position = (1, 0, 0) | u.au
        vorb = relative_orbital_velocity(stars.mass.sum()+earth.mass, 
                                         earth.position.sum())
        earth.velocity = (0, 1, 0) * vorb
        jupiter = planets[1]
        jupiter.name = "Jupiter"
        jupiter.mass = 1 | u.MJupiter
        jupiter.position = (5.2, 0, 0) | u.au
        vorb = relative_orbital_velocity(stars.mass.sum()+jupiter.mass, 
                                         jupiter.position.sum())
        jupiter.velocity = (0, 1, 0) * vorb
        
        
        asteroids = Particles(1)
        asteroid = asteroids[0]
        asteroid.name = "asteroid"
        asteroid.mass = 7.34767309e+22 | u.kg
        asteroid.position = (384400, 0, 0) | u.km
        vorb = relative_orbital_velocity(earth.mass + asteroid.mass, 
                                         asteroid.position.sum())
        asteroid.velocity = (0, 1, 0) * vorb
        '''
        # merge these categories later, when they get evolved
    
    def get_directions():
        '''
        Given a coordinate in space, returns the vector direction to the observer and the star in the system.

        Inputs:
        - coord: an [x, y, z] coordinate in space.

        Returns:
        - obs_direction: the direction to the observer.
        - star_direction: the direction to the star. 
        '''
        # these should be normalised!
    
    def calculate_flux():
        '''
        Given an observable object, calculates the resultant total flux as observed by the observer thorough calling get_direction
        and subsequently Asteroid.get_flux. Stores the values in light_curve.

        Inputs:
        - observable: an Asteroid class observable object.

        Returns:
        - flux: the total flux observed by the observer.
        '''

    def get_gravity_at_point(): 
        '''
        The function to be called by the bridge to calculate accelerations on all observable objects. 
        Proceed as:
        Shape the inputs to [x, y, z] per observable.
        Update asteroid position with the coordinates.
        Calculate acceleration per asteroid in [ax, ay, az] format using Asteroid.get_acceleration per asteroid.
        Reshape to individual ax, ay, az and return.

        Inputs:
        - eps : a length parameter coming from the gravity solver, dummy variable.
        - x : the x-coordinates of observables coming from the gravity solver.
        - y : the y-coordinates of observables coming from the gravity solver.
        - z : the z-coordinates of observables coming from the gravity solver.

        Returns:
        - ax: the x-accelerations of all objects.
        - ay: the y-accelerations of all objects.
        - az: the z-accelerations of all objects.
        '''

system = System()



