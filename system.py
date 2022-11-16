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

    def __init__(self, system_info):
        '''
        Initializes the solar system with all of its objects. Assumes coplanar orbits.

        Inputs:
        - config: as defined in the class input.

        Returns:
        - System: a system class object.
        '''
                
        def cart2pol(x, y):
            rho = np.sqrt(x**2 + y**2)
            phi = np.arctan2(y, x)
            return(rho, phi)
        
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
        
        def relative_orbital_velocity(mass, distance):
            return (c.G*mass/distance).sqrt()
        
        if len(stars) > 1:
            # does it matter if a single star doesn't have a velocity?
            for star in stars:
                norm = (np.linalg.norm(star.position)).in_(u.au)
                vorb = relative_orbital_velocity(stars.mass.sum(), 
                                                 star.position.sum())
                star.vx = -star.y * vorb / norm
                star.vy =  star.x * vorb / norm
                star.vz =  star.z * vorb / norm
            
        for planet in planets:
            # this assumes that the system's center of mass is at [0,0,0]
            norm = (np.linalg.norm(planet.position)).in_(u.au)
            vorb = relative_orbital_velocity(stars.mass.sum()+planet.mass, 
                                             planet.position.sum())
                                             # something still feels off with "distance" here 
            planet.vx = -planet.y * vorb / norm
            planet.vy =  planet.x * vorb / norm
            planet.vz =  planet.z * vorb / norm
        # merge these three categories later, when they get evolved
    
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
            "semimajor_axis" : 5 | u.au,
            "orbital_phase" : np.pi/6
            }
        }
    }
system = System(system_info)



