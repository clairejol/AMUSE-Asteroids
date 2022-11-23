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

from asteroid import Asteroid

import matplotlib.pyplot as plt
from amuse.plot import plot, scatter

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def pol2cart(r, phi):
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    return(x.in_(u.au), y.in_(u.au), 0|u.au)

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
        self.stars = Particles(len(system_info["stars"]))
        self.planets = Particles(len(system_info["planets"]))
        self.asteroids = [Asteroid(system_info["asteroids"]["Bennu"]["radius"])]
        
        particle_sets = {"stars" : self.stars, 
                         "planets" : self.planets,
                         "asteroids" : self.asteroids}
        
        # properties that all bodies have
        for idx in system_info:
            bodies, bodies_info = particle_sets[idx], system_info[idx]
            for i, body in enumerate(bodies_info):
                bodies[i].name   = bodies_info[body]["name"]
                bodies[i].mass   = bodies_info[body]["mass"]
                bodies[i].radius = bodies_info[body]["radius"]
                bodies[i].position = pol2cart(bodies_info[body]["semimajor_axis"], 
                                              bodies_info[body]["orbital_phase"])
        
        def relative_orbital_velocity(mass, distance):
            return (c.G*mass/distance).sqrt()
        
            # does it matter if a single star doesn't have a velocity?
        for star in self.stars:
            if len(self.stars) > 1:
                norm = (np.linalg.norm(star.position)).in_(u.au)
                vorb = relative_orbital_velocity(self.stars.mass.sum(), 
                                                 star.position.sum())
                star.vx = -star.y * vorb / norm
                star.vy =  star.x * vorb / norm
                star.vz =  star.z * vorb / norm
            star.luminosity = system_info["stars"][star.name]["luminosity"]
        
        for planet in self.planets:
            # this assumes that the system's center of mass is at [0,0,0]
            norm = (np.linalg.norm(planet.position)).in_(u.au)
            vorb = relative_orbital_velocity(self.stars.mass.sum()+planet.mass, 
                                             planet.position.sum())
                                             # something still feels off with "distance" here 
            planet.vx = -planet.y * vorb / norm
            planet.vy =  planet.x * vorb / norm
            planet.vz =  planet.z * vorb / norm
        
        for asteroid in self.asteroids:
            # this assumes that the system's center of mass is at [0,0,0]
            # it also is the exact same as for planets, at the moment
            norm = (np.linalg.norm(asteroid.position)).in_(u.au)
            vorb = relative_orbital_velocity(self.stars.mass.sum()+asteroid.mass, 
                                             asteroid.position.sum())
                                             # something still feels off with "distance" here 
            asteroid.vx = -asteroid.y * vorb / norm
            asteroid.vy =  asteroid.x * vorb / norm
            asteroid.vz =  asteroid.z * vorb / norm
        # merge these three categories later, when they get evolved
        # should we move them to the center of mass somewhere? do the three need to be merged for that?
                
        self.observer = self.planets[0] # assume the first planet in the list is the observer
        
        self.calculate_flux(self.asteroids[0])
        self.get_gravity_at_point(0,0,0,0)
    
    def get_directions(self, coord):
        '''
        Given a coordinate in space, returns the vector direction to the observer 
        and the star in the system. Assumes a single star.
        Inputs:
        - coord: an [x, y, z] coordinate in space.
        Returns:
        - obs_direction: the direction to the observer, normalised.
        - star_direction: the direction to the star, normalised. 
        '''
        try:
            coord = coord.in_(u.au)
        except:
            raise ValueError("The coordinate should have a distance unit.")
        
        obs_direction = coord - self.observer.position
        star_direction = coord - self.stars[0].position
        return obs_direction / (np.linalg.norm(obs_direction)).in_(u.au), \
              star_direction / (np.linalg.norm(star_direction)).in_(u.au)
    
    def calculate_flux(self, observable):
        '''
        Given an observable object, calculates the resultant total flux as 
        observed by the observer thorough calling get_direction and 
        subsequently Asteroid.get_flux. Stores the values in light_curve.
        Inputs:
        - observable: an Asteroid class observable object.
        Returns:
        - flux: the total flux observed by the observer.
        '''
        obs_direction, star_direction = self.get_directions(observable.position)
        flux = observable.get_flux(obs_direction, star_direction, self.observer, self.stars[0])
        # currently returns nan, because emissivity in Asteroid is set to 0

        return flux
        

    def get_gravity_at_point(self, eps, x, y, z): 
        '''
        The function to be called by the bridge to calculate accelerations on all observable objects. 
        Proceed as:
        Shape the inputs to [x, y, z] per observable.
        Update asteroid position with the coordinates.
        Calculate acceleration per asteroid in [ax, ay, az] format using 
            Asteroid.get_acceleration per asteroid.
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
        observable = self.asteroids[0]
        obs_direction, star_direction = self.get_directions(observable.position)
        acc = observable.get_acceleration(star_direction, self.stars[0])
        return acc[0], acc[1], acc[2]
        # take actual gravity into account as well!

system_info = {
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
            "orbital_phase" : 0,
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
#system = System(system_info)
#print(system.get_directions((.5,.5,0.)|u.au))
