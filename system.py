"""
File housing the System class which houses the individual particle sets for local identification.
Also calculates the flux observed and the YORP-Yarkovsky forces for asteroids through an augmented bridge.

@author: meijer, nazli
"""

from amuse.units import units as u
from amuse.units import constants as c
from amuse.lab import Particles
import numpy as np

from asteroid import Asteroid

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
    - system_info: a configuration .json file that details the particular system to be created.

    Returns:

    Variables:
    - stars: the AMUSE particles objects belonging to the stars in the solar system.
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
        - system_info: as defined in the class input.

        Returns:
        - System: a system class object.
        '''
        self.stars = Particles(len(system_info["stars"]))
        self.planets = Particles(len(system_info["planets"]))
        self.asteroids = [Asteroid(system_info["asteroids"]["Bennu"]["radius"])]
        
        particle_sets = {"stars" : self.stars, 
                         "planets" : self.planets,
                         "asteroids" : self.asteroids}

        self.light_curve = []
        
        # Update properties that are shared by all bodies.
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
        
        for star in self.stars:
            if len(self.stars) > 1:
                norm = (np.linalg.norm(star.position)).in_(u.au)
                vorb = relative_orbital_velocity(self.stars.mass.sum(), 
                                                 star.position.length())
                star.vx = -star.y * vorb / norm
                star.vy =  star.x * vorb / norm
                star.vz =  star.z * vorb / norm
            else:
                star.velocity = (0,0,0) | u.m * u.s**-1
            star.luminosity = system_info["stars"][star.name]["luminosity"]
        
        for planet in self.planets:
            #This assumes that the system's center of mass is at [0,0,0]
            norm = (np.linalg.norm(planet.position)).in_(u.au)
            vorb = relative_orbital_velocity(self.stars.mass.sum()+planet.mass, 
                                             planet.position.length())
                                        
            planet.vx = -planet.y * vorb / norm
            planet.vy =  planet.x * vorb / norm
            planet.vz =  planet.z * vorb / norm
        
        for asteroid in self.asteroids:
            #This assumes that the system's center of mass is at [0,0,0]
            norm = (np.linalg.norm(asteroid.position)).in_(u.au)
            vorb = relative_orbital_velocity(self.stars.mass.sum()+asteroid.mass, 
                                             asteroid.position.length())
                                             
            asteroid.vx = -asteroid.y * vorb / norm
            asteroid.vy =  asteroid.x * vorb / norm
            asteroid.vz =  asteroid.z * vorb / norm
        #Merge these three categories later, when they get evolved
        # should we move them to the center of mass somewhere? do the three need to be merged for that?
        
        #Assume the first planet is the observer by default.
        self.observer = self.planets[0] 
        self.observable = self.asteroids[0]
    
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
        
        obs_direction = self.observer.position - coord
        star_direction = self.stars[0].position - coord

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
        '''
        obs_direction, star_direction = self.get_directions(observable.position)
        flux = observable.get_flux(obs_direction, star_direction, self.observer, self.stars[0])
        self.light_curve.append(flux)

        return
        

    def get_gravity_at_point(self, eps, x, y, z): #Major rewrite necessary.
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
        #Update the location of the position coming from the simulation.
        self.observable.position = (x, y, z)

        #Calculate directions and acceleration.
        _, star_direction = self.get_directions(self.observable.position)
        acc = self.observable.get_acceleration(star_direction, self.stars[0])

        return acc[0], acc[1], acc[2]