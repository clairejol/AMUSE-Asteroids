'''
File housing the Asteroid class, a subclass of AMUSE Particle with added functionality.
Calculates its own YORP and Yarkovsky forces in the System.

@author: nazli
'''

import numpy as np
from amuse.lab import Particle
from amuse.lab import constants
from amuse.units import units as u 

from cube_sphere import Sphere

class Asteroid(Particle):
    '''
    A subclass of the AMUSE Particle with added functionality for asteroids.

    Inputs:

    Variables:
    (from amuse.lab.Particle)
    - position : [x, y, z] spatial coordinates of the object.

    - tesselations: a list of patches [center, normal, area, temperature, albedo] defining the surface tesselation of the asteroid.
    '''
    def __init__(self, radius, definition = 10, key=None, particles_set=None, set_index=None, set_version=-1, **keyword_arguments):
        super().__init__(key, particles_set, set_index, set_version, **keyword_arguments)        
        self.radius = radius

        tessellations = Sphere(self.radius.value_in(u.m), definition = definition).get_sphere_areas()
        for face in tessellations:
            face += [0, 0, 0] #To store temperature, albedo, emissivity.
        tessellations = np.array(tessellations, dtype=object)

        object.__setattr__(self, "tessellations", tessellations)
        
        self.get_albedo_emissivity()
    
    def get_albedo_emissivity(self, albedo_array = None, emissivity_array = None):
        '''
        Either automatically initialize the surface emissivity and albedo, or import them from array.

        Inputs:
        - albedo_array: an array of albedo values with same length as the number of tessellations. Manual input or auto-generated.
        - emissivity_array: an array of emissivity values with same length as the number of tessellations. Manual input or auto-generated.

        Returns:
        None
        '''
        if not albedo_array: 
            self.tessellations[:,4] = np.full(len(self.tessellations[:,4]), 1) #np.random.uniform(0,1,len(self.tessellations[:,4]))
        else:
            self.tessellations[:,4] = albedo_array   
        if not emissivity_array:
            self.tessellations[:,5] = np.array([0.65]*len(self.tessellations[:,5]))
        else:
            self.tessellations[:,5] = emissivity_array

    def get_acceleration(self, star_direction, star):
        '''
        Given the direction to and the star in the system, calculates the Yarkovsky (and YORP) forces on the asteroid through iterating over 
        its own tesselated surface and returns the accelerations per spatial coordinate.

        Inputs:
        - star_direction: as defined in System.get_directions.
        - star: as defined in System variables.
        
        Returns:
        - acceleration: the total accelerations ax, ay, az on the asteroid.
        '''
        total_force = np.zeros(3) | u.m * u.kg * u.s**-2
        for patch in self.tessellations:
            #Define the variables. Ugly...
            center = patch[0]
            normal = patch[1]
            area = patch[2] | u.m**2
            temperature = None
            albedo = patch[4]
            emissivity = patch[5]
            
            #Calculate off-axis factor. If this number is less than 0, it means the patch is invisible from the star.
            mu_0 = np.dot(star_direction, normal)
            if mu_0 < 0:
                mu_0 = 0

            #Calculate incident flux.
            star_dist = np.linalg.norm(star.position-self.position) #norm is magnitude
            incident_flux = star.luminosity / (4*np.pi*star_dist**2)

            #Calculate and log the patch temperature.
            temperature = (((1-albedo)*mu_0*incident_flux)/(emissivity*constants.Stefan_hyphen_Boltzmann_constant))**(1/4)
            patch[3] = temperature
            
            #Calculate the scattering force.
            scattering_force = -(2/3) * (mu_0 * albedo * incident_flux / constants.c) * area * normal
            total_force += scattering_force

            #Calculate the thermal force.
            thermal_force = -(2/3) * (emissivity * constants.Stefan_hyphen_Boltzmann_constant \
                                      * temperature**4 / constants.c) * area * normal
            total_force += thermal_force

            #Return accelerations. Account for the the zero point float.
            zero_float = 5.45e-14 | u.m * u.s**-2
            zero = 0 | u.m * u.s**-2

            ax = total_force[0]/self.mass if np.abs(total_force[0]/self.mass) > zero_float else zero 
            ay = total_force[1]/self.mass if np.abs(total_force[1]/self.mass) > zero_float else zero 
            az = total_force[2]/self.mass if np.abs(total_force[2]/self.mass) > zero_float else zero 

            factor = 1e0 #Debugging tool.
            
        return factor*ax, factor*ay, factor*az 

    def get_flux(self, obs_direction, star_direction, observer, star):
        '''
        Given the directions and the observer and star in the system, calculates the flux observed by the observer through iterating 
        over its own tesselated surface.

        Inputs:
        - obs_direction: as defined in System.get_directions.
        - star_direction: as defined in System.get_directions.
        - observer: as defined in System variables.
        - star: as defined in System variables.

        Returns:
        - flux: the total flux observed by the observer.
        '''
        total_flux = 0 | u.kg * u.s**-3 
        for patch in self.tessellations:
            #Define the variables. Ugly...
            center = patch[0]
            normal = patch[1]
            area = patch[2] | u.m**2
            temperature = None
            albedo = patch[4]
            emissivity = patch[5]
            
            #Calculate off-axis factor. If this number is less than 0, it means the patch is invisible from the star.
            mu_0 = np.dot(star_direction, normal)
            if mu_0 < 0:
                mu_0 = 0 

            #Calculate the incident flux.
            star_dist = np.linalg.norm(star.position-self.position)
            incident_flux = star.luminosity / (4*np.pi*star_dist**2)

            #Calculate and log the patch temperature.
            temperature = (((1-albedo)*mu_0*incident_flux)/(emissivity*constants.Stefan_hyphen_Boltzmann_constant))**(1/4)
            patch[3] = temperature
            
            #Calculate the received flux due to re-emitted radiation. 
            obs_dist = np.linalg.norm(observer.position-self.position)
            re_emitted_flux = constants.Stefan_hyphen_Boltzmann_constant * area * temperature**4 / (2*np.pi*obs_dist**2)
            re_emitted_flux = re_emitted_flux * np.dot(obs_direction, normal)
            total_flux += re_emitted_flux
            
            #Define the reflection reception criterion.
            def reflection_reception(v1, v2, observer):
                observer_size = 2000 | u.REarth 
                obs_dist = np.linalg.norm(observer.position-self.position)
                angle = np.arctan(observer_size/obs_dist)
                return bool(np.arccos(np.dot(v1, v2)) < angle)

            #Calculate the received flux due to reflected radiation.
            reflected_flux = incident_flux * albedo
            reflected_direction = star_direction - 2*np.dot(star_direction,normal)*normal 
            if reflection_reception(obs_direction, reflected_direction, observer):
                total_flux += reflected_flux
                
        return total_flux