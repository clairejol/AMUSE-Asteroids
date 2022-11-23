import os 
import numpy as np
import itertools
from amuse.lab import Particle, Particles
from amuse.lab import constants
from amuse.units import units as u 

from cube_sphere import Sphere

class Asteroid(Particle):
    '''
    A modified Particle-class to add extra functionality for asteroid with tessellated faces.
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
        '''
        if not albedo_array: 
            self.tessellations[:,4] = np.random.uniform(0,1,len(self.tessellations[:,4]))
        else:
            self.tessellations[:,4] = albedo_array   
        if not emissivity_array:
            self.tessellations[:,5] = np.array([0.65]*len(self.tessellations[:,5]))
        else:
            self.tessellations[:,5] = emissivity_array

    def get_acceleration(self, star_direction, star):
        '''
        Calculate the acceleration on a given asteroid due to YORP and Yarkovsky forces.
        '''
        total_force = np.zeros(3) | u.m * u.kg * u.s**-2 # it needs units for the adding to work
        for patch in self.tessellations:
            #Define the variables. Ugly...
            center = patch[0]
            normal = patch[1]
            area = patch[2] | u.m**2
            temperature = None
            albedo = patch[4]
            emissivity = patch[5]
            print()
            print("patch center:", np.sign(center))
            print("patch normal:", np.sign(normal))
            
            #Calculate off-axis factor.
            mu_0 = np.dot(star_direction, normal)
            if mu_0 < 0:
                mu_0 = 0
                # might need different solution if a heated patch rotates out of the sunlight
                # because the temperature does not instantly go to 0 then

            #Calculate incident flux.
            star_dist = np.linalg.norm(star.position-self.position)
            incident_flux = star.luminosity / (4*np.pi*star_dist**2)

            #Calculate and log the patch temperature.
            temperature = (((1-albedo)*mu_0*incident_flux)/(emissivity*constants.Stefan_hyphen_Boltzmann_constant))**(1/4)
            patch[3] = temperature
            
            #Calculate the scattering force.
            scattering_force = -(2/3) * (mu_0 * albedo * incident_flux / constants.c) * area * normal
            total_force += scattering_force
            print("scattering_force:",(scattering_force))
            #print("scattering_force:",np.sign(scattering_force))
            
            #Calculate the thermal force.
            thermal_force = -(2/3) * (emissivity * constants.Stefan_hyphen_Boltzmann_constant \
                                      * temperature**4 / constants.c) * area * normal
            total_force += thermal_force
            #print("thermal_force:",np.sign(thermal_force))
            print("thermal_force:",(thermal_force))
            
        return -total_force[0]/self.mass, -total_force[1]/self.mass, -total_force[2]/self.mass 
               # i put minus signs in front of these by hand, because it puts out 
               # a force towards the star, not away from it as we expect it to,
               # but there's probably a better solution

    def get_flux(self, obs_direction, star_direction, observer, star):
        '''
        Calculate the total flux observable from the observer for the asteroid.
        '''
        total_flux = 0 | u.kg * u.s**-3 # it needs units for the adding of fluxes to work
        i = 0
        for patch in self.tessellations:
            #Define the variables. Ugly...
            center = patch[0]
            normal = patch[1]
            area = patch[2] | u.m**2
            temperature = None
            albedo = patch[4]
            emissivity = patch[5]
            
            #Calculate off-axis factor.
            mu_0 = np.dot(star_direction, normal) 

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
                return bool(np.arctan(np.dot(v1, v2)) < angle)

            #Calculate the received flux due to reflected radiation.
            reflected_flux = incident_flux * albedo
            reflected_direction = star_direction - 2*np.dot(-star_direction,normal)*normal 
            if reflection_reception(obs_direction, reflected_direction, observer):
                i+=1
                total_flux += reflected_flux
                
        return total_flux