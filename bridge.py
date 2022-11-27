import numpy as np
from amuse.datamodel import ParticlesSuperset
from amuse.lab import nbody_system
from amuse.couple import bridge
from amuse.community.hermite.interface import Hermite

from amuse.units import units as u
from amuse.units import constants as c

def Evolve(system, timestep, endtime):
    '''
    A function to evolve a given system with the forces relevant in the project.

    Inputs:
    - system: a System object to be evolved.
    - timestep: the time interval for simulation steps.
    - endtime: the duration of time to be evolved.

    Returns:
    - system: the evolved system.

    Details to be determined, but assembles the n-body and gravity solvers, and assembles the bridge to incorporate the YORP and Yarkovsky
    forces on the system observables. Calls the system method to calculate and store observable fluxes.
    The bridge coupling should be a hierarchical bridge with coupling between (star + planets) <--> (asteroids) <--> (system).
    '''

    #Build the solar system first.
    solar_system = ParticlesSuperset([system.stars, system.planets])
    solar_converter = nbody_system.nbody_to_si(solar_system.mass.sum(), solar_system.position.length())
    
    ss_gravity_code = Hermite(solar_converter)
    ss_gravity_code.particles.add_particles(solar_system)
    ch_gravity2solar = ss_gravity_code.particles.new_channel_to(solar_system)

    #Now build the asteroid system.
    asteroid_sys = system.asteroids[0].particles_set
    ast_converter = nbody_system.nbody_to_si(asteroid_sys.mass.sum(),asteroid_sys.position.length())

    ast_gravity_code = Hermite(ast_converter)
    ast_gravity_code.particles.add_particles(asteroid_sys)
    ch_gravity2ast = ast_gravity_code.particles.new_channel_to(asteroid_sys)

    bridgey_bridge = bridge.Bridge()
    bridgey_bridge.add_system(ast_gravity_code, (ss_gravity_code,) )
    bridgey_bridge.add_code(ast_gravity_code, (system))
    bridgey_bridge.timestep = timestep | u.yr

    evolution_times = np.arange(0, endtime, timestep) | u.yr

    for time in evolution_times:
        #Evolve the model.
        bridgey_bridge.evolve_model(time)
        
        #Copy the properties back.
        ch_gravity2solar.copy()
        ch_gravity2ast.copy()

        #Calculate the flux through the system.
        system.calculate_flux()
        
    bridgey_bridge.stop()