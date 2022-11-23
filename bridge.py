import numpy
from amuse.datamodel import ParticlesSuperset
from matplotlib import pyplot
from amuse.units import units, constants
from amuse.lab import Particles
from amuse.lab import nbody_system
from amuse.couple import bridge
from amuse.community.hermite.interface import Hermite

from amuse.units import units as u
from amuse.units import constants as c
from amuse.lab import Particles
import numpy as np

#from asteroid import Asteroid

import matplotlib.pyplot as plt
from amuse.plot import plot, scatter

def Evolve(system, timestep):
    '''
    A function to evolve a given system with the forces relevant in the project.
    Inputs:
    - system: a System object to be evolved.
    - time: the duration of time to be evolved.
    Returns:
    - system: the evolved system.
    Details to be determined, but assembles the n-body and gravity solvers, and assembles the bridge to incorporate the YORP and Yarkovsky
    forces on the system observables. Calls the system method to calculate and store observable fluxes.
    The bridge coupling should be a hierarchical bridge with coupling between (star + planets) <--> (asteroids) <--> (system).
    '''

    planet_system = ParticlesSuperset([system.stars, system.planets])

    ss_converter=nbody_system.nbody_to_si(planet_system.mass.sum(),
                                       planet_system.position.length())
    ss_gravity_code = Hermite(ss_converter)
    ss_gravity_code.particles.add_particles(planet_system)
    ch_g2s = ss_gravity_code.particles.new_channel_to(planet_system)

    asteroid_sys = system.asteroids[0].particles_set

    ast_converter = nbody_system.nbody_to_si(asteroid_sys.mass.sum(),
                                       asteroid_sys.position.length())
    ast_gravity_code = Hermite(ast_converter)
    ast_gravity_code.particles.add_particles(asteroid_sys)
    ch_g2a = ast_gravity_code.particles.new_channel_to(asteroid_sys)

    gravity = bridge.Bridge()
    gravity.add_system(ast_gravity_code, (ss_gravity_code,) )
    gravity.add_code(ast_gravity_code, (INPUT))
    gravity.timestep = timestep | units.yr

    times = numpy.arange(0, 30, timestep) | units.yr
    x = [] | units.pc
    y = [] | units.pc
    for time in times:
        gravity.evolve_model(time)
        ch_g2l.copy()
        x.append(star[0].x)
        y.append(star[0].y)
    gravity.stop()
