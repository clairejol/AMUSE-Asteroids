import numpy
from amuse.lab import *

class Plummer_potential(object):

    def __init__(self, mass=1.0e+6 | units.MSun, R_Pl = 1|units.parsec):
        self.mass= mass
        self.radius = R_Pl
        self.epsilon2 = (0.1 | units.parsec)**2

    def get_potential_at_point(self,eps,x,y,z):
        r=(x**2+y**2+z**2 + self.radius**2 + self.epsilon2)**0.5
        potential = 2 * constants.G * self.mass/r
        return potential

    def get_gravity_at_point(self, eps, x,y,z):
        phi_0 = self.get_potential_at_point(eps, x,y,z)
        grav = AdaptingVectorQuantity()
        dpos = 0.001*(x**2+y**2+z**2).sqrt()
        phi_dx = self.get_potential_at_point(0,x+dpos,y,z) - phi_0
        phi_dy = self.get_potential_at_point(0,x,y+dpos,z) - phi_0
        phi_dz = self.get_potential_at_point(0,x,y, z+dpos) - phi_0
        return phi_dx/dpos, phi_dy/dpos, phi_dz/dpos

    def mass_in(self, r):
        return self.mass * r**3/(r**2 + self.radius**2)**(3./2.)

    def density(self, r):
        return 3.*self.mass/(4*numpy.pi*self.radius**3) * (1 + (r/self.radius)**2)**(-5./2.)

    def velocity_dispersion(self, R):
        r=(R**2 + self.radius**2)**0.5
        return numpy.sqrt(constants.G*self.mass_in(r)/(6*r))

    def circular_velocity(self, R):
        return numpy.sqrt(constants.G*self.mass_in(R)/R)
    
    def stop(self):
        return

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-M", unit=units.MSun,
                      dest="M_cl", type="float",default = 1.0e+6 | units.MSun,
                      help="cluster mass [%default]")
    result.add_option("-R", unit= units.parsec,
                      dest="R_Pl", type="float",default = 1.0 | units.parsec,
                      help="Plummer radius [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    Plummer =  Plummer_potential(o.M_cl, o.R_Pl)
    print("Mass in(", o.R_Pl.in_(units.parsec), ") = ", Plummer.mass_in(o.R_Pl).in_(units.MSun))
    print("velocity disp in(", o.R_Pl.in_(units.parsec), ") = ", Plummer.velocity_dispersion(o.R_Pl).in_(units.kms))


    
