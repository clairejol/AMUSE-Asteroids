'''
This document aims to outline the overall structure of all the codebase that is developed to create and simulate a system of asteroids affected
by the YORP and Yarkovsky effects of their host star and possible other planets.
'''

def polygon_area(): ## Implemented ##
    '''
    Given a polygon and the normal of the plane defined by it, returns its area.

    Inputs:
    - polygon: a list of [x, y, z] coordinates that define the vertices of the polygon.
    - normal: a vector [x, y, z] that defines the direction of the surface normal.

    Returns: 
    - area: the total area of the polygon.
    '''

class Cube: ## Implemented ##
    '''
    A class for a cube object to be used as a base for all other shapes of objects.

    Inputs:
    - radius: the half-length of the edge of the resultant cube.
    - definition: the number of patches to be created per side edge of the cube. The number of patches per face of cube is given by
    (definition)^2.

    Returns:
    cube: if called, returns the list of [x, y, z] coordinates of its of point cloud.

    Variables:
    cube: as defined in the return.
    faces: a list of faces of the cube, each face containing the points on the given face.
    normals: the directions of the normal vectors per patch.
    '''
    
    @staticmethod
    def generate_cube(): ## Implemented ##
        '''
        Given a radius (aka. half of side length), generates a cube around (0,0,0).
        Can return the full cube or a list of inidividual faces. 

        Inputs: 
        - radius: as defined in the init for the class.
        - definition: as defined in the init for the class.
        - return_faces (boolean): if True, additionally retuns the faces of the cube as defined in variables.

        Returns:
        - cube: as defined in the returns of the class.
        - faces: as defined in the variables of the class.   
        '''
    
    def get_face_normals(): ## Implemented ##
        '''
        Given a face of a cube, returns the bases and the tips of all normal vectors on that face.

        Inputs: 
        - face: a face as defined in faces in the variables of the class.

        Returns:
        - norm_centers: the centers of the inidividual patches on the face, also the base of the normal vectors.
        - norm_tips: the tip points of all the normal vectors per patch on the face.
        '''

    def get_normals(): ## Implemented ##
        '''
        Given a cube, returns all of the normal vectors for all of the individual facets. 
        Can additionally return the base points of all normals.

        Inputs: 
        - return_points (boolean): if True, returns the center points and tips of all normals of the cube.

        Returns:
        - normals: the normal vectors for all the patches on all faces of the input cube.
        - centers, tips (optional): as defined in get_face_normals.
        '''
    
    def get_patches(): ## Implemented ## 
        '''
        Given a face of a cube, returns the arrays of corners of inidividual patches.

        Inputs:
        - face: a face as defined in the faces in the variables of the class.

        Returns:
        - patches: a list of patches on the input face, defined by the vertices [p0, p1, p2, p3] of a given square patch, 
        where every point is an [x, y, z] coordinate.
        '''
    
    def get_areas(): ## Implemented ##
        '''
        Given a a cube, return the areas of the individual polygons on all faces and the normal associated.

        Inputs:
        
        Returns:
        - result: a list of tuples [area, patch_normal] for all patches on the cube.
        '''

class Sphere(Cube): ## Implemented ##
    '''
    A sphere object built from a cube of same radius.

    Inputs:
    - radius: the radius of the resultant sphere.
    - definition: as defined in the parent Cube class.

    Returns:
    - sphere: if called, returns the list of points in its own sphere point cloud.

    Variables:
    - sphere: as defined in the return of the call.
    - normals: overwritten from the Cube class to be the normals of the patches on the sphere patches.
    '''

    @staticmethod
    def cube2sphere(): ## Implemented ##
        '''
        Convert from the surface of a cube to a sphere of same radius.
        Radius is determined automatically.

        Inputs:
        - vector: a point or collection of points in [x, y, z] coordinates assumed to be on a cube with automatically determined
        radius.

        Returns:
        - vector: the points converted to the corresponding points on the sphere with same radius.
        '''
    
    def generate_sphere(): ## Implemented ##
        '''
        Generates a sphere given a radius and definition through creating a cube and morphing into sphere.

        Inputs:
        - radius: as defined in the class inputs.
        - definition: as defined in the class inputs.

        Returns:
        - sphere: as defined in the class call returns.
        '''
    
    def get_sphere_normals(): ## Implemented ##
        '''
        Given a sphere, returns all of the normal vectors for the all of the individual facets. 
        Can additionally return the base points of all normals.

        Inputs:
        - cube_face: if not None, the function calculates the normals for only the inputted cube face through calling get_face_normals.
                     if None, the function calculates the normals through calling get_normals with return_points = True to get all centers
                     and tips of all normal vectors of all patches on the cube.
        - return_centers (boolean): if True, additionally returns the base points of all the normal vectors.
        
        Returns:
        - sphere_normals: the list of normal vectors of all the patches (or cube face) of the sphere.
        - sphere_centers (optional): the base points of all the normal vectors returned by the function.
        '''
    
    def get_sphere_areas(): ## Implemented ##
        '''
        The area calculation counterpart for the sphere.
        
        Inputs:

        Returns:
        - result: the list of tuples [sphere_patch_center, sphere_patch_normal, sphere_patch_area] for all the patches on the sphere.
        '''

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

    def __init__():
        '''
        Initializes the solar system with all of its objects.

        Inputs:
        - config: as defined in the class input.

        Returns:
        - System: a system class object.
        '''
    
    def get_directions():
        '''
        Given a coordinate in space, returns the vector direction to the observer and the star in the system.

        Inputs:
        - coord: an [x, y, z] coordinate in space.

        Returns:
        - obs_direction: the direction to the observer.
        - star_direction: the direction to the star. 
        '''
    
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

class Asteroid(amuse.lab.Particle):
    '''
    A subclass of the AMUSE Particle with added functionality for asteroids.

    Inputs:

    Variables:
    (from amuse.lab.Particle)
    - position : [x, y, z] spatial coordinates of the object.

    - tesselations: a list of patches [center, normal, area, temperature, albedo] defining the surface tesselation of the asteroid.
    '''

    def get_flux():
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
    
    def get_acceleration():
        '''
        Given the direction to and the star in the system, calculates the Yarkovsky (and YORP) forces on the asteroid through iterating over 
        its own tesselated surface and returns the accelerations per spatial coordinate.

        Inputs:
        - star_direction: as defined in System.get_directions.
        - star: as defined in System variables.
        
        Returns:
        - acceleration: the total acceleration vector [ax, ay, az] on the asteroid.
        '''
    

def Evolve():
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