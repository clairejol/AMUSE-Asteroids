import numpy as np
import itertools 

def arange(l, u, definition = 5):
    '''
    Modified version of the np.arange funvtion be inclusive and generate 
    points by a given definition.
    '''
    increment = (u-l)/definition
    return np.arange(l, u+increment-1e-5, increment)


def polygon_area(polygon, normal):
    '''
    Given a polygon and the normal of the plane defined by it, returns its area.
    '''
    total = [0, 0, 0]
    if isinstance(polygon, list):
        polygon = np.array(polygon)
    
    polygon = polygon[1:] - polygon[0:1]
    N = len(polygon)
    
    for i in range(N):
        vi1 = polygon[i]
        vi2 = polygon[(i+1) % N]
        prod = np.cross(vi1, vi2)
        total[0] += prod[0]
        total[1] += prod[1]
        total[2] += prod[2]
    result = np.dot(total, normal)
    return abs(result)



class Cube:
    '''
    A class for a cube object to be used as a base for all other shapes of objects.
    '''
    def __init__(self, radius, definition = 5):
        self.cube, self.faces = Cube.generate_cube(radius = radius, definition = definition, return_faces = True)
        self.normals = self.get_normals()

    def __call__(self):
        return self.cube
    
    @staticmethod
    def generate_cube(radius, definition, return_faces = False):
        '''
        Given a radius (aka. half of side length), generates a cube around (0,0,0).
        Can return the full cube or a list of inidividual faces.
        '''
        _constant = [radius]
        _neg_constant = [-radius]
        _range = arange(-1*radius, radius, definition)
        
        #Constant x.
        face_1 = np.array(list(itertools.product(_constant, _range, _range)))
        face_2 = np.array(list(itertools.product(_neg_constant, _range, _range)))
        #Constant y.
        face_3 = np.array(list(itertools.product(_range, _constant, _range)))
        face_4 = np.array(list(itertools.product(_range, _neg_constant, _range)))
        #Constant z.
        face_5 = np.array(list(itertools.product(_range, _range, _constant)))
        face_6 = np.array(list(itertools.product(_range, _range, _neg_constant)))

        cube = np.vstack((face_1, face_2, face_3, face_4, face_5, face_6))

        if return_faces:
            return cube, [face_1, face_2, face_3, face_4, face_5, face_6]
        
        else:
            return cube

    def get_face_normals(self, face): 
        '''
        Given a face of a cube, returns the bases and the tips of all normal vectors on that face.
        '''
        x = np.unique(face.T[0])
        y = np.unique(face.T[1])
        z = np.unique(face.T[2])

        bin_size = (x[1] - x[0])/2 if len(x) > 1 else (y[1] - y[0])/2

        if len(x) == 1:
            norm_centers = np.array(list(itertools.product(x, y[:-1]+bin_size, z[:-1]+bin_size)))
            norm_tips = norm_centers + np.tile(np.array([x[0]/np.abs(x[0]), 0, 0]), (norm_centers.shape[0], 1))

        elif len(y) == 1:
            norm_centers = np.array(list(itertools.product(x[:-1]+bin_size, y, z[:-1]+bin_size)))
            norm_tips = norm_centers + np.tile(np.array([0, y[0]/np.abs(y[0]), 0]), (norm_centers.shape[0], 1))
        
        elif len(z) == 1:
            norm_centers = np.array(list(itertools.product(x[:-1]+bin_size, y[:-1]+bin_size, z)))
            norm_tips = norm_centers + np.tile(np.array([0, 0, z[0]/np.abs(z[0])]), (norm_centers.shape[0], 1)) 
        
        return norm_centers, norm_tips

    
    def get_normals(self, return_points = False):
        '''
        Given a cube, returns all of the normal vectors for all of the individual facets. 
        Can additionally return the base points of all normals.
        '''
        _centers = []
        _tips = []

        for face in self.faces:
            centers, tips = self.get_face_normals(face)
            _centers.append(centers)
            _tips.append(tips)

        _centers = np.array(_centers).reshape(np.array(_centers).shape[0]*np.array(_centers).shape[1], np.array(_centers).shape[2])
        _tips = np.array(_tips).reshape(np.array(_tips).shape[0]*np.array(_tips).shape[1], np.array(_tips).shape[2])

        normals = _tips - _centers

        row_sums = np.array([np.linalg.norm(row) for row in normals])
        normals = np.array(normals / row_sums[:, np.newaxis])

        if return_points:
            return _centers, _tips
        
        else:
            return normals

    def get_patches(self, face):
        '''
        Given a face of a cube, returns the arrays of corners of inidividual patches.
        '''
        x = np.unique(face.T[0])
        y = np.unique(face.T[1])
        z = np.unique(face.T[2])

        bin_size = (x[1] - x[0]) if len(x) > 1 else (y[1] - y[0])

        if len(x) == 1:
            patches = []
            for point in itertools.product(x, y[:-1], z[:-1]):
                _p0 = np.array(list(point))
                _p1 = _p0 + np.array([0, bin_size, 0])
                _p2 = _p0 + np.array([0, 0, bin_size])
                _p3 = _p0 + np.array([0, bin_size, bin_size])
                patches.append([_p0, _p1, _p2, _p3])
        
        if len(y) == 1:
            patches = []
            for point in itertools.product(x[:-1], y, z[:-1]):
                _p0 = np.array(list(point))
                _p1 = _p0 + np.array([bin_size, 0, 0])
                _p2 = _p0 + np.array([0, 0, bin_size])
                _p3 = _p0 + np.array([bin_size, 0, bin_size])
                patches.append([_p0, _p1, _p2, _p3])
        
        if len(z) == 1:
            patches = []
            for point in itertools.product(x[:-1], y[:-1], z):
                _p0 = np.array(list(point))
                _p1 = _p0 + np.array([bin_size, 0, 0])
                _p2 = _p0 + np.array([0, bin_size, 0])
                _p3 = _p0 + np.array([bin_size, bin_size, 0])
                patches.append([_p0, _p1, _p2, _p3])
        
        return patches  

    def get_areas(self):
        '''
        Given a a cube, return the areas of the individual polygons on all faces and the normal associated.
        '''
        result = []

        for face in self.faces:
            patches = self.get_patches(face)
            for patch in patches:
                patch = np.array(patch)
                patch_normal_center, patch_normal_tip = self.get_face_normals(patch)
                patch_normal = patch_normal_tip - patch_normal_center
                patch_area = polygon_area(patch, patch_normal.squeeze())
                result.append([patch_area, patch_normal])
        
        return result



class Sphere(Cube):
    '''
    A sphere object built from a cube of same radius.
    '''
    def __init__(self, radius, definition = 5):
        self.sphere = self.generate_sphere(radius, definition)
        self.normals = self.get_sphere_normals()
        
    def __call__(self):
        return self.sphere
    
    @staticmethod
    def cube2sphere(vector):
        '''
        Convert from the surface of a cube to a sphere of same radius.
        Radius is determined automatically.
        '''
        x = vector.T[0]
        y = vector.T[1]
        z = vector.T[2]

        r_x = max(abs(x))
        r_y = max(abs(y))
        r_z = max(abs(z))
        r = max(r_x, r_y, r_z)

        #    r_x = abs(max(x) - min(x))
        #    r_y = abs(max(y) - min(y))
        #    r_z = abs(max(z) - min(z))
        #    r = max([r_x, r_y, r_z])/2

        x_ = r**(-1) * x * (r**2 - (y**2/2) - (z**2/2) + (y**2*z**2/(3*r**2)))**0.5
        y_ = r**(-1) * y * (r**2 - (x**2/2) - (z**2/2) + (x**2*z**2/(3*r**2)))**0.5
        z_ = r**(-1) * z * (r**2 - (x**2/2) - (y**2/2) + (x**2*y**2/(3*r**2)))**0.5

        return np.vstack((x_, y_, z_)).T

    def generate_sphere(self, radius, definition):
        self.cube, self.faces = Cube.generate_cube(radius, definition, return_faces = True)
        return Sphere.cube2sphere(self.cube)

    def get_sphere_normals(self, cube_face = None, return_centers = False):
        '''
        Given a sphere, returns all of the normal vectors for the all of the individual facets. 
        Can additionally return the base points of all normals. 
        '''
        if cube_face is None:
            cube_centers, cube_tips = self.get_normals(return_points = True)
        else:
            cube_centers, cube_tips = self.get_face_normals(cube_face)

        sphere_centers = self.cube2sphere(cube_centers)
        sphere_tips = self.cube2sphere(cube_tips)

        #sphere_normals = sphere_tips - sphere_centers #This does not work as the normal is not preserved.
        sphere_normals = sphere_centers

        row_sums = np.array([np.linalg.norm(row) for row in sphere_normals])
        sphere_normals = np.array(sphere_normals / row_sums[:, np.newaxis])

        if return_centers:
            return sphere_normals, sphere_centers
        
        else:
            return sphere_normals 

    def get_sphere_areas(self):
        '''
        The area calculation counterpart for the sphere.
        '''
        result = []

        for face in self.faces:
            patches = self.get_patches(face)
            for patch in patches:
                patch = np.array(patch)
                sphere_patch = Sphere.cube2sphere(patch)
                sphere_patch_normal, sphere_patch_center = self.get_sphere_normals(cube_face = patch, return_centers = True).squeeze()
                sphere_patch_area = polygon_area(sphere_patch, sphere_patch_normal)
                result.append([sphere_patch_center, sphere_patch_normal, sphere_patch_area])

        return result 