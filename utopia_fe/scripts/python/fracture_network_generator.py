import numpy as np 
import os
import array
from scanf import scanf
from ast import literal_eval

class CubItScript: 

    def __init__(self, path):
        self.n_vertices = 0
        self.n_curves = 0
        self.n_surfaces = 0
        self.n_volumes = 0
    
        if os.path.exists(path):
            os.remove(path)
        self.file = open(path, "w+")
        self.file.write("reset \n")

    def __del__(self):
        self.file.close()    

    def create_vertex(self, x,y,z):
        self.file.write("create vertex %f %f %f color \n" % (x,y,z) )
        
        self.n_vertices += 1
        id_vertex = self.n_vertices
        return (id_vertex) 

    def create_curve(self, vertices):
        n = len(vertices)

        self.file.write("create curve vertex ")
        
        for i in range(1, n+1):
            if i != n:
                self.file.write("%s," % vertices[i-1])
            else:
                self.file.write("%s" % vertices[i-1])

        self.file.write("\n")

        self.n_curves += 1 
        id_curve = self.n_curves
        return(id_curve) 

    def create_surface(self, vertices, sidesets=None):
        n = len(vertices) 

        self.file.write("create surface vertex ")
        
        for i in range(1, n+1):
            if i != n:
                self.file.write("%s," % vertices[i-1])
            else:
                self.file.write("%s" % vertices[i-1])

        self.file.write("\n")

        self.n_surfaces = self.n_surfaces+ 1 
        id_surface = self.n_surfaces


        if sidesets is not None:
            ns = len(sidesets)
            for i in range(0, ns):
                self.file.write("Sideset %d curve %d\n" %(sidesets[i][1], self.n_curves + sidesets[i][0]))


        self.n_curves += n
        return(id_surface) 

    def create_volume(self, coords, surface, sidesets = None, perpendicular=False):
        n = len(coords) 
        if perpendicular == True:
              self.file.write("sweep surface %d perpendicular distance 0.01\n" % (surface)) # perpendicular version, modify the desired distance here
        else: 

            self.file.write("sweep surface %d vector " % surface) # along vector version
        

            for i in range(1, n+1):
                if i != n:
                    self.file.write("%s," % coords[i-1])
                else:
                    self.file.write("%s\n" % coords[i-1])


        self.n_volumes = self.n_volumes+ 1 
        self.n_surfaces = self.n_surfaces+5
        return self.n_volumes

            

    def unite_surfaces(self, surfaces):
        n = len(surfaces) 

        self.file.write("unite volume ")
        
        for i in range(1, n+1):
            self.file.write("%s " % surfaces[i-1])

        self.file.write("\n")

        self.n_surfaces = self.n_surfaces + 1 
        id_surface = self.n_surfaces
        return Surface(id_surface) 

    def surface_block(self, surface, block):
        self.file.write("block %d surface %d\n" % (block, surface))

    def surface_unite_all(self):
        self.file.write("unite volume all\n")
        self.n_surfaces = self.n_surfaces + 1 
        id_surface = self.n_surfaces
        return(id_surface) 

    def surface_all_scale(self, scale_factor):
        self.file.write("surface all scale %f\n" % (scale_factor))

    # def surface_mesh_all(self)
    #     self.file.write("surface all scheme auto\n")
    #     self.file.write("surface all scale %f\n" % (scale_factor))



class Polygon:
    
    def __init__(self, points):
        self.points = points

    def create(self, script, sidesets=None):
        n = len(self.points)

        coord = array.array('l')

        for i in range(0, n):
            coord.append( script.create_vertex(self.points[i][0], self.points[i][1], 0.0) )

        return script.create_surface(coord, sidesets)

    def extrude(self, script, surface, eps, sidesets=None):
        if len(self.points[1]) == 2:
            n = np.array([0, 0, 1])
        else: 
            v0 = np.array([self.points[0][0], self.points[0][1], self.points[0][2]])
            v1 = np.array([self.points[1][0], self.points[1][1], self.points[1][2]])
            v2 = np.array([self.points[2][0], self.points[2][1], self.points[2][2]])

            u = (v1-v0)/np.linalg.norm(v1-v0)
            v = (v2-v0)/np.linalg.norm(v2-v0)

            n = np.cross(u,v)
            n = n/np.linalg.norm(n)

        return script.create_volume(n, surface)
        
class Surface:
    def __init__(self, id):
        self.id = id

    def extrude(self, script, dir, scale_factor = 1.0):
        dir[0] *= scale_factor
        dir[1] *= scale_factor
        dir[2] *= scale_factor

        script.create_volume(dir, self.id, perpendicular=True)
    

class Line: 

    def __init__(self, p0, p1):
        self.p0 = p0
        self.p1 = p1

    def create(self, script):
        v1 = script.create_vertex(p0, script)
        v2 = script.create_vertex(p1, script)
        vertices = np.array([v1,v2])
        return script.create_curves(vertices, script)

    def extrude(self, script, eps, sidesets=None):
        v1 = np.array([self.p0[0], self.p0[1]])
        v2 = np.array([self.p1[0], self.p1[1]])    

        m = (v2-v1)/np.linalg.norm(v2-v1)
        n = np.array([-m[1], m[0]]) 

        p1 = v1 - eps/2*n 
        p2 = v2 - eps/2*n 
        p3 = v2 + eps/2*n 
        p4 = v1 + eps/2*n 

        np.array([p1,p2,p3,p4])

        i1 = script.create_vertex(p1[0], p1[1], 0.0)
        i2 = script.create_vertex(p2[0], p2[1], 0.0)
        i3 = script.create_vertex(p3[0], p3[1], 0.0)
        i4 = script.create_vertex(p4[0], p4[1], 0.0)

        coord = np.array([i1, i2, i3, i4])
        return Surface(script.create_surface(coord, sidesets))


class CubItScriptReader:
    def __init__(self, path):
        file = open(path, "r")
        self.read(file)


    def read(self, file):
        pattern_vertex = 'create vertex %f %f %f color'
        pattern_curve  = 'create curve vertex %d %d'
        pattern_surface = 'create surface vertex %s'
        pattern_unite = 'unite volume %s'

        self.points = []
        self.curves = []
        self.surface = []
        self.unite = [] 
        self.id = []

        for line in file:
            point = scanf(pattern_vertex, line)
            
            if point is not None:
                self.points.append(point)
                
            curve = scanf(pattern_curve, line)
            
            if curve is not None:
                self.curves.append(curve)

            surface = scanf(pattern_surface, line)
          
            if surface is not None:  
                surface = list(map(literal_eval, surface))
                surface = np.array([tuple(i) for i in surface])
                for j in range(0, len(surface)):
                    self.surface.append(surface[j])
            
            
            unite = scanf(pattern_unite, line)
            if unite is not None:
                if unite[0] == 'all':
                    print('Error: \"unite volume all\" not supported, specify volume ids')
                    exit
                else:
                    unite = list(map(literal_eval, unite))
                    unite = np.array([tuple(i) for i in unite])
                    for j in range(0, len(unite)):
                        self.unite.append(unite[j])
                
            

    def extrude(self, cs, eps):
        for point in self.points:
            cs.create_vertex(point[0], point[1], point[2])

        for curve in self.curves:
            p1 = self.points[curve[0]-1]
            p2 = self.points[curve[1]-1]

            line = Line(p1, p2)
            line.extrude(cs, eps)

    def extrude_2D(self, cs):

        for point in self.points:
            cs.create_vertex(point[0], point[1], point[2])
        
        idx = []
        for surface in self.surface:  
            idx.append(cs.create_surface(surface))
        
        for unite in self.unite: 
           s = cs.unite_surfaces(unite)
           for i in range(0, len(unite)):
               idx.remove(unite[i])
           idx.append(s.id)

        for ids in idx:
            vol = Surface(ids)
            vol.extrude(cs, [0,0,1])
        


# cs = CubItScript("complex_matrix.txt")


# p0 = ([0.15,0.9167])
# p1 = ([0.4,0.5])
# l = Line(p0, p1)
# l.extrude(cs, 1e-4)


# p0 = ([0.65,0.8333])
# p1 = ([0.849723,0.167625])
# l = Line(p0, p1)
# l.extrude(cs, 1e-4)
            
            



# cs = CubItScript("fracture-network.txt")

# fracture_network = np.array([1, 2, 3, 4, 5, 6])

# p0 = ([0,0.5])
# p1 = ([1,0.5])
# l = Line(p0, p1)
# l.extrude(cs, 1e-4, [[2, 20], [4, 40]])


# p0 = ([0.5,0.0])
# p1 = ([0.5,1.0])
# l = Line(p0, p1)
# l.extrude(cs, 1e-4)


# p0 = ([0.5,0.75])
# p1 = ([1,0.75])
# l = Line(p0, p1)
# l.extrude(cs, 1e-4, [[2, 20]])


# p0 = ([0.75,0.5])
# p1 = ([0.75,1])
# l = Line(p0, p1)
# l.extrude(cs, 1e-4)



# p0 = ([0.5,0.625])
# p1 = ([0.75,0.625])
# l = Line(p0, p1)
# l.extrude(cs, 1e-4)


# p0 = ([0.625, 0.5])
# p1 = ([0.625, 0.75])
# l = Line(p0, p1)
# l.extrude(cs, 1e-4)


# poly = Polygon([[0, 0,0], [1, 0, 0], [1, 1, 0], [0, 1, 0]])
# surface_matrix = poly.create(cs, [[4, 4], [2, 2]])
# poly.extrude(cs, surface_matrix, 0.5)


# s = cs.unite_surfaces(fracture_network)
# cs.surface_block(s.id, 1)
# s.extrude(cs, [0, 0, 1])





# cs = CubItScript("pourous-matrix.txt")

# conv = CubItScriptReader("input.txt")
# conv.extrude(cs, 1e-2)
# cs.surface_all_scale(1000)
# cs.surface_unite_all();
# cs.surface_all_scale(0.001)


########### Test hydroicoin embedded from 2D to 3D ##################
# cs = CubItScript("hydrocoin_embedded.txt")
# input = CubItScriptReader("hydrocoin_equidim_embeddedfracture_script.txt")
# input.extrude_2D(cs)

########### Test real embedded from 2D to 3D ##################
# cs2 = CubItScript("real_embedded.txt")
# input = CubItScriptReader("fracture-network_script.txt")
# input.extrude_2D(cs2)

########## Test complex embedded from 2D to 3D ##################
cs2 = CubItScript("regular.txt")
input = CubItScriptReader("regular_equidim_nonconforming.txt")
input.extrude_2D(cs2)




