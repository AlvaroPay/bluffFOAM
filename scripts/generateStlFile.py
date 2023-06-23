# Created By  : Álvaro Pay Lozano
# Institution : TU Munich, Department of Aerospace and Geodesy
# Created Date: May 05, 2023
# version ='1.0'
# Description: From bluff body own coding input (5 digit) creates a general bluff body geometry
#              consisting on vertex and faces to export it a .stl file

import os
import numpy as np
from stl import mesh

class generateBluffBodies(object):
    def __init__ (self, caseDir, bluff_code):
        self.caseDir = str(caseDir)
        self.bluff_code = str(bluff_code)
        if len(self.bluff_code) > 5:
            print("Input needs to contain 5 characters: View coding system for more info")
            self.bluff_code = self.bluff_code[:5]

        elif len(self.bluff_code) < 5:
            print("Input needs to contain 5 characters: View coding system for more info")
            return

        self.shape = self.bluff_code[0]
        self.x_ref_length = 1
        self.y_ref_length = int(self.bluff_code[1]+self.bluff_code[2])
        self.AoA = np.radians(int(self.bluff_code[3])*2)
        self.edge = self.bluff_code[4]
        self.nPoints = 400

        if self.edge == "1":
            self.edge = 0.05*int(self.x_ref_length)
        elif self.edge == "2":
            self.edge = 0.1*int(self.x_ref_length)
        elif self.edge == "3":
            self.edge = 0.2*int(self.x_ref_length)
        elif self.edge == "4":
            self.edge = 0.25*int(self.x_ref_length)
        
        if self.y_ref_length <= 5:
            self.y_ref_length = 0.2*self.y_ref_length*self.x_ref_length
        else:
            self.y_ref_length = (self.y_ref_length/2 - 1.5)*self.x_ref_length
        
    def generateVertices(self):
        if self.shape == "1":
            print("circle")
            # Calculate circle points
            angles = np.linspace(0, 2*np.pi, self.nPoints, endpoint=False)
            x = self.x_ref_length * np.cos(angles)
            y = self.y_ref_length * np.sin(angles)
            vertices = np.column_stack((x, y))
            
            self.vertices = np.zeros((2 * vertices.shape[0], 3))
            self.vertices[:vertices.shape[0], 2] = 1
            self.vertices[vertices.shape[0]:, 2] = -1

            self.vertices[:vertices.shape[0], :2] = vertices
            self.vertices[vertices.shape[0]:, :2] = vertices

        elif self.shape == "2":
            print("rectangle")
            # Calculate rectangle points 
            if self.edge == "0":
                p = 2*self.x_ref_length+2*self.y_ref_length
                d = p/self.nPoints

                x = np.linspace(-self.x_ref_length/2, self.x_ref_length/2, round(self.x_ref_length/d), endpoint=False)
                y = np.ones(round(self.x_ref_length/d))*self.y_ref_length/2
                top_vertex = np.column_stack((x, y))

                # BOTTOM
                x = np.linspace(self.x_ref_length/2, -self.x_ref_length/2, round(self.x_ref_length/d), endpoint=False)
                y = np.ones(round(self.x_ref_length/d))*-self.y_ref_length/2
                bottom_vertex = np.column_stack((x, y))

                # RIGHT
                x = np.ones(round(self.y_ref_length/d))*self.x_ref_length/2
                y = np.linspace(self.y_ref_length/2, -self.y_ref_length/2, round(self.y_ref_length/d), endpoint=False)
                right_vertex = np.column_stack((x, y))

                # LEFT
                x = np.ones(round(self.y_ref_length/d))*-self.x_ref_length/2
                y = np.linspace(-self.y_ref_length/2, self.y_ref_length/2, round(self.y_ref_length/d), endpoint=False)
                left_vertex = np.column_stack((x, y))

                vertices = np.vstack((top_vertex, right_vertex, bottom_vertex, left_vertex))
                self.vertices = np.zeros((2 * vertices.shape[0], 3))
                self.vertices[:vertices.shape[0], 2] = 1
                self.vertices[vertices.shape[0]:, 2] = -1

                self.vertices[:vertices.shape[0], :2] = vertices
                self.vertices[vertices.shape[0]:, :2] = vertices
                
            else:
                
                p = 2*self.x_ref_length+2*self.y_ref_length+np.pi*self.edge
                d = p/(self.nPoints)                
                
                angle_r = np.linspace(np.pi/2, 0, round((np.pi*self.edge/4)/d), endpoint=False)
                angle_l = np.linspace(np.pi, np.pi/2, round((np.pi*self.edge/4)/d), endpoint=False)

                # TOP
                x = np.linspace(-self.x_ref_length/2+self.edge/2, self.x_ref_length/2-self.edge/2, round((self.x_ref_length-self.edge)/d), endpoint=False)
                y = np.ones(round((self.x_ref_length-self.edge)/d))*self.y_ref_length/2
                top = np.column_stack((x, y))

                x_r = self.x_ref_length/2-self.edge/2+self.edge/2 * np.cos(angle_r)
                y_r = self.y_ref_length/2-self.edge/2+self.edge/2 * np.sin(angle_r)
                x_l = -self.x_ref_length/2+self.edge/2+self.edge/2 * np.cos(angle_l)
                y_l = self.y_ref_length/2-self.edge/2+self.edge/2 * np.sin(angle_l)
                top_right = np.column_stack((x_r, y_r))
                top_left = np.column_stack((x_l, y_l))

                top_vertex = np.vstack((top_left, top, top_right))

                # BOTTOM
                x = np.linspace(self.x_ref_length/2-self.edge/2, -self.x_ref_length/2+self.edge/2, round((self.x_ref_length-self.edge)/d), endpoint=False)
                y = np.ones(round((self.x_ref_length-self.edge)/d))*-self.y_ref_length/2
                bottom = np.column_stack((x, y))

                angle_r = np.linspace(0, -np.pi/2, round((np.pi*self.edge/4)/d), endpoint=False)
                angle_l = np.linspace(-np.pi/2, -np.pi, round((np.pi*self.edge/4)/d), endpoint=False)

                x_r = self.x_ref_length/2-self.edge/2+self.edge/2 * np.cos(angle_r)
                y_r = -self.y_ref_length/2+self.edge/2+self.edge/2 * np.sin(angle_r)
                x_l = -self.x_ref_length/2+self.edge/2+self.edge/2 * np.cos(angle_l)
                y_l = -self.y_ref_length/2+self.edge/2+self.edge/2 * np.sin(angle_l)
                bottom_right = np.column_stack((x_r, y_r))
                bottom_left = np.column_stack((x_l, y_l))

                bottom_vertex = np.vstack((bottom_right, bottom, bottom_left))

                # RIGHT
                x = np.ones(round((self.y_ref_length-self.edge/2)/d))*self.x_ref_length/2
                y = np.linspace(self.y_ref_length/2-self.edge/2, -self.y_ref_length/2+self.edge/2, round((self.y_ref_length-self.edge/2)/d), endpoint=False)
                right_vertex = np.column_stack((x, y))

                # LEFT
                x = np.ones(round((self.y_ref_length-self.edge/2)/d))*-self.x_ref_length/2
                y = np.linspace(-self.y_ref_length/2+self.edge/2, self.y_ref_length/2-self.edge/2, round((self.y_ref_length-self.edge/2)/d), endpoint=False)
                left_vertex = np.column_stack((x, y))

                vertices = np.vstack((top_vertex, right_vertex, bottom_vertex, left_vertex))
                self.vertices = np.zeros((2 * vertices.shape[0], 3))
                self.vertices[:vertices.shape[0], 2] = 1
                self.vertices[vertices.shape[0]:, 2] = -1

                self.vertices[:vertices.shape[0], :2] = vertices
                self.vertices[vertices.shape[0]:, :2] = vertices

        elif self.shape == "3":
            print("triangle")
            if self.edge == "0":
                p = 2*(np.sqrt((self.x_ref_length)**2+(self.y_ref_length/2)**2))+self.y_ref_length
                d = p/self.nPoints
            
                x = np.linspace(-2*self.x_ref_length/3, self.x_ref_length/3, round(np.sqrt(self.x_ref_length**2+(self.y_ref_length/2)**2)/d), endpoint=False)
                y = np.linspace(0,self.y_ref_length/2, round(np.sqrt(self.x_ref_length**2+(self.y_ref_length/2)**2)/d), endpoint=False)
                first_vertex = np.column_stack((x, y))
            
                x = np.linspace(self.x_ref_length/3, -2*self.x_ref_length/3, round(np.sqrt(self.x_ref_length**2+(self.y_ref_length/2)**2)/d), endpoint=False)
                y = np.linspace(-self.y_ref_length/2, 0, round(np.sqrt(self.x_ref_length**2+(self.y_ref_length/2)**2)/d), endpoint=False)
                second_vertex = np.column_stack((x, y))

                x = np.ones(round(self.y_ref_length/d))*self.x_ref_length/3
                y = np.linspace(self.y_ref_length/2, -self.y_ref_length/2, round(self.y_ref_length/d), endpoint=False)
                third_vertex = np.column_stack((x, y))
            
                vertices = np.vstack((first_vertex, third_vertex, second_vertex))
                self.vertices = np.zeros((2 * vertices.shape[0], 3))
                self.vertices[:vertices.shape[0], 2] = 1
                self.vertices[vertices.shape[0]:, 2] = -1

                self.vertices[:vertices.shape[0], :2] = vertices
                self.vertices[vertices.shape[0]:, :2] = vertices
            
            else:
                p = 2*(np.sqrt((self.x_ref_length)**2+(self.y_ref_length/2)**2))+self.y_ref_length+np.pi*self.edge/4
                d = p/(self.nPoints)  
                theta = np.arctan2(self.y_ref_length/2,self.x_ref_length)
                angle_r = np.linspace(3*np.pi/2-theta, np.pi/2+theta, round((np.pi*self.edge/3)/d), endpoint=False)
            
                x = np.linspace(-2*self.x_ref_length/3+(self.x_ref_length/self.y_ref_length)*self.edge*np.cos(theta), self.x_ref_length/3, round(np.sqrt(self.x_ref_length**2+(self.y_ref_length/2)**2)/d), endpoint=False)
                y = np.linspace((self.x_ref_length/self.y_ref_length)*self.edge*np.sin(theta),self.y_ref_length/2, round(np.sqrt(self.x_ref_length**2+(self.y_ref_length/2)**2)/d), endpoint=False)
                first_vertex = np.column_stack((x, y))
            
                x = np.linspace(self.x_ref_length/3, -2*self.x_ref_length/3+(self.x_ref_length/self.y_ref_length)*self.edge*np.cos(theta), round(np.sqrt(self.x_ref_length**2+(self.y_ref_length/2)**2)/d), endpoint=False)
                y = np.linspace(-self.y_ref_length/2, -(self.x_ref_length/self.y_ref_length)*self.edge*np.sin(theta), round(np.sqrt(self.x_ref_length**2+(self.y_ref_length/2)**2)/d), endpoint=False)
                second_vertex = np.column_stack((x, y))

                x = np.ones(round(self.y_ref_length/d))*self.x_ref_length/3
                y = np.linspace(self.y_ref_length/2, -self.y_ref_length/2, round(self.y_ref_length/d), endpoint=False)
                third_vertex = np.column_stack((x, y))

                x_r = -2*self.x_ref_length/3+self.edge/(2*np.sin(theta))+self.edge/2*np.cos(angle_r)
                y_r = self.edge/2*np.sin(angle_r)
                rounding = np.column_stack((x_r, y_r))

                vertices = np.vstack((first_vertex, second_vertex, third_vertex, rounding))
                self.vertices = np.zeros((2 * vertices.shape[0], 3))
                self.vertices[:vertices.shape[0], 2] = 1
                self.vertices[vertices.shape[0]:, 2] = -1

                self.vertices[:vertices.shape[0], :2] = vertices
                self.vertices[vertices.shape[0]:, :2] = vertices


        elif self.shape == "4":
            print("semistadium")
            
            if self.edge == "0":
                #Semicircle
                p = np.pi*np.sqrt(((self.y_ref_length/2)**2+(self.x_ref_length/2)**2)/2)+2*self.x_ref_length+self.y_ref_length
                d = p/self.nPoints

                angles = np.linspace(np.pi/2, -np.pi/2, round((np.pi*np.sqrt(((self.y_ref_length/2)**2+(self.x_ref_length/2)**2)/2))/d), endpoint=False)
                x = -self.x_ref_length/2 * np.cos(angles) - self.x_ref_length/2 + 4/3*(self.x_ref_length/(2*np.pi))
                y = -self.y_ref_length/2 * np.sin(angles) 
                semi = np.column_stack((x,y))
                
                # TOP
                x = np.linspace(-self.x_ref_length/2 + 4/3*(self.x_ref_length/(2*np.pi)), self.x_ref_length/2 + 4/3*(self.x_ref_length/(2*np.pi)), round(self.x_ref_length/d), endpoint=False)
                y = np.ones(round(self.x_ref_length/d))*self.y_ref_length/2
                top_vertex = np.column_stack((x, y))

                # BOTTOM
                x = np.linspace(self.x_ref_length/2 + 4/3*(self.x_ref_length/(2*np.pi)), -self.x_ref_length/2 + 4/3*(self.x_ref_length/(2*np.pi)), round(self.x_ref_length/d), endpoint=False)
                y = np.ones(round(self.x_ref_length/d))*-self.y_ref_length/2
                bottom_vertex = np.column_stack((x, y))

                # RIGHT
                x = np.ones(round(self.y_ref_length/d))*self.x_ref_length/2 + 4/3*(self.x_ref_length/(2*np.pi))
                y = np.linspace(self.y_ref_length/2, -self.y_ref_length/2, round(self.y_ref_length/d), endpoint=False)
                right_vertex = np.column_stack((x, y))

                vertices = np.vstack((top_vertex, right_vertex, bottom_vertex, semi))
                self.vertices = np.zeros((2 * vertices.shape[0], 3))
                self.vertices[:vertices.shape[0], 2] = 1
                self.vertices[vertices.shape[0]:, 2] = -1

                self.vertices[:vertices.shape[0], :2] = vertices
                self.vertices[vertices.shape[0]:, :2] = vertices
            
            else:
                
                #Semicircle
                p = np.pi*np.sqrt(((self.y_ref_length/2)**2+(self.x_ref_length/2)**2)/2)+2*self.x_ref_length+self.y_ref_length-2*self.edge+np.pi*self.edge/2
                d = p/self.nPoints

                angles = np.linspace(np.pi/2, -np.pi/2, round((np.pi*np.sqrt(((self.y_ref_length/2)**2+(self.x_ref_length/2)**2)/2))/d), endpoint=False)
                x = -self.x_ref_length/2 * np.cos(angles) - self.x_ref_length/2 + 4/3*(self.x_ref_length/(2*np.pi))
                y = -self.y_ref_length/2 * np.sin(angles) 
                semi = np.column_stack((x,y))
                
                # TOP
                x = np.linspace(-self.x_ref_length/2 + 4/3*(self.x_ref_length/(2*np.pi)), self.x_ref_length/2 + 4/3*(self.x_ref_length/(2*np.pi))-self.edge/2, round((self.x_ref_length-self.edge/2)/d), endpoint=False)
                y = np.ones(round((self.x_ref_length-self.edge/2)/d))*self.y_ref_length/2
                top = np.column_stack((x, y))

                # BOTTOM
                x = np.linspace(self.x_ref_length/2-self.edge/2 + 4/3*(self.x_ref_length/(2*np.pi)), -self.x_ref_length/2 + 4/3*(self.x_ref_length/(2*np.pi)), round((self.x_ref_length-self.edge/2)/d), endpoint=False)
                y = np.ones(round((self.x_ref_length-self.edge/2)/d))*-self.y_ref_length/2
                bottom = np.column_stack((x, y))

                # ROUNDING
                angle_t = np.linspace(np.pi/2, 0, round((np.pi*self.edge/4)/d), endpoint=False)
                angle_b = np.linspace(0, -np.pi/2, round((np.pi*self.edge/4)/d), endpoint=False)

                x_t = self.x_ref_length/2-self.edge/2+4/3*(self.x_ref_length/(2*np.pi))+self.edge/2*np.cos(angle_t)
                y_t = self.y_ref_length/2-self.edge/2+self.edge/2*np.sin(angle_t)
                top_r = np.column_stack((x_t, y_t))

                x_b = self.x_ref_length/2-self.edge/2+4/3*(self.x_ref_length/(2*np.pi))+self.edge/2*np.cos(angle_b)
                y_b = -self.y_ref_length/2+self.edge/2+self.edge/2*np.sin(angle_b)
                bottom_r = np.column_stack((x_b, y_b))

                bottom_vertex = np.vstack((bottom_r, bottom))
                top_vertex = np.vstack((top, top_r))

                # RIGHT
                x = np.ones(round((self.y_ref_length-self.edge/2)/d))*self.x_ref_length/2 + 4/3*(self.x_ref_length/(2*np.pi))
                y = np.linspace(self.y_ref_length/2-self.edge/2, -self.y_ref_length/2+self.edge/2, round((self.y_ref_length-self.edge/2)/d), endpoint=False)
                right_vertex = np.column_stack((x, y))

                vertices = np.vstack((top_vertex, right_vertex, bottom_vertex, semi))
                self.vertices = np.zeros((2 * vertices.shape[0], 3))
                self.vertices[:vertices.shape[0], 2] = 1
                self.vertices[vertices.shape[0]:, 2] = -1

                self.vertices[:vertices.shape[0], :2] = vertices
                self.vertices[vertices.shape[0]:, :2] = vertices                
    
        if self.AoA != 0:
            if int(self.shape) != 1 or self.y_ref_length != 5:
                # Rotation of the points
                R = np.array([[np.cos(-self.AoA), -np.sin(-self.AoA), 0],
                [np.sin(-self.AoA), np.cos(-self.AoA),0],
                [0, 0, 1]])
                self.vertices = np.dot(R, self.vertices.T).T
        
        ar = self.y_ref_length 
    
    def generateFaces(self):       
        faces = np.zeros((self.vertices.shape[0], 3))

        for i in range(int(0.5 * faces.shape[0]) - 1):
            faces[2 * i] = [i, i + 1, i + int(0.5 * faces.shape[0])]
            faces[2 * i + 1] = [i + 1, i + 1 + int(0.5 * faces.shape[0]), i + int(0.5 * faces.shape[0])]

        faces[-2] = [int(0.5 * faces.shape[0]) - 1, 0, faces.shape[0] - 1]
        faces[-1] = [0, int(0.5 * faces.shape[0]), faces.shape[0] - 1]

        facesSide = np.zeros((self.vertices.shape[0] - 4, 3))

        for j in range(1, int(0.5 * facesSide.shape[0]) - 4):
            facesSide[2 * j - 1] = [j + 1, j, int(0.5 * self.vertices.shape[0]) - j]
            facesSide[2 * j] = [j + 1, int(0.5 * self.vertices.shape[0]) - j, int(0.5 * self.vertices.shape[0]) - j - 1]

        facesSide[0] = [1, 0, int(0.5 * self.vertices.shape[0]) - 1]
        facesSide[int(0.5 * facesSide.shape[0]) - 1] = [int(0.25 * self.vertices.shape[0]),
                                                        int(0.25 * self.vertices.shape[0]) - 1,
                                                        int(0.25 * self.vertices.shape[0]) + 1]

        facesSide[int(0.5 * facesSide.shape[0]):] = facesSide[:int(0.5 * facesSide.shape[0])] + int(
            0.5 * self.vertices.shape[0])

        faces = np.concatenate((faces, facesSide), axis=0)

        self.faces = faces.astype(int)
        print(self.faces)

    def generateStandardTriangleLanguageFile(self):
        bluff = mesh.Mesh(np.zeros(self.faces.shape[0], dtype=mesh.Mesh.dtype))
        for i, f in enumerate(self.faces):
            for j in range(3):
                bluff.vectors[i][j] = self.vertices[f[j], :]

        folderDir = os.path.join(self.caseDir, 'constant/triSurface')

        if not os.path.exists(folderDir):
            os.makedirs(folderDir)

        stlFileName = "bluff" + str(self.bluff_code) + ".stl"
        saveDir = os.path.join(folderDir, stlFileName)
        bluff.save(saveDir) 
   

       

                
