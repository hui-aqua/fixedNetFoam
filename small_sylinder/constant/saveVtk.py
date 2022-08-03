#!/usr/bin/env python
# -*- coding: utf-8 -*-

from vtk.vtkIOXML import vtkXMLDataSetWriter
import vtk


point=[[0,0,0],
   [1,1,1],
   [1,0,0],
   [0,1,0],
   [0,0,1]]

face=[[0,2,3],
   [0,2,4],
   [0,3,4]]


def __MakeMultiPoint(point_list:list):
    ug = vtk.vtkUnstructuredGrid()
    # A polyvertex is a cell represents a set of 0D vertices
    numberOfVertices = len(point_list)

    points = vtk.vtkPoints()
    polyVertex = vtk.vtkPolyVertex()
    polyVertex.GetPointIds().SetNumberOfIds(numberOfVertices)

    for index, item in enumerate(point_list):
        points.InsertNextPoint(item)
        polyVertex.GetPointIds().SetId(index, index)

    
    ug.SetPoints(points)
    ug.InsertNextCell(polyVertex.GetCellType(), polyVertex.GetPointIds())

    return ug


def __MakeMultiFace(point_list:list,triangle_list:list):
    ug = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()
    triangles = vtk.vtkTriangle()  # 3-point elements
    pixel=vtk.vtkPixel() # 4-point elements
    
    # add points   
    for i in range(len(point_list)):
        points.InsertNextPoint(point_list[i])
    ug.SetPoints(points)
    
    # add elements
    for item in triangle_list:
        if len(item)==3: # add 3-point elements
            for i in range(3): 
                triangles.GetPointIds().SetId(i, item[i])
            ug.InsertNextCell(triangles.GetCellType(), triangles.GetPointIds())    
        elif len(item)==4:  # add 4-point elements
            for i in range(4): 
                pixel.GetPointIds().SetId(i,item[i])
            ug.InsertNextCell(pixel.GetCellType(), pixel.GetPointIds())    
    return ug

def write_vtk(file_name:str):
    u1=__MakeMultiPoint(point)
    u3=__MakeMultiFace(point,face)
    
    writer = vtkXMLDataSetWriter()
    #u1
    writer.SetInputData(u1)
    writer.SetFileName(file_name+'.point.vtu')
    writer.Write()
   
    writer.SetInputData(u3)
    writer.SetFileName(file_name+'.face.vtu')
    writer.Write()    


def save_element(element):
    """
    :param hydro_element: A numpy array of all the net panel element
    :param cwd: work path,
    :return: write the net panels to "constant" folder as a file named "surf"
    """
    head_file = ["FoamFile",
                 "{",
                 "    version     1906;",
                 "    format      ascii;",
                 "    class       dictionary;",
                 "    location    \"constant\";",
                 "    object      element;",
                 "}",
                 "// >>>>>>>>>>>>>>>>>>>>>>>>>>Hui Cheng>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //",
                 "{",
                 "numOfSurf   " + str(len(element)) + ";"]
    with open('constant/netElement', 'w') as output_file:
        for item in head_file:
            output_file.write(item + "\n")
        for index, item in enumerate(element):
            output_file.write(
                "e" + str(index) + "\t( " + str(item[0]) + "\t" + str(item[1]) + "\t" + str(item[2]) + " );\n")
        output_file.write("}\n")
    output_file.close()


def save_point(position):
    """
    :param position: A numpy array of all the nodes' position
    :param cwd: work path,
    :return: write the nodes' positions to "constant" folder as a file named "posi"
    """
    
    head_file = [
                 "FoamFile",
                 "{",
                 "    version     1906;",
                 "    format      ascii;",
                 "    class       dictionary;",
                 "    location    \"constant\";",
                 "    object      position;",
                 "}",
                 "// >>>>>>>>>>>>>>>>>>>>>>>>>>Hui Cheng>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //",
                 "{",
                 "numOfPoint   " + str(len(position)) + ";"]
    with open('constant/netPoint', 'w') as output_file:
        for item in head_file:
            output_file.write(item + "\n")
        for index, item in enumerate(position):
            output_file.write(
                "p" + str(index) + "\t( " + str(item[0]) + "\t" + str(item[1]) + "\t" + str(item[2]) + " );\n")
        output_file.write("}\n")    
    output_file.close()
    