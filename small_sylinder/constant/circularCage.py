"""
Author: Dr. Hui CHeng
Any questions about this code,
please email: hui.cheng@uis.no \n

"""

import numpy as np
import saveVtk as sv
# default global values, can be changed.
cr_top = 0.1/2.0  # [m]
cage_height = 0.1  # [m]
NT = 16  # 64
NN = 5  # 17

# private function
def __gen_points():
    point_one_cage = []
    for j in range(NN+1):
        for i in range(0, NT):
            point_one_cage.append(
                [cr_top * np.cos(i * 2 * np.pi / float(NT)),
                 cr_top * np.sin(i * 2 * np.pi / float(NT)),
                 - j * cage_height / float(NN)])
    return point_one_cage


def __gen_lines():    
    # for horizontal lines
    horizontal_lines = []
    for j in range(NN+1):
        for i in range(NT-1):
            horizontal_lines.append([i+j*NT, 1+i+j*NT])
        horizontal_lines.append([(j+1)*NT-1,j*NT])
    
    # vertical con for netting
    vertical_lines=[]
    for j in range(NN):
        for i in range(NT):
            vertical_lines.append([i+j*NT, i+(1+j)*NT])
    return horizontal_lines ,vertical_lines


def __gen_surfs():
    """ we assume the quad is numbered as [0,1,2,3]:
            0-------2
            |       |
            |       |
            1-------3
        
    """

    surf_element = []
    for j in range(NN):
        for i in range(NT-1):
            surf_element.append([i+j*NT, i+(1+j)*NT,
                               1+i+j*NT, 1+i+(1+j)*NT])
            
        surf_element.append([1+i+j*NT, 1+i+(1+j)*NT,
                                 j*NT, 1+1+i+j*NT])
    return surf_element

# public function
def gen_lines():
    hlines,vlines=__gen_lines()
    return hlines,vlines

def gen_cage():
    points = __gen_points()
    surfs = __gen_surfs()
    tris=__gen__tri(surfs)
    return points, tris

def __gen__tri(quad:list):
    """create a list of triangular based on quad that generated using __gen_surfs()

    Args:
        quad (list): we assume the quad is numbered as [0,1,2,3]:
            0-------2
            |       |
            |       |
            1-------3

    Returns:
        _type_: [[0,1,2],[3,2,1]]
    """
    tri_element=[]
    for item in quad:
        tri_element.append(item[:3])
        tri_element.append([item[3],item[2],item[1]])
    return tri_element
    
    
if __name__=="__main__":
    sv.point,sv.face=gen_cage()
    np.savetxt('points.out', sv.point, delimiter=',',fmt='%1u')
    np.savetxt('surfs.out',  sv.face, delimiter=',', fmt='%1u')
    
    sv.write_vtk("0")
    sv.save_element(sv.face)
    sv.save_point(sv.point)
    
    