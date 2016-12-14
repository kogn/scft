#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
$ mayavi2 -x structured_grid.py

"""

import numpy as np
from numpy import cos, sin, pi,array
import tvtk
from tvtk.api import tvtk
from mayavi.scripts import mayavi2


# dims = (1001, 1001, 2)
# # dims = (128, 201, 33)

# sgrid = tvtk.StructuredGrid(dimensions=dims)
# sgridb = tvtk.StructuredGrid(dimensions=dims)
# sgride = tvtk.StructuredGrid(dimensions=dims)
# sgrid_det = tvtk.StructuredGrid(dimensions=dims)
# ################################################################################ 文件名
# file1 = open(r"polyvertical_t_-1.79_R_24.52_D_4.00_I_1001_J_1001_K_2_eta_1039.3_L21_+0.0_s0_1.00_test_D_4_DD_U_B_point.txt","r")
# ################################################################################
# records1 = file1.readlines()
# file1.close()
# pts,s1,s2,s3,s4,vct1 = [[],[],[],[],[],[]]
# for line in records1:
#     pts.append([float(line.split(" ")[0]),float(line.split(" ")[1]),float(line.split(" ")[2])])
#     vct1.append([float(line.split(" ")[6]),float(line.split(" ")[7]),float(line.split(" ")[8])])
#     lambda3 = float(line.split(" ")[3])
#     lambda2 = float(line.split(" ")[4])
#     lambda1 = float(line.split(" ")[5])  #最大特征值
#     s1.append(float(line.split(" ")[11]))
#     # s1.append(lambda1 - lambda2)
#     s2.append(float(line.split(" ")[9]))
#     # s3.append(float(line.split(" ")[5]))
#     s3.append(float(line.split(" ")[12]))
#     s4.append(lambda1*lambda2*lambda3)

# s = np.array(s1)
# ss = np.array(s2)
# maxeg = np.array(s3)
# vct=np.array(vct1)
# det = np.array(s4)

# sgrid.points = np.array(pts)
# sgrid.point_data.scalars = np.ravel(s.copy())
# sgrid.point_data.scalars.name = 'Cl'
# sgrid.point_data.vectors = vct
# sgrid.point_data.vectors.name = 'vector'

# sgridb.points = np.array(pts)
# sgridb.point_data.scalars = np.ravel(ss.copy())
# sgridb.point_data.scalars.name = 'beta'

# sgride.points = np.array(pts)
# sgride.point_data.scalars = np.ravel(maxeg.copy())
# sgride.point_data.scalars.name = 'MaxEg'

# sgrid_det.points = np.array(pts)
# sgrid_det.point_data.scalars = np.ravel(det.copy())
# sgrid_det.point_data.scalars.name = 'DetQ'

# # Uncomment the next two lines to save the dataset to a VTK XML file.
# w = tvtk.XMLStructuredGridWriter(input=sgridb, file_name='beta.vts')
# w.write()

'''
w = tvtk.XMLStructuredGridWriter(input=sgrid, file_name='Cl.vts')
w.write()
'''

'''
dims = (128, 201, 33)
# dims = (128, 201, 64)

sgrid = tvtk.StructuredGrid(dimensions=dims)
sgridb = tvtk.StructuredGrid(dimensions=dims)
sgride = tvtk.StructuredGrid(dimensions=dims)
sgrid_det = tvtk.StructuredGrid(dimensions=dims)
################################################################################ 文件名
file1 = open(r"polyvertical_t_-1.79_R_24.52_D_2.20_I_201_J_128_K_32_eta_1039.3_L21_+0.0_s0_1.00_test_D_2.2_AD_1_point.txt","r")
################################################################################
records1 = file1.readlines()
file1.close()
pts,s1,s2,s3,s4,vct1 = [[],[],[],[],[],[]]
for line in records1:
    pts.append([float(line.split(" ")[0]),float(line.split(" ")[1]),float(line.split(" ")[2])])
    vct1.append([float(line.split(" ")[6]),float(line.split(" ")[7]),float(line.split(" ")[8])])
    lambda3 = float(line.split(" ")[3])
    lambda2 = float(line.split(" ")[4])
    lambda1 = float(line.split(" ")[5])  #最大特征值
    s1.append(float(line.split(" ")[11]))
    # s1.append(lambda1 - lambda2)
    s2.append(float(line.split(" ")[9]))
    # s3.append(float(line.split(" ")[5]))
    s3.append(float(line.split(" ")[12]))
    s4.append(lambda1*lambda2*lambda3)

s = np.array(s1)
ss = np.array(s2)
maxeg = np.array(s3)
vct=np.array(vct1)
det = np.array(s4)

sgrid.points = np.array(pts)
sgrid.point_data.scalars = np.ravel(s.copy())
sgrid.point_data.scalars.name = 'Cl'
sgrid.point_data.vectors = vct
sgrid.point_data.vectors.name = 'vector'

sgridb.points = np.array(pts)
sgridb.point_data.scalars = np.ravel(ss.copy())
sgridb.point_data.scalars.name = 'beta'

sgride.points = np.array(pts)
sgride.point_data.scalars = np.ravel(maxeg.copy())
sgride.point_data.scalars.name = 'MaxEg'

sgrid_det.points = np.array(pts)
sgrid_det.point_data.scalars = np.ravel(det.copy())
sgrid_det.point_data.scalars.name = 'DetQ'

# Uncomment the next two lines to save the dataset to a VTK XML file.
w = tvtk.XMLStructuredGridWriter(input=sgrid, file_name='Cl.vts')
w.write()
'''

dims = (32, 128, 1)
sgrid_v = tvtk.StructuredGrid(dimensions=dims)

################################################################################文件名
file1 = open("data/kpParam_2_2_0_0_6_1_2_0.5_0.5","r")
################################################################################
records1 = file1.readlines()
file1.close()
vct2,pts2,tensor = [[],[],[]]
i = 0;
for line in records1:
    pts2.append([i/128,i%128,0])
    i = i+1;
    #pts2.append([float(line.split(" ")[0]),float(line.split(" ")[1]),float(line.split(" ")[2])])
    linesp= line.split()
    # vct2.append([float(line.split()[6]),float(line.split()[7]),float(line.split()[8])])
    tensor.append([1.0/3 + float(linesp[1]),float(linesp[4]),float(linesp[5]),float(linesp[4]),1.0/3 + float(linesp[2]),float(linesp[6]),float(linesp[5]),float(linesp[6]), 1.0/3 +float(linesp[3])])


#vector=np.array(vct2)

sgrid_v.points = np.array(pts2)
#sgrid_v.point_data.vectors = vector
#sgrid_v.point_data.vectors.name = 'vector'
sgrid_v.point_data.tensors = tensor
sgrid_v.point_data.tensors.name = 'Q'

w = tvtk.XMLStructuredGridWriter(input=sgrid_v, file_name='script/Q1.vts')
w.write()





