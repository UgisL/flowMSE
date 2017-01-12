#!/usr/bin/python2.7
# ----------------------------------------------- COPYRIGHT --------------------------------------
# Copyright 2016-2017
# Ugis Lacis, ugis.lacis@gmail.com
# Shervin Bagheri, shervin.bagheri@mech.kth.se
# -------------------------------------------- LICENSE LGPLv3 ------------------------------------
# This file is part of Poroelastic_full_strCont.
#
# Poroelastic_full_strCont is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Poroelastic_full_strCont is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Poroelastic_full_strCont. If not, see <http://www.gnu.org/licenses/>.
# ---------------------------------------------- DESCRIPTION -------------------------------------
#
# This software is generating the geometry definition for Gmsh mesher. It defines all necessary
# physical domains for the FreeFem++ solver. It is meant to be executed using Python. It has been
# tested to produce a valid Gmsh geometry file using Python version 2.7, numpy version 1.6.2 and
# scipy version 0.11.0.
#
# Execute this file by typing "./3D_MSE_gen_geom.py" in command line (Unix). Make sure that
# the top line in this file includes valid path to Python and that case file
# <3D_MSE_parameters.in> exists. For more information, consult documentation in <doc/>  

import sys
import numpy as np
import scipy.optimize as opt
import shutil as sh

# Read the case parameters
try:
    InputParam = np.loadtxt("3D_MSE_parameters.in")
except IOError:
    sys.exit("Error: expecting parameter file <3D_MSE_parameters.in>")
ztop   =     InputParam[0]
Nbot   = int(InputParam[1])
eA     =     InputParam[2] # x axis
eB     =     InputParam[3] # y axis
eC     =     InputParam[4] # z axis
cR     =     InputParam[5]
ePhi   =     InputParam[6]
resM   =     InputParam[7]
resM2  =     InputParam[8]
# Display used parameters
print "Finished reading parameters, ztop = "+str(ztop)+", Nbot = "+str(Nbot)

# Auxiliary variable, intersection distance from the center of the unit cell
zmin = - Nbot

# Number of B-spline points for 1/4th of curved circle
Nbsp = 6
# Rescale angle to radians, sign for correct turn direction
ePhi = -ePhi/180.0*np.pi


# -----------------------------------------------------------------------------------------------------
# Construct interior unit cell, below the interface cell
# -----------------------------------------------------------------------------------------------------
GeoFile  = '3D_MSE_gen_geomI.geo'
GeoFileS = '3D_MSE_gen_geomIs.geo'
# Determine the z coordinates
zmin = - Nbot - 1.0
zmax = - Nbot

# Compute the count of points {center ellipse sup} + {border cyl} + {bspline cyl} + {mid-ellipse points} + {cell-corners} + {add ell support}
npnt  = 1+3+2+2+2 + 6*9 + 6*4*Nbsp + 3*4+2*4 + 8 + 2+2+2+2
Coord = np.zeros((npnt,4))
# Write the smallest mesh spacing by default
for i in range(npnt):
    Coord[i,3] = resM;

# Spacing for ellipse segment support points
dsub = 0.025
# Center of the ellipsoid (0,0,0), unmodified
# Aligned coordinate system
Coord[1,0] = dsub; Coord[2,1] = dsub; Coord[3,2] = dsub;
# Rotated coordinate system, ellipseoid rotation
Coord[4,0] = dsub*np.cos(ePhi);    Coord[4,2] =-dsub*np.sin(ePhi)
Coord[5,0] = dsub*np.sin(ePhi);    Coord[5,2] = dsub*np.cos(ePhi)
# Rotated coordinate system, 45 degrees or pi/4 radians around y-axis
Coord[6,0] = dsub*np.cos(np.pi/4); Coord[6,2] = dsub*np.sin(np.pi/4)
Coord[7,0] =-dsub*np.sin(np.pi/4); Coord[7,2] = dsub*np.cos(np.pi/4)
# Rotated coordinate system, 45 degrees or pi/4 radians around z-axis
Coord[8,0] = dsub*np.cos(np.pi/4); Coord[8,1] = dsub*np.sin(np.pi/4)
Coord[9,0] =-dsub*np.sin(np.pi/4); Coord[9,1] = dsub*np.cos(np.pi/4)
# For planes z = +- y
nTES  = npnt - 8
cs45  = np.cos(np.pi/4)
ePhi2 = ePhi
Coord[nTES+0,0] = dsub*np.cos(ePhi2); Coord[nTES+0,1] =-dsub*np.sin(ePhi2)*cs45 ; Coord[nTES+0,2] =-dsub*np.sin(ePhi2)*cs45    # hat{x} coordinates after 2nd rotation
Coord[nTES+1,0] = dsub*np.sin(ePhi2); Coord[nTES+1,1] = dsub*np.cos(ePhi2)*cs45;  Coord[nTES+1,2] = dsub*np.cos(ePhi2)*cs45    # hat{z} coordinates after 2nd rotation
Coord[nTES+2,0] = dsub*np.cos(ePhi2); Coord[nTES+2,1] = dsub*np.sin(ePhi2)*cs45 ; Coord[nTES+2,2] =-dsub*np.sin(ePhi2)*cs45    # hat{x} coordinates after 2nd rotation
Coord[nTES+3,0] = dsub*np.sin(ePhi2); Coord[nTES+3,1] =-dsub*np.cos(ePhi2)*cs45;  Coord[nTES+3,2] = dsub*np.cos(ePhi2)*cs45    # hat{z} coordinates after 2nd rotation
# For planes y = +- x
Coord[nTES+4,0] = dsub*np.cos(ePhi2)*cs45; Coord[nTES+4,1] = dsub*np.cos(ePhi2)*cs45 ; Coord[nTES+4,2] =-dsub*np.sin(ePhi2)    # hat{x} coordinates after 2nd rotation
Coord[nTES+5,0] = dsub*np.sin(ePhi2)*cs45; Coord[nTES+5,1] = dsub*np.sin(ePhi2)*cs45;  Coord[nTES+5,2] = dsub*np.cos(ePhi2)    # hat{z} coordinates after 2nd rotation
Coord[nTES+6,0] = dsub*np.cos(ePhi2)*cs45; Coord[nTES+6,1] =-dsub*np.cos(ePhi2)*cs45 ; Coord[nTES+6,2] =-dsub*np.sin(ePhi2)    # hat{x} coordinates after 2nd rotation
Coord[nTES+7,0] = dsub*np.sin(ePhi2)*cs45; Coord[nTES+7,1] =-dsub*np.sin(ePhi2)*cs45;  Coord[nTES+7,2] = dsub*np.cos(ePhi2)    # hat{z} coordinates after 2nd rotation

npntOffs = 10
# Generate support points for circular curved segments (modify only non-zero components)
Coord[npntOffs+0,0] = -0.5; Coord[npntOffs+1,0] = 0.5;
Coord[npntOffs+2,1] = -0.5; Coord[npntOffs+3,1] = 0.5;
Coord[npntOffs+4,2] = -0.5; Coord[npntOffs+5,2] = 0.5;

# Generate points for side circles, first around x axis, then y, then z, first - number offset
npntOffs = npntOffs+6
for i in range(2):
    for j in range(8):
        Coord[npntOffs+i*8+j,0] = Coord[i+npntOffs-6,0];
        Coord[npntOffs+i*8+j,1] = cR*np.cos(j*np.pi/4);
        Coord[npntOffs+i*8+j,2] = cR*np.sin(j*np.pi/4);
npntOffs = npntOffs+16
for i in range(2):
    for j in range(8):
        Coord[npntOffs+i*8+j,0] = cR*np.cos(j*np.pi/4);
        Coord[npntOffs+i*8+j,1] = Coord[i+npntOffs-20,1];
        Coord[npntOffs+i*8+j,2] = cR*np.sin(j*np.pi/4);
npntOffs = npntOffs+16
for i in range(2):
    for j in range(8):
        Coord[npntOffs+i*8+j,0] = cR*np.cos(j*np.pi/4);
        Coord[npntOffs+i*8+j,1] = cR*np.sin(j*np.pi/4);
        Coord[npntOffs+i*8+j,2] = Coord[i+npntOffs-34,2];
        
# Generate points for inner circles, using B-splines
npntOffs = npntOffs+16
# Constant auxiliary variables
sinP2 = np.sin(ePhi)**2
cosP2 = np.cos(ePhi)**2
scP   = np.sin(ePhi)*np.cos(ePhi)
# Ax = cos^2/a^2+sin^2/c^2; Bx = 2*z*sin*cos/c^2-2*z*sin*cos/a^2; Cx = z^2*sin^2/a^2+y^2/b^2+z^2*cos^2/c^2-1; y = cR*sin(t); z = cR*cos(t)
for i in range(4*Nbsp):
    tp            = 2.0*np.pi*i/(4.0*Nbsp)
    ycur          = cR*np.sin(tp)
    zcur          = cR*np.cos(tp)
    Ax            = cosP2/eA**2 + sinP2/eC**2
    Bx            = 2*zcur*scP*(1.0/eC**2 - 1.0/eA**2)
    Cx            = zcur**2*(sinP2/eA**2 + cosP2/eC**2) + ycur**2/eB**2-1
    Coord[npntOffs+i,0] = (-Bx + np.sqrt(Bx**2-4.0*Ax*Cx) )/(2.0*Ax)
    Coord[npntOffs+i,1] = ycur
    Coord[npntOffs+i,2] = zcur
for i in range(4*Nbsp):
    tp                   = 2.0*np.pi*i/(4.0*Nbsp)
    ycur                 = cR*np.sin(tp)
    zcur                 = cR*np.cos(tp)
    Ax                   = cosP2/eA**2 + sinP2/eC**2
    Bx                   = 2*zcur*scP*(1.0/eC**2 - 1.0/eA**2)
    Cx                   = zcur**2*(sinP2/eA**2 + cosP2/eC**2) + ycur**2/eB**2-1
    Coord[npntOffs+4*Nbsp+i,0] = (-Bx - np.sqrt(Bx**2-4.0*Ax*Cx) )/(2.0*Ax)
    Coord[npntOffs+4*Nbsp+i,1] = ycur
    Coord[npntOffs+4*Nbsp+i,2] = zcur
# Ay = 1/b^2; By = 0; Cy = (x^2*cos^2+z^2*sin^2-2*x*z*sin*cos)/a^2+(x^2*sin^2+z^2*cos^2+2*x*z*sin*cos)/c^2-1; x = cR*sin(t); z = cR*cos(t)
for i in range(4*Nbsp):
    tp                   = 2.0*np.pi*i/(4.0*Nbsp)
    xcur                 = cR*np.sin(tp)
    zcur                 = cR*np.cos(tp)
    Ay                   = 1.0/eB**2
    By                   = 0
    Cy                   = (xcur**2*cosP2 + zcur**2*sinP2 - 2*xcur*zcur*scP)/eA**2  \
                         + (xcur**2*sinP2 + zcur**2*cosP2 + 2*xcur*zcur*scP)/eC**2 - 1.0
    Coord[npntOffs+8*Nbsp+i,0] = xcur
    Coord[npntOffs+8*Nbsp+i,1] = np.sqrt(By**2-4.0*Ay*Cy)/(2.0*Ay)
    Coord[npntOffs+8*Nbsp+i,2] = zcur
for i in range(4*Nbsp):
    tp                    = 2.0*np.pi*i/(4.0*Nbsp)
    xcur                  = cR*np.sin(tp)
    zcur                  = cR*np.cos(tp)
    Ay                    = 1.0/eB**2
    By                    = 0
    Cy                    = (xcur**2*cosP2 + zcur**2*sinP2 - 2*xcur*zcur*scP)/eA**2  \
                          + (xcur**2*sinP2 + zcur**2*cosP2 + 2*xcur*zcur*scP)/eC**2 - 1.0
    Coord[npntOffs+12*Nbsp+i,0] = xcur
    Coord[npntOffs+12*Nbsp+i,1] =-np.sqrt(By**2-4.0*Ay*Cy)/(2.0*Ay)
    Coord[npntOffs+12*Nbsp+i,2] = zcur
# Az = sin^2/a^2+cos^2/c^2; Bz = 2*x*sin*cos/c^2-2*x*sin*cos/a^2; Cz = x^2*cos^2/a^2+y^2/b^2+x^2*sin^2/c^2-1; x = cR*sin(t); y = cR*cos(t)
for i in range(4*Nbsp):
    tp                    = 2.0*np.pi*i/(4.0*Nbsp)
    xcur                  = cR*np.sin(tp)
    ycur                  = cR*np.cos(tp)
    Az                    = sinP2/eA**2 + cosP2/eC**2
    Bz                    = 2*xcur*scP*(1.0/eC**2 - 1.0/eA**2)
    Cz                    = xcur**2*(cosP2/eA**2 + sinP2/eC**2) + ycur**2/eB**2 - 1.0
    Coord[npntOffs+16*Nbsp+i,0] = xcur
    Coord[npntOffs+16*Nbsp+i,1] = ycur
    Coord[npntOffs+16*Nbsp+i,2] = (-Bz + np.sqrt(Bz**2-4.0*Az*Cz) )/(2.0*Az)
for i in range(4*Nbsp):
    tp                    = 2.0*np.pi*i/(4.0*Nbsp)
    xcur                  = cR*np.sin(tp)
    ycur                  = cR*np.cos(tp)
    Az                    = sinP2/eA**2 + cosP2/eC**2
    Bz                    = 2*xcur*scP*(1.0/eC**2 - 1.0/eA**2)
    Cz                    = xcur**2*(cosP2/eA**2 + sinP2/eC**2) + ycur**2/eB**2 - 1.0
    Coord[npntOffs+20*Nbsp+i,0] = xcur
    Coord[npntOffs+20*Nbsp+i,1] = ycur
    Coord[npntOffs+20*Nbsp+i,2] = (-Bz - np.sqrt(Bz**2-4.0*Az*Cz) )/(2.0*Az)

# Mark half points of the connecting aligned arcs
npntOffs = npntOffs+6*4*Nbsp
# z = 0
Ax  = np.sqrt(1.0/(cosP2/eA**2+sinP2/eC**2))
tpm = np.arctan(eB/Ax)
Coord[npntOffs+0,0] = Ax*np.sin(tpm); Coord[npntOffs+0,1] = eB*np.cos(tpm); Coord[npntOffs+0,2] = 0
Coord[npntOffs+1,0] =-Ax*np.sin(tpm); Coord[npntOffs+1,1] = eB*np.cos(tpm); Coord[npntOffs+1,2] = 0
Coord[npntOffs+2,0] = Ax*np.sin(tpm); Coord[npntOffs+2,1] =-eB*np.cos(tpm); Coord[npntOffs+2,2] = 0
Coord[npntOffs+3,0] =-Ax*np.sin(tpm); Coord[npntOffs+3,1] =-eB*np.cos(tpm); Coord[npntOffs+3,2] = 0
XYi = Ax*np.sin(tpm)

# y = 0, x = +- z
npntOffs = npntOffs+4
xzcur = np.sqrt(1.0/( (np.cos(ePhi)-np.sin(ePhi))**2/eA**2 + (np.sin(ePhi)+np.cos(ePhi))**2/eC**2 ))
Coord[npntOffs+0,0] = xzcur; Coord[npntOffs+0,1] = 0; Coord[npntOffs+0,2] = xzcur
Coord[npntOffs+1,0] =-xzcur; Coord[npntOffs+1,1] = 0; Coord[npntOffs+1,2] =-xzcur
xzcur = np.sqrt(1.0/( (np.cos(ePhi)+np.sin(ePhi))**2/eA**2 + (np.sin(ePhi)-np.cos(ePhi))**2/eC**2 ))
Coord[npntOffs+2,0] = xzcur; Coord[npntOffs+2,1] = 0; Coord[npntOffs+2,2] =-xzcur
Coord[npntOffs+3,0] =-xzcur; Coord[npntOffs+3,1] = 0; Coord[npntOffs+3,2] = xzcur

# x = 0
npntOffs = npntOffs+4
Az  = np.sqrt(1.0/(sinP2/eA**2+cosP2/eC**2))
tpm = np.arctan(eB/Az);
Coord[npntOffs+0,0] = 0; Coord[npntOffs+0,1] = eB*np.cos(tpm); Coord[npntOffs+0,2] = Az*np.sin(tpm)
Coord[npntOffs+1,0] = 0; Coord[npntOffs+1,1] =-eB*np.cos(tpm); Coord[npntOffs+1,2] = Az*np.sin(tpm)
Coord[npntOffs+2,0] = 0; Coord[npntOffs+2,1] = eB*np.cos(tpm); Coord[npntOffs+2,2] =-Az*np.sin(tpm)
Coord[npntOffs+3,0] = 0; Coord[npntOffs+3,1] =-eB*np.cos(tpm); Coord[npntOffs+3,2] =-Az*np.sin(tpm)
YZi = Az*np.sin(tpm)

# x = y = z
npntOffs = npntOffs+4
curxyz = np.sqrt( 1.0/( (np.cos(ePhi)-np.sin(ePhi))**2/eA**2 + 1.0/eB**2 \
                      + (np.sin(ePhi)+np.cos(ePhi))**2/eC**2 ) )
Coord[npntOffs+0,0] = curxyz; Coord[npntOffs+0,1] = curxyz; Coord[npntOffs+0,2] = curxyz
Coord[npntOffs+2,0] = curxyz; Coord[npntOffs+2,1] =-curxyz; Coord[npntOffs+2,2] = curxyz
Coord[npntOffs+5,0] =-curxyz; Coord[npntOffs+5,1] = curxyz; Coord[npntOffs+5,2] =-curxyz
Coord[npntOffs+7,0] =-curxyz; Coord[npntOffs+7,1] =-curxyz; Coord[npntOffs+7,2] =-curxyz
curxyz = np.sqrt( 1.0/( (np.cos(ePhi)+np.sin(ePhi))**2/eA**2 + 1.0/eB**2 \
                      + (np.sin(ePhi)-np.cos(ePhi))**2/eC**2 ) )
Coord[npntOffs+1,0] =-curxyz; Coord[npntOffs+1,1] = curxyz; Coord[npntOffs+1,2] = curxyz
Coord[npntOffs+3,0] =-curxyz; Coord[npntOffs+3,1] =-curxyz; Coord[npntOffs+3,2] = curxyz
Coord[npntOffs+4,0] = curxyz; Coord[npntOffs+4,1] = curxyz; Coord[npntOffs+4,2] =-curxyz
Coord[npntOffs+6,0] = curxyz; Coord[npntOffs+6,1] =-curxyz; Coord[npntOffs+6,2] =-curxyz

# Additional points at the corners of the unit cell
npntOffs = npntOffs+8
for i in range(2):
    for j in range(2):
        for k in range(2):
            Coord[npntOffs+i*4+j*2+k,0] = -0.5 + i;
            Coord[npntOffs+i*4+j*2+k,1] = -0.5 + j;
            Coord[npntOffs+i*4+j*2+k,2] = -0.5 + k;
            Coord[npntOffs+i*4+j*2+k,3] = resM2;

# Shifting all the points, to position mesh between zmin and zmax
for i in range(npnt):
    Coord[i,2] = Coord[i,2] + 0.5 + zmin

# Exporting all the points
fidF = open(GeoFile,'w')
for i in range(npnt):
    fidF.write('Point(%d) = {%.10f, %.10f, %.10f, %.10f};\n'%(i+1,Coord[i,0],Coord[i,1],Coord[i,2],Coord[i,3]))

# -------------------------------------------------------------
# Below here, the operations can be copied manually from GUI
# Constructions from loops also allowed, if feasable
# -------------------------------------------------------------
# Generate the connecting circle arcs around x, y and z axis
fidF.write('\
Circle(1) = {33, 13, 34};\n\
Circle(2) = {34, 13, 35};\n\
Circle(3) = {35, 13, 36};\n\
Circle(4) = {36, 13, 37};\n\
Circle(5) = {37, 13, 38};\n\
Circle(6) = {38, 13, 39};\n\
Circle(7) = {39, 13, 40};\n\
Circle(8) = {40, 13, 33};\n\
Circle(9) = {41, 14, 42};\n\
Circle(10) = {42, 14, 43};\n\
Circle(11) = {43, 14, 44};\n\
Circle(12) = {44, 14, 45};\n\
Circle(13) = {45, 14, 46};\n\
Circle(14) = {46, 14, 47};\n\
Circle(15) = {47, 14, 48};\n\
Circle(16) = {48, 14, 41};\n\
Circle(17) = {25, 12, 32};\n\
Circle(18) = {32, 12, 31};\n\
Circle(19) = {31, 12, 30};\n\
Circle(20) = {30, 12, 29};\n\
Circle(21) = {29, 12, 28};\n\
Circle(22) = {28, 12, 27};\n\
Circle(23) = {27, 12, 26};\n\
Circle(24) = {26, 12, 25};\n\
Circle(25) = {21, 11, 22};\n\
Circle(26) = {22, 11, 23};\n\
Circle(27) = {23, 11, 24};\n\
Circle(28) = {24, 11, 17};\n\
Circle(29) = {17, 11, 18};\n\
Circle(30) = {18, 11, 19};\n\
Circle(31) = {19, 11, 20};\n\
Circle(32) = {20, 11, 21};\n\
Circle(33) = {54, 15, 53};\n\
Circle(34) = {53, 15, 52};\n\
Circle(35) = {52, 15, 51};\n\
Circle(36) = {51, 15, 50};\n\
Circle(37) = {50, 15, 49};\n\
Circle(38) = {49, 15, 56};\n\
Circle(39) = {56, 15, 55};\n\
Circle(40) = {55, 15, 54};\n\
Circle(41) = {63, 16, 62};\n\
Circle(42) = {62, 16, 61};\n\
Circle(43) = {61, 16, 60};\n\
Circle(44) = {60, 16, 59};\n\
Circle(45) = {59, 16, 58};\n\
Circle(46) = {58, 16, 57};\n\
Circle(47) = {57, 16, 64};\n\
Circle(48) = {64, 16, 63};\n\
')

# Write Bsplines of the cylinder intersection with ellipsoid
fidF.write('\
BSpline(49) = {143, 142, 141, 140};\n\
BSpline(50) = {140, 139, 138, 137};\n\
BSpline(51) = {137, 160, 159, 158};\n\
BSpline(52) = {158, 157, 156, 155};\n\
BSpline(53) = {155, 154, 153, 152};\n\
BSpline(54) = {152, 151, 150, 149};\n\
BSpline(55) = {149, 148, 147, 146};\n\
BSpline(56) = {146, 145, 144, 143};\n\
BSpline(57) = {128, 129, 130, 131};\n\
BSpline(58) = {131, 132, 133, 134};\n\
BSpline(59) = {134, 135, 136, 113};\n\
BSpline(60) = {113, 114, 115, 116};\n\
BSpline(61) = {116, 117, 118, 119};\n\
BSpline(62) = {119, 120, 121, 122};\n\
BSpline(63) = {122, 123, 124, 125};\n\
BSpline(64) = {125, 126, 127, 128};\n\
BSpline(65) = {107, 108, 109, 110};\n\
BSpline(66) = {110, 110, 111, 112, 89};\n\
BSpline(67) = {89, 90, 91, 92};\n\
BSpline(68) = {92, 93, 94, 95};\n\
BSpline(69) = {95, 96, 97, 98};\n\
BSpline(70) = {98, 99, 100, 101};\n\
BSpline(71) = {101, 102, 103, 104};\n\
BSpline(72) = {104, 105, 106, 107};\n\
BSpline(73) = {71, 70, 69, 68};\n\
BSpline(74) = {68, 67, 66, 65};\n\
BSpline(75) = {65, 88, 87, 86};\n\
BSpline(76) = {86, 85, 84, 83};\n\
BSpline(77) = {83, 82, 81, 80};\n\
BSpline(78) = {80, 79, 78, 77};\n\
BSpline(79) = {77, 76, 75, 74};\n\
BSpline(80) = {74, 73, 72, 71};\n\
BSpline(81) = {197, 198, 199, 200};\n\
BSpline(82) = {200, 201, 202, 203};\n\
BSpline(83) = {203, 204, 205, 206};\n\
BSpline(84) = {206, 207, 208, 185};\n\
BSpline(85) = {185, 186, 187, 188};\n\
BSpline(86) = {188, 189, 190, 191};\n\
BSpline(87) = {191, 192, 193, 194};\n\
BSpline(88) = {194, 195, 196, 197};\n\
BSpline(89) = {161, 184, 183, 182};\n\
BSpline(90) = {182, 181, 180, 179};\n\
BSpline(91) = {179, 178, 177, 176};\n\
BSpline(92) = {176, 175, 174, 173};\n\
BSpline(93) = {173, 172, 171, 170};\n\
BSpline(94) = {170, 169, 168, 167};\n\
BSpline(95) = {167, 166, 165, 164};\n\
BSpline(96) = {164, 163, 162, 161};\n\
')

# Write connecting lines for the side cylinders
fidF.write('\
Line(97) = {21, 107};\n\
Line(98) = {20, 110};\n\
Line(99) = {19, 89};\n\
Line(100) = {18, 92};\n\
Line(101) = {17, 95};\n\
Line(102) = {24, 98};\n\
Line(103) = {23, 101};\n\
Line(104) = {22, 104};\n\
Line(105) = {27, 65};\n\
Line(106) = {26, 68};\n\
Line(107) = {25, 71};\n\
Line(108) = {32, 74};\n\
Line(109) = {31, 77};\n\
Line(110) = {30, 80};\n\
Line(111) = {29, 83};\n\
Line(112) = {28, 86};\n\
Line(113) = {44, 134};\n\
Line(114) = {45, 131};\n\
Line(115) = {46, 128};\n\
Line(116) = {47, 125};\n\
Line(117) = {48, 122};\n\
Line(118) = {41, 119};\n\
Line(119) = {116, 42};\n\
Line(120) = {43, 113};\n\
Line(121) = {35, 137};\n\
Line(122) = {34, 140};\n\
Line(123) = {33, 143};\n\
Line(124) = {40, 146};\n\
Line(125) = {39, 149};\n\
Line(126) = {38, 152};\n\
Line(127) = {37, 155};\n\
Line(128) = {36, 158};\n\
Line(129) = {61, 179};\n\
Line(130) = {60, 182};\n\
Line(131) = {59, 161};\n\
Line(132) = {58, 164};\n\
Line(133) = {57, 167};\n\
Line(134) = {64, 170};\n\
Line(135) = {63, 173};\n\
Line(136) = {62, 176};\n\
Line(137) = {53, 203};\n\
Line(138) = {54, 200};\n\
Line(139) = {55, 197};\n\
Line(140) = {56, 194};\n\
Line(141) = {49, 191};\n\
Line(142) = {50, 188};\n\
Line(143) = {51, 185};\n\
Line(144) = {52, 206};\n\
')

# Write the ellipse arcs around the ellipsoid
# Plane z = 0
fidF.write('\
Ellipse(145) = {143, 1, 2, 211};\n\
Ellipse(146) = {211, 1, 2, 83};\n\
Ellipse(147) = {71, 1, 2, 209};\n\
Ellipse(148) = {209, 1, 2, 119};\n\
Ellipse(149) = {131, 1, 2, 210};\n\
Ellipse(150) = {210, 1, 2, 95};\n\
Ellipse(151) = {107, 1, 2, 212};\n\
Ellipse(152) = {212, 1, 2, 155};\n\
')
# Plane y = 0
fidF.write('\
Ellipse(153) = {65, 1, 5, 213};\n\
Ellipse(154) = {213, 1, 5, 167};\n\
Ellipse(155) = {179, 1, 5, 216};\n\
Ellipse(156) = {216, 1, 5, 89};\n\
Ellipse(157) = {101, 1, 5, 214};\n\
Ellipse(158) = {214, 1, 5, 203};\n\
Ellipse(159) = {191, 1, 5, 215};\n\
Ellipse(160) = {215, 1, 5, 77};\n\
')
# Plane x = 0
fidF.write('\
Ellipse(161) = {137, 1, 3, 218};\n\
Ellipse(162) = {218, 1, 3, 173};\n\
Ellipse(163) = {161, 1, 3, 217};\n\
Ellipse(164) = {217, 1, 3, 113};\n\
Ellipse(165) = {125, 1, 3, 219};\n\
Ellipse(166) = {219, 1, 3, 185};\n\
Ellipse(167) = {197, 1, 3, 220};\n\
Ellipse(168) = {220, 1, 3, 149};\n\
')

# Plane z = x
fidF.write('\
Ellipse(169) = {140, 1, 7, 223};\n\
Ellipse(170) = {223, 1, 7, 213};\n\
Ellipse(171) = {213, 1, 7, 221};\n\
Ellipse(172) = {221, 1, 7, 116};\n\
Ellipse(173) = {128, 1, 7, 226};\n\
Ellipse(174) = {226, 1, 7, 214};\n\
Ellipse(175) = {214, 1, 7, 228};\n\
Ellipse(176) = {228, 1, 7, 152};\n\
')

# Plane z = -x
fidF.write('\
Ellipse(177) = {158, 1, 3, 224};\n\
Ellipse(178) = {224, 1, 3, 216};\n\
Ellipse(179) = {216, 1, 3, 222};\n\
Ellipse(180) = {222, 1, 3, 134};\n\
Ellipse(181) = {122, 1, 3, 225};\n\
Ellipse(182) = {225, 1, 3, 215};\n\
Ellipse(183) = {215, 1, 3, 227};\n\
Ellipse(184) = {227, 1, 3, 146};\n\
')

# Plane y = x
maP = nTES+5
fidF.write('\
Ellipse(185) = {164, 1, %d, 221};\n\
Ellipse(186) = {221, 1, %d, 209};\n\
Ellipse(187) = {209, 1, %d, 225};\n\
Ellipse(188) = {225, 1, %d, 188};\n\
Ellipse(189) = {200, 1, %d, 228};\n\
Ellipse(190) = {228, 1, %d, 212};\n\
Ellipse(191) = {212, 1, %d, 224};\n\
Ellipse(192) = {224, 1, %d, 176};\n\
'%(maP,maP,maP,maP,maP,maP,maP,maP))

# Plane y = -x
maP = nTES+7
fidF.write('\
Ellipse(193) = {182, 1, %d, 222};\n\
Ellipse(194) = {222, 1, %d, 210};\n\
Ellipse(195) = {210, 1, %d, 226};\n\
Ellipse(196) = {226, 1, %d, 206};\n\
Ellipse(197) = {194, 1, %d, 227};\n\
Ellipse(198) = {227, 1, %d, 211};\n\
Ellipse(199) = {211, 1, %d, 223};\n\
Ellipse(200) = {223, 1, %d, 170};\n\
'%(maP,maP,maP,maP,maP,maP,maP,maP))

# Plane z = y
maP = nTES+1
fidF.write('\
Ellipse(201) = {68, 1, %d, 221};\n\
Ellipse(202) = {221, 1, %d, 217};\n\
Ellipse(203) = {222, 1, %d, 92};\n\
Ellipse(204) = {217, 1, %d, 222};\n\
Ellipse(205) = {104, 1, %d, 228};\n\
Ellipse(206) = {228, 1, %d, 220};\n\
Ellipse(207) = {220, 1, %d, 227};\n\
Ellipse(208) = {227, 1, %d, 80};\n\
'%(maP,maP,maP,maP,maP,maP,maP,maP))

# Plane z = -y
maP = nTES+3
fidF.write('\
Ellipse(209) = {110, 1, %d, 224};\n\
Ellipse(210) = {224, 1, %d, 218};\n\
Ellipse(211) = {218, 1, %d, 223};\n\
Ellipse(212) = {223, 1, %d, 86};\n\
Ellipse(213) = {74, 1, %d, 225};\n\
Ellipse(214) = {225, 1, %d, 219};\n\
Ellipse(215) = {219, 1, %d, 226};\n\
Ellipse(216) = {226, 1, %d, 98};\n\
'%(maP,maP,maP,maP,maP,maP,maP,maP))


# Lines, denoting the outer boundary of the cell
fidF.write('\
Line(217) = {234, 236};\n\
Line(218) = {236, 232};\n\
Line(219) = {232, 230};\n\
Line(220) = {230, 234};\n\
Line(221) = {234, 233};\n\
Line(222) = {236, 235};\n\
Line(223) = {232, 231};\n\
Line(224) = {231, 229};\n\
Line(225) = {230, 229};\n\
Line(226) = {229, 233};\n\
Line(227) = {233, 235};\n\
Line(228) = {235, 231};\n\
')

# Add interior surfaces
fidF.write('\
Line Loop(230) = {21, 112, 76, -111};\n\
Line Loop(232) = {22, 105, 75, -112};\n\
Line Loop(234) = {23, 106, 74, -105};\n\
Line Loop(236) = {106, -73, -107, -24};\n\
Line Loop(238) = {80, -107, 17, 108};\n\
Line Loop(240) = {108, -79, -109, -18};\n\
Line Loop(242) = {109, -78, -110, -19};\n\
Line Loop(244) = {110, -77, -111, -20};\n\
Line Loop(246) = {66, -99, 31, 98};\n\
Line Loop(248) = {32, 97, 65, -98};\n\
Line Loop(250) = {97, -72, -104, -25};\n\
Line Loop(252) = {104, -71, -103, -26};\n\
Line Loop(254) = {103, -70, -102, -27};\n\
Line Loop(256) = {28, 101, 69, -102};\n\
Line Loop(258) = {29, 100, 68, -101};\n\
Line Loop(260) = {30, 99, 67, -100};\n\
Line Loop(262) = {4, 127, -52, -128};\n\
Line Loop(264) = {5, 126, -53, -127};\n\
Line Loop(266) = {126, 54, -125, -6};\n\
Line Loop(268) = {125, 55, -124, -7};\n\
Line Loop(270) = {124, 56, -123, -8};\n\
Line Loop(272) = {1, 122, -49, -123};\n\
Line Loop(274) = {122, 50, -121, -2};\n\
Line Loop(276) = {3, 128, -51, -121};\n\
Line Loop(278) = {58, -113, 12, 114};\n\
Line Loop(280) = {59, -120, 11, 113};\n\
Line Loop(282) = {120, 60, 119, 10};\n\
Line Loop(284) = {119, -9, 118, -61};\n\
Line Loop(286) = {118, 62, -117, 16};\n\
Line Loop(288) = {15, 117, 63, -116};\n\
Line Loop(290) = {14, 116, 64, -115};\n\
Line Loop(292) = {13, 115, 57, -114};\n\
Line Loop(294) = {144, 84, -143, -35};\n\
Line Loop(296) = {143, 85, -142, -36};\n\
Line Loop(298) = {86, -141, -37, 142};\n\
Line Loop(300) = {141, 87, -140, -38};\n\
Line Loop(302) = {140, 88, -139, -39};\n\
Line Loop(304) = {139, 81, -138, -40};\n\
Line Loop(306) = {33, 137, -82, -138};\n\
Line Loop(308) = {34, 144, -83, -137};\n\
Line Loop(310) = {90, -129, 43, 130};\n\
Line Loop(312) = {130, -89, -131, -44};\n\
Line Loop(314) = {131, -96, -132, -45};\n\
Line Loop(316) = {132, -95, -133, -46};\n\
Line Loop(318) = {94, -133, 47, 134};\n\
Line Loop(320) = {134, -93, -135, -48};\n\
Line Loop(322) = {41, 136, 92, -135};\n\
Line Loop(324) = {136, -91, -129, -42};\n\
Line Loop(326) = {155, 179, -193, 90};\n\
Line Loop(328) = {193, -204, -163, 89};\n\
Line Loop(330) = {179, 203, -67, -156};\n\
Line Loop(332) = {68, -150, -194, 203};\n\
Line Loop(334) = {194, -149, 58, -180};\n\
Line Loop(336) = {180, 59, -164, 204};\n\
Line Loop(338) = {164, 60, -172, 202};\n\
Line Loop(340) = {163, -202, -185, 96};\n\
Line Loop(342) = {185, -171, 154, 95};\n\
Line Loop(344) = {153, 171, -201, 74};\n\
Line Loop(346) = {201, 186, -147, 73};\n\
Line Loop(348) = {186, 148, -61, -172};\n\
Line Loop(350) = {155, -178, 192, -91};\n\
Line Loop(352) = {156, -66, 209, 178};\n\
Line Loop(354) = {209, -191, -151, 65};\n\
Line Loop(356) = {152, -52, 177, -191};\n\
Line Loop(358) = {177, 210, -161, 51};\n\
Line Loop(360) = {210, 162, -92, -192};\n\
Line Loop(362) = {162, 93, -200, -211};\n\
Line Loop(364) = {211, -169, 50, 161};\n\
Line Loop(366) = {169, -199, -145, 49};\n\
Line Loop(368) = {170, 154, -94, -200};\n\
Line Loop(370) = {170, -153, 75, -212};\n\
Line Loop(372) = {212, 76, -146, 199};\n\
Line Loop(374) = {158, 83, -196, 174};\n\
Line Loop(376) = {174, -157, -70, -216};\n\
Line Loop(378) = {195, 216, -69, -150};\n\
Line Loop(380) = {149, 195, -173, 57};\n\
Line Loop(382) = {173, -215, -165, 64};\n\
Line Loop(384) = {215, 196, 84, -166};\n\
Line Loop(386) = {165, -214, -181, 63};\n\
Line Loop(388) = {214, 166, 85, -188};\n\
Line Loop(390) = {181, -187, 148, 62};\n\
Line Loop(392) = {187, -213, 80, 147};\n\
Line Loop(394) = {213, 182, 160, 79};\n\
Line Loop(396) = {182, -159, -86, -188};\n\
Line Loop(398) = {160, -78, -208, -183};\n\
Line Loop(400) = {183, -197, -87, 159};\n\
Line Loop(402) = {208, -77, -146, -198};\n\
Line Loop(404) = {145, -198, 184, 56};\n\
Line Loop(406) = {184, -55, -168, 207};\n\
Line Loop(408) = {207, -197, 88, 167};\n\
Line Loop(410) = {168, -54, -176, 206};\n\
Line Loop(412) = {206, -167, 81, 189};\n\
Line Loop(414) = {176, -53, -152, -190};\n\
Line Loop(416) = {151, -190, -205, 72};\n\
Line Loop(418) = {205, -175, -157, 71};\n\
Line Loop(420) = {175, -189, 82, -158};\n\
')
# Interior surface indices - from 230 to 420 (step 2), inclusive, loop
for i in range(96):
    fidF.write('Ruled Surface(%d) = {%d};\n'%(230+i*2,230+i*2))

# Add closing surfaces, fluid surfaces
fidF.write('\
Line Loop(1001) = {225, -224, -223, 219, -32, -31, -30, -29, -28, -27, -26, -25};\n\
Line Loop(1002) = {221, 227, -222, -217, -17, -24, -23, -22, -21, -20, -19, -18};\n\
Line Loop(1003) = {226, -221, -220, 225, -8, -7, -6, -5, -4, -3, -2, -1};\n\
Line Loop(1004) = {218, 223, -228, -222, -16, -15, -14, -13, -12, -11, -10, -9};\n\
Line Loop(1005) = {227, 228, 224, 226, -36, -35, -34, -33, -40, -39, -38, -37};\n\
Line Loop(1006) = {218, 219, 220, 217, -44, -43, -42, -41, -48, -47, -46, -45};\n\
')
# Define plane surfaces
for i in range(6):
    fidF.write('Plane Surface(%d) = {%d};\n'%(1001+i,1001+i))

# Add closing surfaces, solid surfaces, 11006 - will be closed later
fidF.write('\
Line Loop(11001) = {32, 25, 26, 27, 28, 29, 30, 31};\n\
Line Loop(11002) = {23, 24, 17, 18, 19, 20, 21, 22};\n\
Line Loop(11003) = {1, 2, 3, 4, 5, 6, 7, 8};\n\
Line Loop(11004) = {9, 10, 11, 12, 13, 14, 15, 16};\n\
Line Loop(11005) = {38, 39, 40, 33, 34, 35, 36, 37};\n\
Line Loop(11006) = {46, 47, 48, 41, 42, 43, 44, 45};\n\
')
# Define plane surfaces
for i in range(6):
    fidF.write('Plane Surface(%d) = {%d};\n'%(11001+i,11001+i))

# Define the volume, fluid
fidF.write('Surface Loop(1301) = {')
for i in range(6):
    fidF.write('%d, '%(1001+i))
for i in range(95):
    curInd = 230+i*2
    fidF.write('%d, '%curInd)
fidF.write('%d};\n'%(curInd+2))
fidF.write('Volume(1401) = {1301};\n')
# Define the volume, solid
fidF.write('Surface Loop(1302) = {')
fidF.write('%d, '%(11001))
fidF.write('%d, '%(11002))
fidF.write('%d, '%(11003))
fidF.write('%d, '%(11004))
fidF.write('%d, '%(11005))
fidF.write('%d, '%(11006))
for i in range(95):
    curInd = 230+i*2
    fidF.write('%d, '%curInd)
fidF.write('%d};\n'%(curInd+2))
fidF.write('Volume(1402) = {1302};\n')

# Define the periodic surfaces, direction vector to translate the right surface for matching
# In this case, left is top surface, right is bottom surface (must be raised by one unit)
fidF.write('Periodic Surface { 1002}  = { 1001} Translate {1,0,0};\n');
fidF.write('Periodic Surface {11002}  = {11001} Translate {1,0,0};\n');
fidF.write('Periodic Surface { 1004}  = { 1003} Translate {0,1,0};\n');
fidF.write('Periodic Surface {11004}  = {11003} Translate {0,1,0};\n');
fidF.write('Periodic Surface { 1006}  = { 1005} Translate {0,0,1};\n');
fidF.write('Periodic Surface {11006}  = {11005} Translate {0,0,1};\n');

# Make sure that we output only entities denoted using physical labels
fidF.write('Mesh.SaveAll = 0;\n')
# Set the mesh resolution
fidF.write('Mesh.CharacteristicLengthMax = %.4f;\n'%(resM2))

fidF.close()

# Copy the fluid file to the solid file
sh.copyfile(GeoFile,GeoFileS)

# Append physical definitions
# FLUID
fidF = open(GeoFile,'a')
# Define the physical surfaces, periodic boundaries
fidF.write('Physical Surface(3001) = {1001};\n' )
fidF.write('Physical Surface(3002) = {1002};\n' )
fidF.write('Physical Surface(3003) = {1003};\n' )
fidF.write('Physical Surface(3004) = {1004};\n' )
fidF.write('Physical Surface(3005) = {1005};\n' )
fidF.write('Physical Surface(3006) = {1006};\n' )
# Boundaries with solid
fidF.write('Physical Surface(3007) = {')
for i in range(95):
    curInd = 230+i*2
    fidF.write('%d, '%(curInd))
fidF.write('%d};\n'%(curInd+2))
# Define the physical volume
fidF.write('Physical Volume(4001) = {1401};\n\n' )
fidF.close()

# Solid
fidS = open(GeoFileS,'a')
# Define the physical points, SOLID
# (2001, 2002) - on X axis, (2002, 2003) - on Z axis
fidS.write('Physical Point(2001) = {23};\n')
fidS.write('Physical Point(2002) = {31};\n')
fidS.write('Physical Point(2003) = {27};\n')
# Define the physical surfaces, SOLID
# Periodic boundaries, including the translated ones (need to have separate
# pysical labels, because they are disconnected)
# LEFT BOUNDARY
fidS.write('Physical Surface(3001) = {11001};\n')
# RIGHT BOUNDARY
fidS.write('Physical Surface(3002) = {11002};\n')
# FRONT BOUNDARY
fidS.write('Physical Surface(3003) = {11003};\n')
# BACK BOUNDARY
fidS.write('Physical Surface(3004) = {11004};\n')
# BOTTOM BOUNDARY
fidS.write('Physical Surface(3005) = {11005};\n' )
# TOP BOUNDARY
fidS.write('Physical Surface(3006) = {11006};\n')
# SOLID BOUNDARY
fidS.write('Physical Surface(3007) = {')
for i in range(95):
    curInd = 230+i*2
    fidS.write('%d, '%(curInd))
fidS.write('%d};\n'%(curInd+2))
# Define the physical volume, SOLID
fidS.write('Physical Volume(4001) = {1402};\n\n' )
# Close the file
fidS.close()
            
# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------------------------
# Construct interface cell, bottom part
# -----------------------------------------------------------------------------------------------------
GeoFile  = '3D_MSE_gen_geomB.geo'
# Determine the z coordinates
zmin = -Nbot
zmax =  0.0

# Generate FE space definition files, fluid
fid = open('3D_MSE_gen_geomB.FEM','w')
fid.write('fespace UUUPhB (Th3B, [P2,P2,P2,P1],periodic=[\n')
fid.write('[3001 ,y,z],[3002 ,y,z],\n')
for i in range(Nbot-1):
    fid.write('[3001%d,y,z],[3002%d,y,z],\n'%(i+1,i+1))
fid.write('[3003 ,x,z],[3004 ,x,z],\n')
for i in range(Nbot-2):
    fid.write('[3003%d,x,z],[3004%d,x,z],\n'%(i+1,i+1))
fid.write('[3003%d,x,z],[3004%d,x,z]]);\n'%(i+2,i+2))
fid.close()

# Compute the count of points {center ellipse sup} + {border cyl} + {bspline cyl} + {mid-ellipse points} + {cell-corners} + {add ell support}
npnt  = 1+3+2+2+2 + 6*9 + 6*4*Nbsp + 3*4+2*4 + 8 + 2+2+2+2
Coord = np.zeros((npnt,4))
# Write the smallest mesh spacing by default
for i in range(npnt):
    Coord[i,3] = resM;

# Spacing for ellipse segment support points
dsub = 0.025
# Center of the ellipsoid (0,0,0), unmodified
# Aligned coordinate system
Coord[1,0] = dsub; Coord[2,1] = dsub; Coord[3,2] = dsub;
# Rotated coordinate system, ellipseoid rotation
Coord[4,0] = dsub*np.cos(ePhi);    Coord[4,2] =-dsub*np.sin(ePhi)
Coord[5,0] = dsub*np.sin(ePhi);    Coord[5,2] = dsub*np.cos(ePhi)
# Rotated coordinate system, 45 degrees or pi/4 radians around y-axis
Coord[6,0] = dsub*np.cos(np.pi/4); Coord[6,2] = dsub*np.sin(np.pi/4)
Coord[7,0] =-dsub*np.sin(np.pi/4); Coord[7,2] = dsub*np.cos(np.pi/4)
# Rotated coordinate system, 45 degrees or pi/4 radians around z-axis
Coord[8,0] = dsub*np.cos(np.pi/4); Coord[8,1] = dsub*np.sin(np.pi/4)
Coord[9,0] =-dsub*np.sin(np.pi/4); Coord[9,1] = dsub*np.cos(np.pi/4)
# For planes z = +- y
nTES = npnt - 8
Coord[nTES+0,0] = dsub*np.cos(ePhi2); Coord[nTES+0,1] =-dsub*np.sin(ePhi2)*cs45 ; Coord[nTES+0,2] =-dsub*np.sin(ePhi2)*cs45    # hat{x} coordinates after 2nd rotation
Coord[nTES+1,0] = dsub*np.sin(ePhi2); Coord[nTES+1,1] = dsub*np.cos(ePhi2)*cs45;  Coord[nTES+1,2] = dsub*np.cos(ePhi2)*cs45    # hat{z} coordinates after 2nd rotation
Coord[nTES+2,0] = dsub*np.cos(ePhi2); Coord[nTES+2,1] = dsub*np.sin(ePhi2)*cs45 ; Coord[nTES+2,2] =-dsub*np.sin(ePhi2)*cs45    # hat{x} coordinates after 2nd rotation
Coord[nTES+3,0] = dsub*np.sin(ePhi2); Coord[nTES+3,1] =-dsub*np.cos(ePhi2)*cs45;  Coord[nTES+3,2] = dsub*np.cos(ePhi2)*cs45    # hat{z} coordinates after 2nd rotation
# For planes y = +- x
Coord[nTES+4,0] = dsub*np.cos(ePhi2)*cs45; Coord[nTES+4,1] = dsub*np.cos(ePhi2)*cs45 ; Coord[nTES+4,2] =-dsub*np.sin(ePhi2)    # hat{x} coordinates after 2nd rotation
Coord[nTES+5,0] = dsub*np.sin(ePhi2)*cs45; Coord[nTES+5,1] = dsub*np.sin(ePhi2)*cs45;  Coord[nTES+5,2] = dsub*np.cos(ePhi2)    # hat{z} coordinates after 2nd rotation
Coord[nTES+6,0] = dsub*np.cos(ePhi2)*cs45; Coord[nTES+6,1] =-dsub*np.cos(ePhi2)*cs45 ; Coord[nTES+6,2] =-dsub*np.sin(ePhi2)    # hat{x} coordinates after 2nd rotation
Coord[nTES+7,0] = dsub*np.sin(ePhi2)*cs45; Coord[nTES+7,1] =-dsub*np.sin(ePhi2)*cs45;  Coord[nTES+7,2] = dsub*np.cos(ePhi2)    # hat{z} coordinates after 2nd rotation

npntOffs = 10
# Generate support points for circular curved segments (modify only non-zero components)
Coord[npntOffs+0,0] = -0.5; Coord[npntOffs+1,0] = 0.5;
Coord[npntOffs+2,1] = -0.5; Coord[npntOffs+3,1] = 0.5;
Coord[npntOffs+4,2] = -0.5; Coord[npntOffs+5,2] = 0.5;

# Generate points for side circles, first around x axis, then y, then z, first - number offset
npntOffs = npntOffs+6
for i in range(2):
    for j in range(8):
        Coord[npntOffs+i*8+j,0] = Coord[i+npntOffs-6,0];
        Coord[npntOffs+i*8+j,1] = cR*np.cos(j*np.pi/4);
        Coord[npntOffs+i*8+j,2] = cR*np.sin(j*np.pi/4);
npntOffs = npntOffs+16
for i in range(2):
    for j in range(8):
        Coord[npntOffs+i*8+j,0] = cR*np.cos(j*np.pi/4);
        Coord[npntOffs+i*8+j,1] = Coord[i+npntOffs-20,1];
        Coord[npntOffs+i*8+j,2] = cR*np.sin(j*np.pi/4);
npntOffs = npntOffs+16
for i in range(2):
    for j in range(8):
        Coord[npntOffs+i*8+j,0] = cR*np.cos(j*np.pi/4);
        Coord[npntOffs+i*8+j,1] = cR*np.sin(j*np.pi/4);
        Coord[npntOffs+i*8+j,2] = Coord[i+npntOffs-34,2];
        
# Generate points for inner circles, using B-splines
npntOffs = npntOffs+16
# Constant auxiliary variables
sinP2 = np.sin(ePhi)**2
cosP2 = np.cos(ePhi)**2
scP   = np.sin(ePhi)*np.cos(ePhi)
# Ax = cos^2/a^2+sin^2/c^2; Bx = 2*z*sin*cos/c^2-2*z*sin*cos/a^2; Cx = z^2*sin^2/a^2+y^2/b^2+z^2*cos^2/c^2-1; y = cR*sin(t); z = cR*cos(t)
for i in range(4*Nbsp):
    tp            = 2.0*np.pi*i/(4.0*Nbsp)
    ycur          = cR*np.sin(tp)
    zcur          = cR*np.cos(tp)
    Ax            = cosP2/eA**2 + sinP2/eC**2
    Bx            = 2*zcur*scP*(1.0/eC**2 - 1.0/eA**2)
    Cx            = zcur**2*(sinP2/eA**2 + cosP2/eC**2) + ycur**2/eB**2-1
    Coord[npntOffs+i,0] = (-Bx + np.sqrt(Bx**2-4.0*Ax*Cx) )/(2.0*Ax)
    Coord[npntOffs+i,1] = ycur
    Coord[npntOffs+i,2] = zcur
for i in range(4*Nbsp):
    tp                   = 2.0*np.pi*i/(4.0*Nbsp)
    ycur                 = cR*np.sin(tp)
    zcur                 = cR*np.cos(tp)
    Ax                   = cosP2/eA**2 + sinP2/eC**2
    Bx                   = 2*zcur*scP*(1.0/eC**2 - 1.0/eA**2)
    Cx                   = zcur**2*(sinP2/eA**2 + cosP2/eC**2) + ycur**2/eB**2-1
    Coord[npntOffs+4*Nbsp+i,0] = (-Bx - np.sqrt(Bx**2-4.0*Ax*Cx) )/(2.0*Ax)
    Coord[npntOffs+4*Nbsp+i,1] = ycur
    Coord[npntOffs+4*Nbsp+i,2] = zcur
# Ay = 1/b^2; By = 0; Cy = (x^2*cos^2+z^2*sin^2-2*x*z*sin*cos)/a^2+(x^2*sin^2+z^2*cos^2+2*x*z*sin*cos)/c^2-1; x = cR*sin(t); z = cR*cos(t)
for i in range(4*Nbsp):
    tp                   = 2.0*np.pi*i/(4.0*Nbsp)
    xcur                 = cR*np.sin(tp)
    zcur                 = cR*np.cos(tp)
    Ay                   = 1.0/eB**2
    By                   = 0
    Cy                   = (xcur**2*cosP2 + zcur**2*sinP2 - 2*xcur*zcur*scP)/eA**2  \
                         + (xcur**2*sinP2 + zcur**2*cosP2 + 2*xcur*zcur*scP)/eC**2 - 1.0
    Coord[npntOffs+8*Nbsp+i,0] = xcur
    Coord[npntOffs+8*Nbsp+i,1] = np.sqrt(By**2-4.0*Ay*Cy)/(2.0*Ay)
    Coord[npntOffs+8*Nbsp+i,2] = zcur
for i in range(4*Nbsp):
    tp                    = 2.0*np.pi*i/(4.0*Nbsp)
    xcur                  = cR*np.sin(tp)
    zcur                  = cR*np.cos(tp)
    Ay                    = 1.0/eB**2
    By                    = 0
    Cy                    = (xcur**2*cosP2 + zcur**2*sinP2 - 2*xcur*zcur*scP)/eA**2  \
                          + (xcur**2*sinP2 + zcur**2*cosP2 + 2*xcur*zcur*scP)/eC**2 - 1.0
    Coord[npntOffs+12*Nbsp+i,0] = xcur
    Coord[npntOffs+12*Nbsp+i,1] =-np.sqrt(By**2-4.0*Ay*Cy)/(2.0*Ay)
    Coord[npntOffs+12*Nbsp+i,2] = zcur
# Az = sin^2/a^2+cos^2/c^2; Bz = 2*x*sin*cos/c^2-2*x*sin*cos/a^2; Cz = x^2*cos^2/a^2+y^2/b^2+x^2*sin^2/c^2-1; x = cR*sin(t); y = cR*cos(t)
for i in range(4*Nbsp):
    tp                    = 2.0*np.pi*i/(4.0*Nbsp)
    xcur                  = cR*np.sin(tp)
    ycur                  = cR*np.cos(tp)
    Az                    = sinP2/eA**2 + cosP2/eC**2
    Bz                    = 2*xcur*scP*(1.0/eC**2 - 1.0/eA**2)
    Cz                    = xcur**2*(cosP2/eA**2 + sinP2/eC**2) + ycur**2/eB**2 - 1.0
    Coord[npntOffs+16*Nbsp+i,0] = xcur
    Coord[npntOffs+16*Nbsp+i,1] = ycur
    Coord[npntOffs+16*Nbsp+i,2] = (-Bz + np.sqrt(Bz**2-4.0*Az*Cz) )/(2.0*Az)
for i in range(4*Nbsp):
    tp                    = 2.0*np.pi*i/(4.0*Nbsp)
    xcur                  = cR*np.sin(tp)
    ycur                  = cR*np.cos(tp)
    Az                    = sinP2/eA**2 + cosP2/eC**2
    Bz                    = 2*xcur*scP*(1.0/eC**2 - 1.0/eA**2)
    Cz                    = xcur**2*(cosP2/eA**2 + sinP2/eC**2) + ycur**2/eB**2 - 1.0
    Coord[npntOffs+20*Nbsp+i,0] = xcur
    Coord[npntOffs+20*Nbsp+i,1] = ycur
    Coord[npntOffs+20*Nbsp+i,2] = (-Bz - np.sqrt(Bz**2-4.0*Az*Cz) )/(2.0*Az)

# Mark half points of the connecting aligned arcs
npntOffs = npntOffs+6*4*Nbsp
# z = 0
Ax  = np.sqrt(1.0/(cosP2/eA**2+sinP2/eC**2))
tpm = np.arctan(eB/Ax)
Coord[npntOffs+0,0] = Ax*np.sin(tpm); Coord[npntOffs+0,1] = eB*np.cos(tpm); Coord[npntOffs+0,2] = 0
Coord[npntOffs+1,0] =-Ax*np.sin(tpm); Coord[npntOffs+1,1] = eB*np.cos(tpm); Coord[npntOffs+1,2] = 0
Coord[npntOffs+2,0] = Ax*np.sin(tpm); Coord[npntOffs+2,1] =-eB*np.cos(tpm); Coord[npntOffs+2,2] = 0
Coord[npntOffs+3,0] =-Ax*np.sin(tpm); Coord[npntOffs+3,1] =-eB*np.cos(tpm); Coord[npntOffs+3,2] = 0
XYi = Ax*np.sin(tpm)

# y = 0, x = +- z
npntOffs = npntOffs+4
xzcur = np.sqrt(1.0/( (np.cos(ePhi)-np.sin(ePhi))**2/eA**2 + (np.sin(ePhi)+np.cos(ePhi))**2/eC**2 ))
Coord[npntOffs+0,0] = xzcur; Coord[npntOffs+0,1] = 0; Coord[npntOffs+0,2] = xzcur
Coord[npntOffs+1,0] =-xzcur; Coord[npntOffs+1,1] = 0; Coord[npntOffs+1,2] =-xzcur
xzcur = np.sqrt(1.0/( (np.cos(ePhi)+np.sin(ePhi))**2/eA**2 + (np.sin(ePhi)-np.cos(ePhi))**2/eC**2 ))
Coord[npntOffs+2,0] = xzcur; Coord[npntOffs+2,1] = 0; Coord[npntOffs+2,2] =-xzcur
Coord[npntOffs+3,0] =-xzcur; Coord[npntOffs+3,1] = 0; Coord[npntOffs+3,2] = xzcur

# x = 0
npntOffs = npntOffs+4
Az  = np.sqrt(1.0/(sinP2/eA**2+cosP2/eC**2))
tpm = np.arctan(eB/Az);
Coord[npntOffs+0,0] = 0; Coord[npntOffs+0,1] = eB*np.cos(tpm); Coord[npntOffs+0,2] = Az*np.sin(tpm)
Coord[npntOffs+1,0] = 0; Coord[npntOffs+1,1] =-eB*np.cos(tpm); Coord[npntOffs+1,2] = Az*np.sin(tpm)
Coord[npntOffs+2,0] = 0; Coord[npntOffs+2,1] = eB*np.cos(tpm); Coord[npntOffs+2,2] =-Az*np.sin(tpm)
Coord[npntOffs+3,0] = 0; Coord[npntOffs+3,1] =-eB*np.cos(tpm); Coord[npntOffs+3,2] =-Az*np.sin(tpm)
YZi = Az*np.sin(tpm)

# x = y = z
npntOffs = npntOffs+4
curxyz = np.sqrt( 1.0/( (np.cos(ePhi)-np.sin(ePhi))**2/eA**2 + 1.0/eB**2 \
                      + (np.sin(ePhi)+np.cos(ePhi))**2/eC**2 ) )
Coord[npntOffs+0,0] = curxyz; Coord[npntOffs+0,1] = curxyz; Coord[npntOffs+0,2] = curxyz
Coord[npntOffs+2,0] = curxyz; Coord[npntOffs+2,1] =-curxyz; Coord[npntOffs+2,2] = curxyz
Coord[npntOffs+5,0] =-curxyz; Coord[npntOffs+5,1] = curxyz; Coord[npntOffs+5,2] =-curxyz
Coord[npntOffs+7,0] =-curxyz; Coord[npntOffs+7,1] =-curxyz; Coord[npntOffs+7,2] =-curxyz
curxyz = np.sqrt( 1.0/( (np.cos(ePhi)+np.sin(ePhi))**2/eA**2 + 1.0/eB**2 \
                      + (np.sin(ePhi)-np.cos(ePhi))**2/eC**2 ) )
Coord[npntOffs+1,0] =-curxyz; Coord[npntOffs+1,1] = curxyz; Coord[npntOffs+1,2] = curxyz
Coord[npntOffs+3,0] =-curxyz; Coord[npntOffs+3,1] =-curxyz; Coord[npntOffs+3,2] = curxyz
Coord[npntOffs+4,0] = curxyz; Coord[npntOffs+4,1] = curxyz; Coord[npntOffs+4,2] =-curxyz
Coord[npntOffs+6,0] = curxyz; Coord[npntOffs+6,1] =-curxyz; Coord[npntOffs+6,2] =-curxyz

# Additional points at the corners of the unit cell
npntOffs = npntOffs+8
for i in range(2):
    for j in range(2):
        for k in range(2):
            Coord[npntOffs+i*4+j*2+k,0] = -0.5 + i;
            Coord[npntOffs+i*4+j*2+k,1] = -0.5 + j;
            Coord[npntOffs+i*4+j*2+k,2] = -0.5 + k;
            Coord[npntOffs+i*4+j*2+k,3] = resM2;

# Shifting points such that z=0 corresponds to top of the structure
for i in range(npnt):
    Coord[i,2] = Coord[i,2] - 0.5

# Exporting all the points
fidF = open(GeoFile,'w')
for i in range(npnt):
    fidF.write('Point(%d) = {%.10f, %.10f, %.10f, %.10f};\n'%(i+1,Coord[i,0],Coord[i,1],Coord[i,2],Coord[i,3]))

# -------------------------------------------------------------
# Below here, the operations can be copied manually from GUI
# Constructions from loops also allowed, if feasable
# -------------------------------------------------------------
# Generate the connecting circle arcs around x, y and z axis
fidF.write('\
Circle(1) = {33, 13, 34};\n\
Circle(2) = {34, 13, 35};\n\
Circle(3) = {35, 13, 36};\n\
Circle(4) = {36, 13, 37};\n\
Circle(5) = {37, 13, 38};\n\
Circle(6) = {38, 13, 39};\n\
Circle(7) = {39, 13, 40};\n\
Circle(8) = {40, 13, 33};\n\
Circle(9) = {41, 14, 42};\n\
Circle(10) = {42, 14, 43};\n\
Circle(11) = {43, 14, 44};\n\
Circle(12) = {44, 14, 45};\n\
Circle(13) = {45, 14, 46};\n\
Circle(14) = {46, 14, 47};\n\
Circle(15) = {47, 14, 48};\n\
Circle(16) = {48, 14, 41};\n\
Circle(17) = {25, 12, 32};\n\
Circle(18) = {32, 12, 31};\n\
Circle(19) = {31, 12, 30};\n\
Circle(20) = {30, 12, 29};\n\
Circle(21) = {29, 12, 28};\n\
Circle(22) = {28, 12, 27};\n\
Circle(23) = {27, 12, 26};\n\
Circle(24) = {26, 12, 25};\n\
Circle(25) = {21, 11, 22};\n\
Circle(26) = {22, 11, 23};\n\
Circle(27) = {23, 11, 24};\n\
Circle(28) = {24, 11, 17};\n\
Circle(29) = {17, 11, 18};\n\
Circle(30) = {18, 11, 19};\n\
Circle(31) = {19, 11, 20};\n\
Circle(32) = {20, 11, 21};\n\
Circle(33) = {54, 15, 53};\n\
Circle(34) = {53, 15, 52};\n\
Circle(35) = {52, 15, 51};\n\
Circle(36) = {51, 15, 50};\n\
Circle(37) = {50, 15, 49};\n\
Circle(38) = {49, 15, 56};\n\
Circle(39) = {56, 15, 55};\n\
Circle(40) = {55, 15, 54};\n\
Circle(41) = {63, 16, 62};\n\
Circle(42) = {62, 16, 61};\n\
Circle(43) = {61, 16, 60};\n\
Circle(44) = {60, 16, 59};\n\
Circle(45) = {59, 16, 58};\n\
Circle(46) = {58, 16, 57};\n\
Circle(47) = {57, 16, 64};\n\
Circle(48) = {64, 16, 63};\n\
')

# Write Bsplines of the cylinder intersection with ellipsoid
fidF.write('\
BSpline(49) = {143, 142, 141, 140};\n\
BSpline(50) = {140, 139, 138, 137};\n\
BSpline(51) = {137, 160, 159, 158};\n\
BSpline(52) = {158, 157, 156, 155};\n\
BSpline(53) = {155, 154, 153, 152};\n\
BSpline(54) = {152, 151, 150, 149};\n\
BSpline(55) = {149, 148, 147, 146};\n\
BSpline(56) = {146, 145, 144, 143};\n\
BSpline(57) = {128, 129, 130, 131};\n\
BSpline(58) = {131, 132, 133, 134};\n\
BSpline(59) = {134, 135, 136, 113};\n\
BSpline(60) = {113, 114, 115, 116};\n\
BSpline(61) = {116, 117, 118, 119};\n\
BSpline(62) = {119, 120, 121, 122};\n\
BSpline(63) = {122, 123, 124, 125};\n\
BSpline(64) = {125, 126, 127, 128};\n\
BSpline(65) = {107, 108, 109, 110};\n\
BSpline(66) = {110, 110, 111, 112, 89};\n\
BSpline(67) = {89, 90, 91, 92};\n\
BSpline(68) = {92, 93, 94, 95};\n\
BSpline(69) = {95, 96, 97, 98};\n\
BSpline(70) = {98, 99, 100, 101};\n\
BSpline(71) = {101, 102, 103, 104};\n\
BSpline(72) = {104, 105, 106, 107};\n\
BSpline(73) = {71, 70, 69, 68};\n\
BSpline(74) = {68, 67, 66, 65};\n\
BSpline(75) = {65, 88, 87, 86};\n\
BSpline(76) = {86, 85, 84, 83};\n\
BSpline(77) = {83, 82, 81, 80};\n\
BSpline(78) = {80, 79, 78, 77};\n\
BSpline(79) = {77, 76, 75, 74};\n\
BSpline(80) = {74, 73, 72, 71};\n\
BSpline(81) = {197, 198, 199, 200};\n\
BSpline(82) = {200, 201, 202, 203};\n\
BSpline(83) = {203, 204, 205, 206};\n\
BSpline(84) = {206, 207, 208, 185};\n\
BSpline(85) = {185, 186, 187, 188};\n\
BSpline(86) = {188, 189, 190, 191};\n\
BSpline(87) = {191, 192, 193, 194};\n\
BSpline(88) = {194, 195, 196, 197};\n\
BSpline(89) = {161, 184, 183, 182};\n\
BSpline(90) = {182, 181, 180, 179};\n\
BSpline(91) = {179, 178, 177, 176};\n\
BSpline(92) = {176, 175, 174, 173};\n\
BSpline(93) = {173, 172, 171, 170};\n\
BSpline(94) = {170, 169, 168, 167};\n\
BSpline(95) = {167, 166, 165, 164};\n\
BSpline(96) = {164, 163, 162, 161};\n\
')

# Write connecting lines for the side cylinders
fidF.write('\
Line(97) = {21, 107};\n\
Line(98) = {20, 110};\n\
Line(99) = {19, 89};\n\
Line(100) = {18, 92};\n\
Line(101) = {17, 95};\n\
Line(102) = {24, 98};\n\
Line(103) = {23, 101};\n\
Line(104) = {22, 104};\n\
Line(105) = {27, 65};\n\
Line(106) = {26, 68};\n\
Line(107) = {25, 71};\n\
Line(108) = {32, 74};\n\
Line(109) = {31, 77};\n\
Line(110) = {30, 80};\n\
Line(111) = {29, 83};\n\
Line(112) = {28, 86};\n\
Line(113) = {44, 134};\n\
Line(114) = {45, 131};\n\
Line(115) = {46, 128};\n\
Line(116) = {47, 125};\n\
Line(117) = {48, 122};\n\
Line(118) = {41, 119};\n\
Line(119) = {116, 42};\n\
Line(120) = {43, 113};\n\
Line(121) = {35, 137};\n\
Line(122) = {34, 140};\n\
Line(123) = {33, 143};\n\
Line(124) = {40, 146};\n\
Line(125) = {39, 149};\n\
Line(126) = {38, 152};\n\
Line(127) = {37, 155};\n\
Line(128) = {36, 158};\n\
Line(129) = {61, 179};\n\
Line(130) = {60, 182};\n\
Line(131) = {59, 161};\n\
Line(132) = {58, 164};\n\
Line(133) = {57, 167};\n\
Line(134) = {64, 170};\n\
Line(135) = {63, 173};\n\
Line(136) = {62, 176};\n\
Line(137) = {53, 203};\n\
Line(138) = {54, 200};\n\
Line(139) = {55, 197};\n\
Line(140) = {56, 194};\n\
Line(141) = {49, 191};\n\
Line(142) = {50, 188};\n\
Line(143) = {51, 185};\n\
Line(144) = {52, 206};\n\
')

# Write the ellipse arcs around the ellipsoid
# Plane z = 0
fidF.write('\
Ellipse(145) = {143, 1, 2, 211};\n\
Ellipse(146) = {211, 1, 2, 83};\n\
Ellipse(147) = {71, 1, 2, 209};\n\
Ellipse(148) = {209, 1, 2, 119};\n\
Ellipse(149) = {131, 1, 2, 210};\n\
Ellipse(150) = {210, 1, 2, 95};\n\
Ellipse(151) = {107, 1, 2, 212};\n\
Ellipse(152) = {212, 1, 2, 155};\n\
')
# Plane y = 0
fidF.write('\
Ellipse(153) = {65, 1, 5, 213};\n\
Ellipse(154) = {213, 1, 5, 167};\n\
Ellipse(155) = {179, 1, 5, 216};\n\
Ellipse(156) = {216, 1, 5, 89};\n\
Ellipse(157) = {101, 1, 5, 214};\n\
Ellipse(158) = {214, 1, 5, 203};\n\
Ellipse(159) = {191, 1, 5, 215};\n\
Ellipse(160) = {215, 1, 5, 77};\n\
')
# Plane x = 0
fidF.write('\
Ellipse(161) = {137, 1, 3, 218};\n\
Ellipse(162) = {218, 1, 3, 173};\n\
Ellipse(163) = {161, 1, 3, 217};\n\
Ellipse(164) = {217, 1, 3, 113};\n\
Ellipse(165) = {125, 1, 3, 219};\n\
Ellipse(166) = {219, 1, 3, 185};\n\
Ellipse(167) = {197, 1, 3, 220};\n\
Ellipse(168) = {220, 1, 3, 149};\n\
')

# Plane z = x
fidF.write('\
Ellipse(169) = {140, 1, 7, 223};\n\
Ellipse(170) = {223, 1, 7, 213};\n\
Ellipse(171) = {213, 1, 7, 221};\n\
Ellipse(172) = {221, 1, 7, 116};\n\
Ellipse(173) = {128, 1, 7, 226};\n\
Ellipse(174) = {226, 1, 7, 214};\n\
Ellipse(175) = {214, 1, 7, 228};\n\
Ellipse(176) = {228, 1, 7, 152};\n\
')

# Plane z = -x
fidF.write('\
Ellipse(177) = {158, 1, 3, 224};\n\
Ellipse(178) = {224, 1, 3, 216};\n\
Ellipse(179) = {216, 1, 3, 222};\n\
Ellipse(180) = {222, 1, 3, 134};\n\
Ellipse(181) = {122, 1, 3, 225};\n\
Ellipse(182) = {225, 1, 3, 215};\n\
Ellipse(183) = {215, 1, 3, 227};\n\
Ellipse(184) = {227, 1, 3, 146};\n\
')

# Plane y = x
maP = nTES+5
fidF.write('\
Ellipse(185) = {164, 1, %d, 221};\n\
Ellipse(186) = {221, 1, %d, 209};\n\
Ellipse(187) = {209, 1, %d, 225};\n\
Ellipse(188) = {225, 1, %d, 188};\n\
Ellipse(189) = {200, 1, %d, 228};\n\
Ellipse(190) = {228, 1, %d, 212};\n\
Ellipse(191) = {212, 1, %d, 224};\n\
Ellipse(192) = {224, 1, %d, 176};\n\
'%(maP,maP,maP,maP,maP,maP,maP,maP))

# Plane y = -x
maP = nTES+7
fidF.write('\
Ellipse(193) = {182, 1, %d, 222};\n\
Ellipse(194) = {222, 1, %d, 210};\n\
Ellipse(195) = {210, 1, %d, 226};\n\
Ellipse(196) = {226, 1, %d, 206};\n\
Ellipse(197) = {194, 1, %d, 227};\n\
Ellipse(198) = {227, 1, %d, 211};\n\
Ellipse(199) = {211, 1, %d, 223};\n\
Ellipse(200) = {223, 1, %d, 170};\n\
'%(maP,maP,maP,maP,maP,maP,maP,maP))

# Plane z = y
maP = nTES+1
fidF.write('\
Ellipse(201) = {68, 1, %d, 221};\n\
Ellipse(202) = {221, 1, %d, 217};\n\
Ellipse(203) = {222, 1, %d, 92};\n\
Ellipse(204) = {217, 1, %d, 222};\n\
Ellipse(205) = {104, 1, %d, 228};\n\
Ellipse(206) = {228, 1, %d, 220};\n\
Ellipse(207) = {220, 1, %d, 227};\n\
Ellipse(208) = {227, 1, %d, 80};\n\
'%(maP,maP,maP,maP,maP,maP,maP,maP))

# Plane z = -y
maP = nTES+3
fidF.write('\
Ellipse(209) = {110, 1, %d, 224};\n\
Ellipse(210) = {224, 1, %d, 218};\n\
Ellipse(211) = {218, 1, %d, 223};\n\
Ellipse(212) = {223, 1, %d, 86};\n\
Ellipse(213) = {74, 1, %d, 225};\n\
Ellipse(214) = {225, 1, %d, 219};\n\
Ellipse(215) = {219, 1, %d, 226};\n\
Ellipse(216) = {226, 1, %d, 98};\n\
'%(maP,maP,maP,maP,maP,maP,maP,maP))


# Lines, denoting the outer boundary of the cell
fidF.write('\
Line(217) = {234, 236};\n\
Line(218) = {236, 232};\n\
Line(219) = {232, 230};\n\
Line(220) = {230, 234};\n\
Line(221) = {234, 233};\n\
Line(222) = {236, 235};\n\
Line(223) = {232, 231};\n\
Line(224) = {231, 229};\n\
Line(225) = {230, 229};\n\
Line(226) = {229, 233};\n\
Line(227) = {233, 235};\n\
Line(228) = {235, 231};\n\
')

# Add interior surfaces
fidF.write('\
Line Loop(230) = {21, 112, 76, -111};\n\
Line Loop(232) = {22, 105, 75, -112};\n\
Line Loop(234) = {23, 106, 74, -105};\n\
Line Loop(236) = {106, -73, -107, -24};\n\
Line Loop(238) = {80, -107, 17, 108};\n\
Line Loop(240) = {108, -79, -109, -18};\n\
Line Loop(242) = {109, -78, -110, -19};\n\
Line Loop(244) = {110, -77, -111, -20};\n\
Line Loop(246) = {66, -99, 31, 98};\n\
Line Loop(248) = {32, 97, 65, -98};\n\
Line Loop(250) = {97, -72, -104, -25};\n\
Line Loop(252) = {104, -71, -103, -26};\n\
Line Loop(254) = {103, -70, -102, -27};\n\
Line Loop(256) = {28, 101, 69, -102};\n\
Line Loop(258) = {29, 100, 68, -101};\n\
Line Loop(260) = {30, 99, 67, -100};\n\
Line Loop(262) = {4, 127, -52, -128};\n\
Line Loop(264) = {5, 126, -53, -127};\n\
Line Loop(266) = {126, 54, -125, -6};\n\
Line Loop(268) = {125, 55, -124, -7};\n\
Line Loop(270) = {124, 56, -123, -8};\n\
Line Loop(272) = {1, 122, -49, -123};\n\
Line Loop(274) = {122, 50, -121, -2};\n\
Line Loop(276) = {3, 128, -51, -121};\n\
Line Loop(278) = {58, -113, 12, 114};\n\
Line Loop(280) = {59, -120, 11, 113};\n\
Line Loop(282) = {120, 60, 119, 10};\n\
Line Loop(284) = {119, -9, 118, -61};\n\
Line Loop(286) = {118, 62, -117, 16};\n\
Line Loop(288) = {15, 117, 63, -116};\n\
Line Loop(290) = {14, 116, 64, -115};\n\
Line Loop(292) = {13, 115, 57, -114};\n\
Line Loop(294) = {144, 84, -143, -35};\n\
Line Loop(296) = {143, 85, -142, -36};\n\
Line Loop(298) = {86, -141, -37, 142};\n\
Line Loop(300) = {141, 87, -140, -38};\n\
Line Loop(302) = {140, 88, -139, -39};\n\
Line Loop(304) = {139, 81, -138, -40};\n\
Line Loop(306) = {33, 137, -82, -138};\n\
Line Loop(308) = {34, 144, -83, -137};\n\
Line Loop(310) = {90, -129, 43, 130};\n\
Line Loop(312) = {130, -89, -131, -44};\n\
Line Loop(314) = {131, -96, -132, -45};\n\
Line Loop(316) = {132, -95, -133, -46};\n\
Line Loop(318) = {94, -133, 47, 134};\n\
Line Loop(320) = {134, -93, -135, -48};\n\
Line Loop(322) = {41, 136, 92, -135};\n\
Line Loop(324) = {136, -91, -129, -42};\n\
Line Loop(326) = {155, 179, -193, 90};\n\
Line Loop(328) = {193, -204, -163, 89};\n\
Line Loop(330) = {179, 203, -67, -156};\n\
Line Loop(332) = {68, -150, -194, 203};\n\
Line Loop(334) = {194, -149, 58, -180};\n\
Line Loop(336) = {180, 59, -164, 204};\n\
Line Loop(338) = {164, 60, -172, 202};\n\
Line Loop(340) = {163, -202, -185, 96};\n\
Line Loop(342) = {185, -171, 154, 95};\n\
Line Loop(344) = {153, 171, -201, 74};\n\
Line Loop(346) = {201, 186, -147, 73};\n\
Line Loop(348) = {186, 148, -61, -172};\n\
Line Loop(350) = {155, -178, 192, -91};\n\
Line Loop(352) = {156, -66, 209, 178};\n\
Line Loop(354) = {209, -191, -151, 65};\n\
Line Loop(356) = {152, -52, 177, -191};\n\
Line Loop(358) = {177, 210, -161, 51};\n\
Line Loop(360) = {210, 162, -92, -192};\n\
Line Loop(362) = {162, 93, -200, -211};\n\
Line Loop(364) = {211, -169, 50, 161};\n\
Line Loop(366) = {169, -199, -145, 49};\n\
Line Loop(368) = {170, 154, -94, -200};\n\
Line Loop(370) = {170, -153, 75, -212};\n\
Line Loop(372) = {212, 76, -146, 199};\n\
Line Loop(374) = {158, 83, -196, 174};\n\
Line Loop(376) = {174, -157, -70, -216};\n\
Line Loop(378) = {195, 216, -69, -150};\n\
Line Loop(380) = {149, 195, -173, 57};\n\
Line Loop(382) = {173, -215, -165, 64};\n\
Line Loop(384) = {215, 196, 84, -166};\n\
Line Loop(386) = {165, -214, -181, 63};\n\
Line Loop(388) = {214, 166, 85, -188};\n\
Line Loop(390) = {181, -187, 148, 62};\n\
Line Loop(392) = {187, -213, 80, 147};\n\
Line Loop(394) = {213, 182, 160, 79};\n\
Line Loop(396) = {182, -159, -86, -188};\n\
Line Loop(398) = {160, -78, -208, -183};\n\
Line Loop(400) = {183, -197, -87, 159};\n\
Line Loop(402) = {208, -77, -146, -198};\n\
Line Loop(404) = {145, -198, 184, 56};\n\
Line Loop(406) = {184, -55, -168, 207};\n\
Line Loop(408) = {207, -197, 88, 167};\n\
Line Loop(410) = {168, -54, -176, 206};\n\
Line Loop(412) = {206, -167, 81, 189};\n\
Line Loop(414) = {176, -53, -152, -190};\n\
Line Loop(416) = {151, -190, -205, 72};\n\
Line Loop(418) = {205, -175, -157, 71};\n\
Line Loop(420) = {175, -189, 82, -158};\n\
')
# Interior surface indices - from 230 to 420 (step 2), inclusive, loop
for i in range(96):
    fidF.write('Ruled Surface(%d) = {%d};\n'%(230+i*2,230+i*2))

# Add closing surfaces
fidF.write('\
Line Loop(1001) = {225, -224, -223, 219, -32, -31, -30, -29, -28, -27, -26, -25};\n\
Line Loop(1002) = {221, 227, -222, -217, -17, -24, -23, -22, -21, -20, -19, -18};\n\
Line Loop(1003) = {226, -221, -220, 225, -8, -7, -6, -5, -4, -3, -2, -1};\n\
Line Loop(1004) = {218, 223, -228, -222, -16, -15, -14, -13, -12, -11, -10, -9};\n\
Line Loop(1005) = {227, 228, 224, 226, -36, -35, -34, -33, -40, -39, -38, -37};\n\
Line Loop(1006) = {218, 219, 220, 217, -44, -43, -42, -41, -48, -47, -46, -45};\n\
')
# Define plane surfaces
for i in range(6):
    fidF.write('Plane Surface(%d) = {%d};\n'%(1001+i,1001+i))

# Add closing surfaces, solid surfaces
fidF.write('\
Line Loop(11001) = {32, 25, 26, 27, 28, 29, 30, 31};\n\
Line Loop(11002) = {23, 24, 17, 18, 19, 20, 21, 22};\n\
Line Loop(11003) = {1, 2, 3, 4, 5, 6, 7, 8};\n\
Line Loop(11004) = {9, 10, 11, 12, 13, 14, 15, 16};\n\
Line Loop(11005) = {38, 39, 40, 33, 34, 35, 36, 37};\n\
Line Loop(11006) = {46, 47, 48, 41, 42, 43, 44, 45};\n\
')
# Define plane surfaces
for i in range(6):
    fidF.write('Plane Surface(%d) = {%d};\n'%(11001+i,11001+i))


# Do extension of unit cell down to the bottom of the coating
for i in range(Nbot-1):
    fidF.write('FreeSurface%d[] = Translate {0,0,-%d} {Duplicata{ Surface{'%(i+1,i+1) )
    for j in range(95):
        curInd = 230+j*2
        fidF.write('%d, '%(curInd))
    fidF.write('%d}; }};\n'%(curInd+2))
for i in range(Nbot-1):
    for j in range(4):
        fidF.write('P%dSurface%d[] = Translate {0,0,-%d} {Duplicata{ Surface{%d}; }};\n'%(1001+j,i+1,i+1,1001+j))
    for j in range(4):
        fidF.write('P%dSurface%d[] = Translate {0,0,-%d} {Duplicata{ Surface{%d}; }};\n'%(11001+j,i+1,i+1,11001+j))

# Translate the closing surface
fidF.write('BSurface[] = Translate {0,0,-%d} {Duplicata{ Surface{ 1005}; }};\n'%(Nbot-1))
fidF.write('BSurfaceS[]= Translate {0,0,-%d} {Duplicata{ Surface{11005}; }};\n'%(Nbot-1))
fidF.write('Delete { Surface{ 1005}; }\n')
fidF.write('Delete { Surface{11005}; }\n')

# Set finer mesh spacing at the interface
fidF.write('Characteristic Length {230, 232, 234, 236} = %.10f;\n'%resM)
    
# Define the volume, add the translated surface entities
fidF.write('Surface Loop(1301) = {')
for i in range(4):
    fidF.write('%d, '%(1001+i))
fidF.write('%d, '%(1006))
for i in range(96):
    curInd = 230+i*2
    fidF.write('%d, '%curInd)
for i in range(Nbot-1):
    for j in range(4):
        fidF.write('P%dSurface%d[], '%(1001+j,i+1))
for i in range(Nbot-1):
    fidF.write('FreeSurface%d[], '%(i+1))
fidF.write('BSurface[]};\n')
fidF.write('Volume(1401) = {1301};\n')

# Define the volume, add the translated surface entities, solid
fidF.write('Surface Loop(1302) = {')
for i in range(4):
    fidF.write('%d, '%(11001+i))
fidF.write('%d, '%(11006))
for i in range(96):
    curInd = 230+i*2
    fidF.write('%d, '%curInd)
for i in range(Nbot-1):
    for j in range(4):
        fidF.write('P%dSurface%d[], '%(11001+j,i+1))
for i in range(Nbot-1):
    fidF.write('FreeSurface%d[], '%(i+1))
fidF.write('BSurfaceS[]};\n')
fidF.write('Volume(1402) = {1302};\n')

# Define the periodic surfaces, direction vector to translate the right surface for matching
# In this case, left is top surface, right is bottom surface (must be raised by one unit)
fidF.write('Periodic Surface { 1002} = { 1001} Translate {1,0,0};\n');
fidF.write('Periodic Surface {11002} = {11001} Translate {1,0,0};\n');
fidF.write('Periodic Surface { 1004} = { 1003} Translate {0,1,0};\n');
fidF.write('Periodic Surface {11004} = {11003} Translate {0,1,0};\n');
# Translated periodic surfaces
for i in range(Nbot-1):
    fidF.write('Periodic Surface {P1002Surface%d[] } = {P1001Surface%d[] } Translate {1,0,0};\n'%(i+1,i+1));
    fidF.write('Periodic Surface {P11002Surface%d[]} = {P11001Surface%d[]} Translate {1,0,0};\n'%(i+1,i+1));
    fidF.write('Periodic Surface {P1004Surface%d[] } = {P1003Surface%d[] } Translate {0,1,0};\n'%(i+1,i+1));
    fidF.write('Periodic Surface {P11004Surface%d[]} = {P11003Surface%d[]} Translate {0,1,0};\n'%(i+1,i+1));

# Make sure that we output only entities denoted using physical labels
fidF.write('Mesh.SaveAll = 0;\n')
# Set the mesh resolution
fidF.write('Mesh.CharacteristicLengthMax = %.4f;\n'%(resM2))


# Append physical definitions
fidF = open(GeoFile,'a')
# Define the physical surfaces
# Periodic boundaries, including the translated ones (need to have separate
# pysical labels, because they are disconnected)
fidF.write('Physical Surface(3001) = {1001};\n')
for i in range(Nbot-1):
    fidF.write('Physical Surface(3001%d) = {P1001Surface%d[]};\n'%(i+1,i+1))
fidF.write('Physical Surface(3002) = {1002};\n')
for i in range(Nbot-1):
    fidF.write('Physical Surface(3002%d) = {P1002Surface%d[]};\n'%(i+1,i+1))
fidF.write('Physical Surface(3003) = {1003};\n')
for i in range(Nbot-1):
    fidF.write('Physical Surface(3003%d) = {P1003Surface%d[]};\n'%(i+1,i+1))
fidF.write('Physical Surface(3004) = {1004};\n')
for i in range(Nbot-1):
    fidF.write('Physical Surface(3004%d) = {P1004Surface%d[]};\n'%(i+1,i+1))
fidF.write('Physical Surface(3005) = {BSurface[]};\n' )
fidF.write('Physical Surface(3006) = {1006};\n' )
# Free boundaries, including the translated ones
fidF.write('Physical Surface(3007) = {')
for i in range(96):
    curInd = 230+i*2
    fidF.write('%d, '%(curInd))
for i in range(Nbot-2):
    fidF.write('FreeSurface%d[], '%(i+1))
fidF.write('FreeSurface%d[]};\n'%(Nbot-1))
# Define the physical volume
fidF.write('Physical Volume(4001) = {1401};\n\n' )
fidF.close()

# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------------------------
# Construct interface cell, top part
# -----------------------------------------------------------------------------------------------------
GeoFile = '3D_MSE_gen_geomT.geo'
# Determine the z coordinates
zmin = 0.0
zmax = ztop

# Compute the count of points {cell-corners} + {circle points surf}
npnt  = 12 + 3

Coord = np.zeros((npnt,4))
# Write the smallest mesh spacing by default
for i in range(npnt):
    Coord[i,3] = resM;

# Set coordinates for corner points
Coord[0,0] = -0.5; Coord[0,1] = -0.5; Coord[0,2] = zmin;
Coord[1,0] = -0.5; Coord[1,1] =  0.5; Coord[1,2] = zmin;
Coord[2,0] =  0.5; Coord[2,1] = -0.5; Coord[2,2] = zmin;
Coord[3,0] =  0.5; Coord[3,1] =  0.5; Coord[3,2] = zmin;
Coord[4,0] = -0.5; Coord[4,1] = -0.5; Coord[4,2] = zmin+1.0;   Coord[4,3]  = resM2;
Coord[5,0] = -0.5; Coord[5,1] =  0.5; Coord[5,2] = zmin+1.0;   Coord[5,3]  = resM2;
Coord[6,0] =  0.5; Coord[6,1] = -0.5; Coord[6,2] = zmin+1.0;   Coord[6,3]  = resM2;
Coord[7,0] =  0.5; Coord[7,1] =  0.5; Coord[7,2] = zmin+1.0;   Coord[7,3]  = resM2;
Coord[8, 0] = -0.5; Coord[8, 1] = -0.5; Coord[8, 2] = ztop;    Coord[8, 3] = resM2;
Coord[9, 0] = -0.5; Coord[9, 1] =  0.5; Coord[9, 2] = ztop;    Coord[9, 3] = resM2;
Coord[10,0] =  0.5; Coord[10,1] = -0.5; Coord[10,2] = ztop;    Coord[10,3] = resM2;
Coord[11,0] =  0.5; Coord[11,1] =  0.5; Coord[11,2] = ztop;    Coord[11,3] = resM2;
# Set coordinates for circle points
Coord[12,0] =  0.0; Coord[12,1] =  0.0; Coord[12,2] = zmin;
Coord[13,0] =  cR;  Coord[13,1] =  0.0; Coord[13,2] = zmin;
Coord[14,0] = -cR;  Coord[14,1] =  0.0; Coord[14,2] = zmin;

# Exporting all the points
fidF = open(GeoFile,'w')
for i in range(npnt):
    fidF.write('Point(%d) = {%.10f, %.10f, %.10f, %.10f};\n'%(i+1,Coord[i,0],Coord[i,1],Coord[i,2],Coord[i,3]))

# Define lines
fidF.write('\
Circle(1) = {14, 13, 15};\n\
Circle(2) = {15, 13, 14};\n\
')
fidF.write('\
Line(7) = {1, 3};\n\
Line(8) = {3, 4};\n\
Line(9) = {4, 2};\n\
Line(10) = {2, 1};\n\
Line(11) = {11, 12};\n\
Line(12) = {12, 10};\n\
Line(13) = {10, 9};\n\
Line(14) = {9, 11};\n\
Line(15) = {3, 7};\n\
Line(16) = {7, 11};\n\
Line(17) = {4, 8};\n\
Line(18) = {8, 12};\n\
Line(19) = {2, 6};\n\
Line(20) = {6, 10};\n\
Line(21) = {1, 5};\n\
Line(22) = {5, 9};\n\
')
# Define surfaces
fidF.write('\
Line Loop(24) = {2, 1};\n\
Line Loop(31) = {10, 7, 8, 9, -2, -1};\n\
Plane Surface(24) = {24};\n\
')
fidF.write('\
Line Loop(33) = {8, 17, 18, -11, -16, -15};\n\
Line Loop(35) = {7, 15, 16, -14, -22, -21};\n\
Line Loop(37) = {10, 21, 22, -13, -20, -19};\n\
Line Loop(39) = {12, -20, -19, -9, 17, 18};\n\
Line Loop(41) = {11, 12, 13, 14};\n\
Plane Surface(31) = {31};\n\
Plane Surface(33) = {33};\n\
Plane Surface(35) = {35};\n\
Plane Surface(37) = {37};\n\
Plane Surface(39) = {39};\n\
Plane Surface(41) = {41};\n\
')

# Define volume
fidF.write('Surface Loop(43) = {41, 33, 31, 37, 35, 39, 24};\n')
fidF.write('Volume(43) = {43};\n')

# Define the physical surfaces
# Periodic boundaries, including the translated ones (need to have separate
# pysical labels, because they are disconnected)
fidF.write('Physical Surface(3001) = {37};\n')        # (y,z)-plane, x-
fidF.write('Physical Surface(3002) = {33};\n')        # (y,z)-plane, x+
fidF.write('Physical Surface(3003) = {35};\n')        # (x,z)-plane, y-
fidF.write('Physical Surface(3004) = {39};\n')        # (x,z)-plane, y+
fidF.write('Physical Surface(3005) = {31};\n')        # Bottom
fidF.write('Physical Surface(3006) = {41};\n')        # Top
fidF.write('Physical Surface(3007) = {24};\n' )       # No-slip walls

# Define the physical volume
fidF.write('Physical Volume(4001) = {43};\n\n' )

# Define the periodic surfaces, direction vector to translate the right surface for matching
# In this case, left is top surface, right is bottom surface (must be raised by one unit)
fidF.write('Periodic Surface {33} = {37} Translate {1,0,0};\n');
fidF.write('Periodic Surface {39} = {35} Translate {0,1,0};\n');

# Set the mesh resolution
fidF.write('Mesh.CharacteristicLengthMax = %.4f;\n'%(resM2))

fidF.close()


# -----------------------------------------------------------------------------------------------------
# Construct Lagrange multiplier space
# -----------------------------------------------------------------------------------------------------
GeoFile = '3D_MSE_gen_geomL.geo'
# Set half the resolution
resM  = 0.5*resM
resM2 = 0.5*resM2
# Determine the z coordinates
zmin = 0.0
zmax = resM

# Compute the count of points {cell-corners} + {circle points bot}
npnt  = 8 + 3

Coord = np.zeros((npnt,4))
# Write the smallest mesh spacing by default
for i in range(npnt):
    Coord[i,3] = resM;
# Set coordinates for corner points
Coord[0,0] = -0.5; Coord[0,1] = -0.5; Coord[0,2] = zmin;
Coord[1,0] = -0.5; Coord[1,1] =  0.5; Coord[1,2] = zmin;
Coord[2,0] =  0.5; Coord[2,1] = -0.5; Coord[2,2] = zmin;
Coord[3,0] =  0.5; Coord[3,1] =  0.5; Coord[3,2] = zmin;
Coord[4,0] = -0.5; Coord[4,1] = -0.5; Coord[4,2] = zmin+resM;
Coord[5,0] = -0.5; Coord[5,1] =  0.5; Coord[5,2] = zmin+resM;
Coord[6,0] =  0.5; Coord[6,1] = -0.5; Coord[6,2] = zmin+resM;
Coord[7,0] =  0.5; Coord[7,1] =  0.5; Coord[7,2] = zmin+resM;
# Set coordinates for circles points
Coord[8, 0] =  0.0; Coord[8, 1] =  0.0; Coord[8, 2] = zmin;
Coord[9, 0] =  cR;  Coord[9, 1] =  0.0; Coord[9, 2] = zmin;
Coord[10,0] = -cR;  Coord[10,1] =  0.0; Coord[10,2] = zmin;

# Exporting all the points
fidF = open(GeoFile,'w')
for i in range(npnt):
    fidF.write('Point(%d) = {%.10f, %.10f, %.10f, %.10f};\n'%(i+1,Coord[i,0],Coord[i,1],Coord[i,2],Coord[i,3]))

# Define lines
fidF.write('\
Line(1) = {1,2};\n\
Line(2) = {1,3};\n\
Line(3) = {3,4};\n\
Line(4) = {2,4};\n\
Line(5) = {5,6};\n\
Line(6) = {5,7};\n\
Line(7) = {7,8};\n\
Line(8) = {6,8};\n\
Line(9) = {1,5};\n\
Line(10) = {2,6};\n\
Line(11) = {3,7};\n\
Line(12) = {4,8};\n\
')
fidF.write('Circle(13) = {10, 9, 11};\n')
fidF.write('Circle(14) = {11, 9, 10};\n')
# Define surfaces
fidF.write('\
Line Loop(19) = {13, 14};\n\
Line Loop(20) = {13, 14, 1,-2,-3,4};\n\
Plane Surface(19) = {19};\n\
')
fidF.write('\
Line Loop(21) = {5,-6,-7, 8};\n\
Line Loop(22) = {1, 10, -5, -9};\n\
Line Loop(23) = {2, 11, -6, -9};\n\
Line Loop(24) = {7, -12, -3, 11};\n\
Line Loop(25) = {12, -8, -10, 4};\n\
Plane Surface(20) = {20};\n\
Plane Surface(21) = {21};\n\
Plane Surface(22) = {22};\n\
Plane Surface(23) = {23};\n\
Plane Surface(24) = {24};\n\
Plane Surface(25) = {25};\n\
')

# Define volume
fidF.write('Surface Loop(43) = {20,21,22,23,24,25,19};\n')
fidF.write('Volume(43) = {43};\n')

# Define the physical surfaces
# Periodic boundaries, including the translated ones (need to have separate
# pysical labels, because they are disconnected)
fidF.write('Physical Surface(3001) = {22};\n')        # (y,z)-plane, x-
fidF.write('Physical Surface(3002) = {24};\n')        # (y,z)-plane, x+
fidF.write('Physical Surface(3003) = {23};\n')        # (x,z)-plane, y-
fidF.write('Physical Surface(3004) = {25};\n')        # (x,z)-plane, y+
fidF.write('Physical Surface(3005) = {20};\n')        # Bottom
fidF.write('Physical Surface(3006) = {21};\n')        # Top
fidF.write('Physical Surface(3007) = {19};\n')        # Interior walls, not needed for Lagrange multipliers
# Define the physical volume
fidF.write('Physical Volume(4001) = {43};\n\n' )
# Define the periodic surfaces, direction vector to translate the right surface for matching
# In this case, left is top surface, right is bottom surface (must be raised by one unit)
fidF.write('Periodic Surface {24} = {22} Translate {1,0,0};\n');
fidF.write('Periodic Surface {25} = {23} Translate {0,1,0};\n');

# Set the mesh resolution
fidF.write('Mesh.CharacteristicLengthMax = %.4f;\n'%(resM2))

fidF.close()
# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------
