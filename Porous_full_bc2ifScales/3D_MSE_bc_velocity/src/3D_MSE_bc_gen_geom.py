#!/usr/bin/python2.7
# ----------------------------------------------- COPYRIGHT --------------------------------------
# Copyright 2016
# Ugis Lacis, ugis.lacis@gmail.com
# Shervin Bagheri, shervin.bagheri@mech.kth.se
# -------------------------------------------- LICENSE LGPLv3 ------------------------------------
# This file is part of Porous_full_bc2ifScales.
#
# Porous_full_bc2ifScales is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Porous_full_bc2ifScales is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Porous_full_bc2ifScales. If not, see <http://www.gnu.org/licenses/>.
# ---------------------------------------------- DESCRIPTION -------------------------------------
#
# This software is generating the geometry definition for Gmsh mesher. It defines all necessary
# physical domains for the FreeFem++ solver. It is meant to be executed using Python. It has been
# tested to produce a valid Gmsh geometry file using Python version 2.7, numpy version 1.6.2 and
# scipy version 0.11.0.
#
# Execute this file by typing "./3D_MSE_bc_gen_geom.py" in command line (Unix). Make sure that
# the top line in this file includes valid path to Python and that case file
# <3D_MSE_bc_parameters.in> exists. For more information, consult documentation in <doc/>

import sys
import numpy as np
import scipy.optimize as opt

# Read the case parameters
try:
    InputParam = np.loadtxt("3D_MSE_bc_parameters.in")
except IOError:
    sys.exit("Error: expecting parameter file <3D_MSE_bc_parameters.in>")
ztop   =     InputParam[0]
Nbot   = int(InputParam[1])
thetas =     InputParam[2]
dt_d   =     InputParam[3]
resM   =     InputParam[4]
resM2  =     InputParam[5]
# Display used parameters
print "Finished reading parameters, ztop = "+str(ztop)+", Nbot = "+str(Nbot)+", thetas = "+str(thetas)+", dt_d = "+str(dt_d)



# Define a target function for getting appropriate volume fraction
def func(hR):
    # Coefficients
    cf1 = dt_d**2
    cf2 = np.sqrt(1 - cf1)
    cf3 = 1 - cf2
    # Difference function
    return 4.0/3.0*np.pi*hR**3 + 6.0*np.pi*cf1*hR**2*(0.5-hR*cf2)-4.0/6.0*np.pi*hR**3*cf3*(3.0*cf1 + cf3**2) - thetas

# Finding radii of hub-sphere and connecting cylinders
hR = opt.bisect(func, 0.0, 0.5)
cR = dt_d*hR
# Display the result to screen
print "Obtained hR = "+str(hR)+", cR = "+str(cR)+", for porosity = "+str(1.0-thetas)

# Auxiliary variable, intersection distance from the center of the unit cell
dR   = np.sqrt(hR**2 - cR**2)
zmin = -dR-0.5 - (Nbot-1.0)

# Save the obtained minimum coordinate to a new parameter file
try:
    InputParam[6] = zmin
except IndexError:
    InputParam = np.append(InputParam, zmin)
np.savetxt("3D_MSE_bc_parameters.txt",InputParam)
# Append the comment line at the end
fid = open("3D_MSE_bc_parameters.txt", "a")
fid.write("# ztop, Nbot, theta, cR/hR, dsMin, dsMax, zmin\n")
fid.close()



# -----------------------------------------------------------------------------------------------------
# Construct interior unit cell, below the interface cell
# -----------------------------------------------------------------------------------------------------
GeoFile = '3D_MSE_bc_gen_geomI.geo'
# Determine the z coordinates
zmin = -dR-0.5 - (Nbot    )
zmax = -dR-0.5 - (Nbot-1.0)

# Compute the amount of points {support} + {outer circles} + {inner circles} + {mid-sphere-points} + {cell-corners}
npnt  = 13 + 6*8 + 6*8 + 2*8+4 + 8
Coord = np.zeros((npnt,4))
# Write the smallest mesh spacing by default
for i in range(npnt):
    Coord[i,3] = resM;
# Generate support points for curved segments (modify only non-zero components)
Coord[1,0] = -0.5; Coord[2,0] = -dR; Coord[3,0] = dR; Coord[4,0] = 0.5;
Coord[5,1] = -0.5; Coord[6,1] = -dR; Coord[7,1] = dR; Coord[8,1] = 0.5;
Coord[9,2] = -0.5; Coord[10,2]= -dR; Coord[11,2]= dR; Coord[12,2]= 0.5;

# Generate points for side circles, first around x axis, then y, then z, first - number offset
npntOffs = 13
for i in range(4):
    for j in range(8):
        Coord[npntOffs+i*8+j,0] = Coord[i+1,0];
        Coord[npntOffs+i*8+j,1] = cR*np.cos(j*np.pi/4);
        Coord[npntOffs+i*8+j,2] = cR*np.sin(j*np.pi/4);
npntOffs = 13+32
for i in range(4):
    for j in range(8):
        Coord[npntOffs+i*8+j,0] = cR*np.cos(j*np.pi/4);
        Coord[npntOffs+i*8+j,1] = Coord[i+5,1];
        Coord[npntOffs+i*8+j,2] = cR*np.sin(j*np.pi/4);
npntOffs = 13+64
for i in range(4):
    for j in range(8):
        Coord[npntOffs+i*8+j,0] = cR*np.cos(j*np.pi/4);
        Coord[npntOffs+i*8+j,1] = cR*np.sin(j*np.pi/4);
        Coord[npntOffs+i*8+j,2] = Coord[i+9,2];

npntOffs = 13+96
for j in range(8):
    Coord[npntOffs+j,0] = hR*np.cos(np.pi/4)*np.cos(j*np.pi/4);
    Coord[npntOffs+j,1] = hR*np.cos(np.pi/4)*np.sin(j*np.pi/4);
    Coord[npntOffs+j,2] = hR*np.sin(np.pi/4);
npntOffs = 13+104
for j in range(4):
    Coord[npntOffs+j,0] = hR*np.cos(np.pi/4+j*np.pi/2);
    Coord[npntOffs+j,1] = hR*np.sin(np.pi/4+j*np.pi/2);
    Coord[npntOffs+j,2] = 0.0;
npntOffs = 13+108
for j in range(8):
    Coord[npntOffs+j,0] = hR*np.cos(np.pi/4)*np.cos(j*np.pi/4);
    Coord[npntOffs+j,1] = hR*np.cos(np.pi/4)*np.sin(j*np.pi/4);
    Coord[npntOffs+j,2] =-hR*np.sin(np.pi/4);

# Additional points at the corners of the unit cell
npntOffs = 13+116
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
fid = open(GeoFile,'w')
for i in range(npnt):
    fid.write('Point(%d) = {%.10f, %.10f, %.10f, %.10f};\n'%(i+1,Coord[i,0],Coord[i,1],Coord[i,2],Coord[i,3]))

# Below here, some operations are copied manually from GUI, but loops are used if found useful
# Generate the connecting circle arcs around x, y and z axis
offsArr = np.array([0, 8, 16, 24, 24, 28, 36, 44, 48, 52, 60, 68])
for i in range(12):
    for j in range(7):
        curLid = i*8 +j+1
        curPid = 13+curLid
        fid.write('Circle(%d) = {%d, %d, %d};\n'%( curLid, curPid , i+2 , curPid+1))
    fid.write('Circle(%d) = {%d, %d, %d};\n'%( curLid+1, curPid+1 , i+2 , curPid-6))

# Write connecting arcs of the sphere (copied from GUI)
fid.write('\
Circle(97) = {127, 1, 91};\n\
Circle(98) = {127, 1, 128};\n\
Circle(99) = {127, 1, 59};\n\
Circle(100) = {127, 1, 120};\n\
Circle(101) = {127, 1, 27};\n\
Circle(102) = {127, 1, 126};\n\
Circle(103) = {126, 1, 90};\n\
Circle(104) = {126, 1, 125};\n\
Circle(106) = {126, 1, 28};\n\
Circle(107) = {125, 1, 89};\n\
Circle(108) = {125, 1, 124};\n\
Circle(109) = {125, 1, 67};\n\
Circle(110) = {125, 1, 119};\n\
Circle(111) = {124, 1, 88};\n\
Circle(112) = {124, 1, 123};\n\
Circle(113) = {124, 1, 68};\n\
Circle(114) = {125, 1, 29};\n\
Circle(115) = {123, 1, 87};\n\
Circle(116) = {123, 1, 122};\n\
Circle(117) = {123, 1, 37};\n\
Circle(118) = {123, 1, 118};\n\
Circle(119) = {123, 1, 69};\n\
Circle(120) = {122, 1, 86};\n\
Circle(121) = {122, 1, 129};\n\
Circle(122) = {122, 1, 36};\n\
Circle(123) = {129, 1, 93};\n\
Circle(124) = {129, 1, 128};\n\
Circle(125) = {129, 1, 61};\n\
Circle(126) = {129, 1, 121};\n\
Circle(127) = {129, 1, 35};\n\
Circle(128) = {128, 1, 92};\n\
Circle(129) = {128, 1, 60};\n\
Circle(130) = {116, 1, 100};\n\
Circle(131) = {116, 1, 117};\n\
Circle(132) = {116, 1, 115};\n\
Circle(133) = {116, 1, 56};\n\
Circle(134) = {117, 1, 101};\n\
Circle(135) = {117, 1, 110};\n\
Circle(136) = {117, 1, 33};\n\
Circle(137) = {117, 1, 121};\n\
Circle(138) = {117, 1, 55};\n\
Circle(139) = {121, 1, 34};\n\
Circle(140) = {121, 1, 54};\n\
Circle(141) = {110, 1, 94};\n\
Circle(142) = {110, 1, 111};\n\
Circle(143) = {110, 1, 32};\n\
Circle(144) = {111, 1, 95};\n\
Circle(145) = {111, 1, 112};\n\
Circle(146) = {111, 1, 63};\n\
Circle(147) = {111, 1, 118};\n\
Circle(148) = {111, 1, 31};\n\
Circle(149) = {118, 1, 62};\n\
Circle(150) = {118, 1, 30};\n\
Circle(151) = {112, 1, 96};\n\
Circle(152) = {112, 1, 64};\n\
Circle(153) = {112, 1, 113};\n\
Circle(154) = {113, 1, 97};\n\
Circle(155) = {113, 1, 114};\n\
Circle(156) = {113, 1, 23};\n\
Circle(157) = {113, 1, 119};\n\
Circle(158) = {113, 1, 65};\n\
Circle(159) = {114, 1, 98};\n\
Circle(160) = {114, 1, 24};\n\
Circle(161) = {114, 1, 115};\n\
Circle(162) = {115, 1, 99};\n\
Circle(163) = {115, 1, 57};\n\
Circle(164) = {115, 1, 25};\n\
Circle(165) = {115, 1, 120};\n\
Circle(166) = {120, 1, 58};\n\
Circle(167) = {120, 1, 26};\n\
Circle(168) = {119, 1, 22};\n\
Circle(169) = {119, 1, 66};\n\
')

# Add the lines for cylinders
for i in range(6):
    nbndOffs = 14+i*16
    npntOffs = 170+i*8
    for j in range(8):
        fid.write('Line(%d) = {%d, %d};\n'%( npntOffs+j, nbndOffs+j , nbndOffs+j+8))

# Add closing surfaces for the fluid problem, define lines, close the surface
# Copied from the GUI output
# Line definitions
fid.write('\
Line(218) = {133, 132};\n\
Line(219) = {132, 130};\n\
Line(220) = {130, 131};\n\
Line(221) = {131, 133};\n\
Line(222) = {133, 137};\n\
Line(223) = {137, 135};\n\
Line(224) = {135, 131};\n\
Line(225) = {137, 136};\n\
Line(226) = {136, 132};\n\
Line(227) = {136, 134};\n\
Line(228) = {134, 130};\n\
Line(229) = {134, 135};\n\
')
# Line loop definitions
fid.write('\
Line Loop(1001) = {218, 219, 220, 221, -2, -1, -8, -7, -6, -5, -4, -3};\n\
Line Loop(1002) = {225, 227, 229, -223, -29, -28, -27, -26, -25, -32, -31, -30};\n\
Line Loop(1003) = {229, 224, -220, -228, -40, -39, -38, -37, -36, -35, -34, -33};\n\
Line Loop(1004) = {225, 226, -218, 222, -57, -64, -63, -62, -61, -60, -59, -58};\n\
Line Loop(1005) = {226, 219, -228, -227, -66, -65, -72, -71, -70, -69, -68, -67};\n\
Line Loop(1006) = {223, 224, 221, 222, -96, -95, -94, -93, -92, -91, -90, -89};\n\
')
# Define plane surfaces
for i in range(6):
    fid.write('Plane Surface(%d) = {%d};\n'%(1001+i,1001+i))

# Rest of surfaces are taken from GUI
# First we do line loops, then surfaces
fid.write('\
Line Loop(1008) = {6, 176, -14, -175};\n\
Line Loop(1010) = {5, 175, -13, -174};\n\
Line Loop(1012) = {174, -12, -173, 4};\n\
Line Loop(1014) = {3, 173, -11, -172};\n\
Line Loop(1016) = {2, 172, -10, -171};\n\
Line Loop(1018) = {1, 171, -9, -170};\n\
Line Loop(1020) = {177, 16, -170, -8};\n\
Line Loop(1022) = {7, 177, -15, -176};\n\
Line Loop(1024) = {183, 30, -184, -22};\n\
Line Loop(1026) = {184, 31, -185, -23};\n\
Line Loop(1028) = {185, 32, -178, -24};\n\
Line Loop(1030) = {178, 25, -179, -17};\n\
Line Loop(1032) = {179, 26, -180, -18};\n\
Line Loop(1034) = {180, 27, -181, -19};\n\
Line Loop(1036) = {181, 28, -182, -20};\n\
Line Loop(1038) = {183, -29, -182, 21};\n\
Line Loop(1040) = {192, 47, -193, -39};\n\
Line Loop(1042) = {193, 48, -186, -40};\n\
Line Loop(1044) = {33, 187, -41, -186};\n\
Line Loop(1046) = {34, 188, -42, -187};\n\
Line Loop(1048) = {35, 189, -43, -188};\n\
Line Loop(1050) = {189, 44, -190, -36};\n\
Line Loop(1052) = {190, 45, -191, -37};\n\
Line Loop(1054) = {191, 46, -192, -38};\n\
Line Loop(1056) = {60, -198, -52, 197};\n\
Line Loop(1058) = {198, 61, -199, -53};\n\
Line Loop(1060) = {62, -200, -54, 199};\n\
Line Loop(1062) = {55, 201, -63, -200};\n\
Line Loop(1064) = {56, 194, -64, -201};\n\
Line Loop(1066) = {194, 57, -195, -49};\n\
Line Loop(1068) = {195, 58, -196, -50};\n\
Line Loop(1070) = {196, 59, -197, -51};\n\
Line Loop(1072) = {67, 205, -75, -204};\n\
Line Loop(1074) = {204, -74, -203, 66};\n\
Line Loop(1076) = {76, -206, -68, 205};\n\
Line Loop(1078) = {77, -207, -69, 206};\n\
Line Loop(1080) = {207, 78, -208, -70};\n\
Line Loop(1082) = {208, 79, -209, -71};\n\
Line Loop(1084) = {209, 80, -202, -72};\n\
Line Loop(1086) = {202, 73, -203, -65};\n\
Line Loop(1088) = {212, 91, -213, -83};\n\
Line Loop(1090) = {213, 92, -214, -84};\n\
Line Loop(1092) = {93, -215, -85, 214};\n\
Line Loop(1094) = {215, 94, -216, -86};\n\
Line Loop(1096) = {216, 95, -217, -87};\n\
Line Loop(1098) = {217, 96, -210, -88};\n\
Line Loop(1100) = {210, 89, -211, -81};\n\
Line Loop(1102) = {212, -90, -211, 82};\n\n\
Line Loop(1104) = {84, -159, -155, 154};\n\
Line Loop(1106) = {159, 85, -162, -161};\n\
Line Loop(1108) = {162, 86, -130, 132};\n\
Line Loop(1110) = {154, -83, -151, 153};\n\
Line Loop(1112) = {130, 87, -134, -131};\n\
Line Loop(1114) = {134, 88, -141, -135};\n\
Line Loop(1116) = {142, 144, -81, -141};\n\
Line Loop(1118) = {144, 82, -151, -145};\n\
Line Loop(1120) = {143, 19, -136, 135};\n\
Line Loop(1122) = {136, 20, -139, -137};\n\
Line Loop(1124) = {137, 140, 41, -138};\n\
Line Loop(1126) = {138, 42, -133, 131};\n\
Line Loop(1128) = {133, 43, -163, -132};\n\
Line Loop(1130) = {163, 44, -166, -165};\n\
Line Loop(1132) = {165, 167, -12, -164};\n\
Line Loop(1134) = {164, -11, -160, 161};\n\
Line Loop(1136) = {160, -10, -156, 155};\n\
Line Loop(1138) = {156, -9, -168, -157};\n\
Line Loop(1140) = {157, 169, -52, -158};\n\
Line Loop(1142) = {158, -51, -152, 153};\n\
Line Loop(1144) = {152, -50, -146, 145};\n\
Line Loop(1146) = {146, -49, -149, -147};\n\
Line Loop(1148) = {147, 150, 17, -148};\n\
Line Loop(1150) = {148, 18, -143, 142};\n\
Line Loop(1152) = {149, -56, -119, 118};\n\
Line Loop(1154) = {119, -55, -113, 112};\n\
Line Loop(1156) = {113, -54, -109, 108};\n\
Line Loop(1158) = {109, -53, -169, -110};\n\
Line Loop(1160) = {110, 168, -16, -114};\n\
Line Loop(1162) = {114, -15, -106, 104};\n\
Line Loop(1164) = {14, -106, -102, 101};\n\
Line Loop(1166) = {101, -13, -167, -100};\n\
Line Loop(1168) = {99, -45, -166, -100};\n\
Line Loop(1170) = {99, 46, -129, -98};\n\
Line Loop(1172) = {47, -125, 124, 129};\n\
Line Loop(1174) = {125, 48, -140, -126};\n\
Line Loop(1176) = {126, 139, 21, -127};\n\
Line Loop(1178) = {22, -122, 121, 127};\n\
Line Loop(1180) = {122, 23, -117, 116};\n\
Line Loop(1182) = {117, 24, -150, -118};\n\
Line Loop(1184) = {115, -73, -120, -116};\n\
Line Loop(1186) = {112, 115, 74, -111};\n\
Line Loop(1188) = {108, 111, 75, -107};\n\
Line Loop(1190) = {104, 107, 76, -103};\n\
Line Loop(1192) = {102, 103, 77, -97};\n\
Line Loop(1194) = {98, 128, -78, -97};\n\
Line Loop(1196) = {128, 79, -123, 124};\n\
Line Loop(1198) = {123, 80, -120, 121};\n\
')

# Now, define the curved surfaces (do a loop here)
for i in range(96):
    curInd = 1008+i*2
    fid.write('Ruled Surface(%d) = {%d};\n'%(curInd,curInd))

# Define the volume
fid.write('Surface Loop(1301) = {')
for i in range(6):
    fid.write('%d, '%(1001+i))
for i in range(95):
    curInd = 1008+i*2
    fid.write('%d, '%curInd)
fid.write('%d};\n'%(curInd+2))
fid.write('Volume(1401) = {1301};\n')

# Define the physical surfaces, periodic boundaries
fid.write('Physical Surface(3001) = {1001};\n' )
fid.write('Physical Surface(3002) = {1002};\n' )
fid.write('Physical Surface(3003) = {1003};\n' )
fid.write('Physical Surface(3004) = {1004};\n' )
fid.write('Physical Surface(3005) = {1005};\n' )
fid.write('Physical Surface(3006) = {1006};\n' )

# Boundaries with solid
fid.write('Physical Surface(3007) = {')
for i in range(95):
    curInd = 1008+i*2
    fid.write('%d, '%(curInd))
fid.write('%d};\n'%(curInd+2))

# Define the physical volume
fid.write('Physical Volume(4001) = {1401};\n\n' )

# Define the periodic surfaces, direction vector to translate the right surface for matching
# In this case, left is top surface, right is bottom surface (must be raised by one unit)
fid.write('Periodic Surface {1002} = {1001} Translate {1,0,0};\n');
fid.write('Periodic Surface {1004} = {1003} Translate {0,1,0};\n');
fid.write('Periodic Surface {1006} = {1005} Translate {0,0,1};\n');

# Set the mesh resolution
fid.write('Mesh.CharacteristicLengthMax = %.4f;\n'%resM2)

fid.close()
# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------------------------
# Construct interface cell, bottom part
# -----------------------------------------------------------------------------------------------------
GeoFile = '3D_MSE_bc_gen_geomB.geo'
# Determine the z coordinates
zmin = -dR-0.5 - (Nbot-1.0)
zmax =     0.0

# Compute the amount of points {support} + {outer circles} + {inner circles} + {mid-sphere-points} + {cell-corners}
npnt  = 13 + 6*8 + 6*8 + 2*8+4 + 8
Coord = np.zeros((npnt,4))
# Write the smallest mesh spacing by default
for i in range(npnt):
    Coord[i,3] = resM;
# Generate support points for curved segments (modify only non-zero components)
Coord[1,0] = -0.5; Coord[2,0] = -dR; Coord[3,0] = dR; Coord[4,0] = 0.5;
Coord[5,1] = -0.5; Coord[6,1] = -dR; Coord[7,1] = dR; Coord[8,1] = 0.5;
Coord[9,2] = -0.5; Coord[10,2]= -dR; Coord[11,2]= dR; Coord[12,2]= 0.5;

# Generate points for side circles, first around x axis, then y, then z, first - number offset
npntOffs = 13
for i in range(4):
    for j in range(8):
        Coord[npntOffs+i*8+j,0] = Coord[i+1,0];
        Coord[npntOffs+i*8+j,1] = cR*np.cos(j*np.pi/4);
        Coord[npntOffs+i*8+j,2] = cR*np.sin(j*np.pi/4);
npntOffs = 13+32
for i in range(4):
    for j in range(8):
        Coord[npntOffs+i*8+j,0] = cR*np.cos(j*np.pi/4);
        Coord[npntOffs+i*8+j,1] = Coord[i+5,1];
        Coord[npntOffs+i*8+j,2] = cR*np.sin(j*np.pi/4);
npntOffs = 13+64
for i in range(4):
    for j in range(8):
        Coord[npntOffs+i*8+j,0] = cR*np.cos(j*np.pi/4);
        Coord[npntOffs+i*8+j,1] = cR*np.sin(j*np.pi/4);
        Coord[npntOffs+i*8+j,2] = Coord[i+9,2];

npntOffs = 13+96
for j in range(8):
    Coord[npntOffs+j,0] = hR*np.cos(np.pi/4)*np.cos(j*np.pi/4);
    Coord[npntOffs+j,1] = hR*np.cos(np.pi/4)*np.sin(j*np.pi/4);
    Coord[npntOffs+j,2] = hR*np.sin(np.pi/4);
npntOffs = 13+104
for j in range(4):
    Coord[npntOffs+j,0] = hR*np.cos(np.pi/4+j*np.pi/2);
    Coord[npntOffs+j,1] = hR*np.sin(np.pi/4+j*np.pi/2);
    Coord[npntOffs+j,2] = 0.0;
npntOffs = 13+108
for j in range(8):
    Coord[npntOffs+j,0] = hR*np.cos(np.pi/4)*np.cos(j*np.pi/4);
    Coord[npntOffs+j,1] = hR*np.cos(np.pi/4)*np.sin(j*np.pi/4);
    Coord[npntOffs+j,2] =-hR*np.sin(np.pi/4);

# Additional points at the corners of the unit cell
npntOffs = 13+116
for i in range(2):
    for j in range(2):
        for k in range(2):
            Coord[npntOffs+i*4+j*2+k,0] = -0.5 + i;
            Coord[npntOffs+i*4+j*2+k,1] = -0.5 + j;
            Coord[npntOffs+i*4+j*2+k,2] = -0.5 + k;
            Coord[npntOffs+i*4+j*2+k,3] = resM2;

# Shifting all the points, such that y=0 corresponds to end of sphere
for i in range(npnt):
    Coord[i,2] = Coord[i,2] - dR

# Exporting all the points
fid = open(GeoFile,'w')
for i in range(npnt):
    fid.write('Point(%d) = {%.10f, %.10f, %.10f, %.10f};\n'%(i+1,Coord[i,0],Coord[i,1],Coord[i,2],Coord[i,3]))

# Below here, some operations are copied manually from GUI, but loops are used if found useful
# Generate the connecting circle arcs around x, y and z axis
offsArr = np.array([0, 8, 16, 24, 24, 28, 36, 44, 48, 52, 60, 68])
for i in range(12):
    for j in range(7):
        curLid = i*8 +j+1
        curPid = 13+curLid
        fid.write('Circle(%d) = {%d, %d, %d};\n'%( curLid, curPid , i+2 , curPid+1))
    fid.write('Circle(%d) = {%d, %d, %d};\n'%( curLid+1, curPid+1 , i+2 , curPid-6))

# Write connecting arcs of the sphere (copied from GUI)
fid.write('\
Circle(97) = {127, 1, 91};\n\
Circle(98) = {127, 1, 128};\n\
Circle(99) = {127, 1, 59};\n\
Circle(100) = {127, 1, 120};\n\
Circle(101) = {127, 1, 27};\n\
Circle(102) = {127, 1, 126};\n\
Circle(103) = {126, 1, 90};\n\
Circle(104) = {126, 1, 125};\n\
Circle(106) = {126, 1, 28};\n\
Circle(107) = {125, 1, 89};\n\
Circle(108) = {125, 1, 124};\n\
Circle(109) = {125, 1, 67};\n\
Circle(110) = {125, 1, 119};\n\
Circle(111) = {124, 1, 88};\n\
Circle(112) = {124, 1, 123};\n\
Circle(113) = {124, 1, 68};\n\
Circle(114) = {125, 1, 29};\n\
Circle(115) = {123, 1, 87};\n\
Circle(116) = {123, 1, 122};\n\
Circle(117) = {123, 1, 37};\n\
Circle(118) = {123, 1, 118};\n\
Circle(119) = {123, 1, 69};\n\
Circle(120) = {122, 1, 86};\n\
Circle(121) = {122, 1, 129};\n\
Circle(122) = {122, 1, 36};\n\
Circle(123) = {129, 1, 93};\n\
Circle(124) = {129, 1, 128};\n\
Circle(125) = {129, 1, 61};\n\
Circle(126) = {129, 1, 121};\n\
Circle(127) = {129, 1, 35};\n\
Circle(128) = {128, 1, 92};\n\
Circle(129) = {128, 1, 60};\n\
Circle(130) = {116, 1, 100};\n\
Circle(131) = {116, 1, 117};\n\
Circle(132) = {116, 1, 115};\n\
Circle(133) = {116, 1, 56};\n\
Circle(134) = {117, 1, 101};\n\
Circle(135) = {117, 1, 110};\n\
Circle(136) = {117, 1, 33};\n\
Circle(137) = {117, 1, 121};\n\
Circle(138) = {117, 1, 55};\n\
Circle(139) = {121, 1, 34};\n\
Circle(140) = {121, 1, 54};\n\
Circle(141) = {110, 1, 94};\n\
Circle(142) = {110, 1, 111};\n\
Circle(143) = {110, 1, 32};\n\
Circle(144) = {111, 1, 95};\n\
Circle(145) = {111, 1, 112};\n\
Circle(146) = {111, 1, 63};\n\
Circle(147) = {111, 1, 118};\n\
Circle(148) = {111, 1, 31};\n\
Circle(149) = {118, 1, 62};\n\
Circle(150) = {118, 1, 30};\n\
Circle(151) = {112, 1, 96};\n\
Circle(152) = {112, 1, 64};\n\
Circle(153) = {112, 1, 113};\n\
Circle(154) = {113, 1, 97};\n\
Circle(155) = {113, 1, 114};\n\
Circle(156) = {113, 1, 23};\n\
Circle(157) = {113, 1, 119};\n\
Circle(158) = {113, 1, 65};\n\
Circle(159) = {114, 1, 98};\n\
Circle(160) = {114, 1, 24};\n\
Circle(161) = {114, 1, 115};\n\
Circle(162) = {115, 1, 99};\n\
Circle(163) = {115, 1, 57};\n\
Circle(164) = {115, 1, 25};\n\
Circle(165) = {115, 1, 120};\n\
Circle(166) = {120, 1, 58};\n\
Circle(167) = {120, 1, 26};\n\
Circle(168) = {119, 1, 22};\n\
Circle(169) = {119, 1, 66};\n\
')

# Add the lines for cylinders
for i in range(6):
    nbndOffs = 14+i*16
    npntOffs = 170+i*8
    for j in range(8):
        fid.write('Line(%d) = {%d, %d};\n'%( npntOffs+j, nbndOffs+j , nbndOffs+j+8))

# Add closing surfaces for the fluid problem, define lines, close the surface
# Copied from the GUI output
# Line definitions
fid.write('\
Line(218) = {133, 132};\n\
Line(219) = {132, 130};\n\
Line(220) = {130, 131};\n\
Line(221) = {131, 133};\n\
Line(222) = {133, 137};\n\
Line(223) = {137, 135};\n\
Line(224) = {135, 131};\n\
Line(225) = {137, 136};\n\
Line(226) = {136, 132};\n\
Line(227) = {136, 134};\n\
Line(228) = {134, 130};\n\
Line(229) = {134, 135};\n\
')
# Line loop definitions
fid.write('\
Line Loop(1001) = {218, 219, 220, 221, -2, -1, -8, -7, -6, -5, -4, -3};\n\
Line Loop(1002) = {225, 227, 229, -223, -29, -28, -27, -26, -25, -32, -31, -30};\n\
Line Loop(1003) = {229, 224, -220, -228, -40, -39, -38, -37, -36, -35, -34, -33};\n\
Line Loop(1004) = {225, 226, -218, 222, -57, -64, -63, -62, -61, -60, -59, -58};\n\
Line Loop(1005) = {226, 219, -228, -227, -66, -65, -72, -71, -70, -69, -68, -67};\n\
')
# Define plane surfaces
for i in range(5):
    fid.write('Plane Surface(%d) = {%d};\n'%(1001+i,1001+i))

# Rest of surfaces are taken from GUI
# First we do line loops, then surfaces
fid.write('\
Line Loop(1008) = {6, 176, -14, -175};\n\
Line Loop(1010) = {5, 175, -13, -174};\n\
Line Loop(1012) = {174, -12, -173, 4};\n\
Line Loop(1014) = {3, 173, -11, -172};\n\
Line Loop(1016) = {2, 172, -10, -171};\n\
Line Loop(1018) = {1, 171, -9, -170};\n\
Line Loop(1020) = {177, 16, -170, -8};\n\
Line Loop(1022) = {7, 177, -15, -176};\n\
Line Loop(1024) = {183, 30, -184, -22};\n\
Line Loop(1026) = {184, 31, -185, -23};\n\
Line Loop(1028) = {185, 32, -178, -24};\n\
Line Loop(1030) = {178, 25, -179, -17};\n\
Line Loop(1032) = {179, 26, -180, -18};\n\
Line Loop(1034) = {180, 27, -181, -19};\n\
Line Loop(1036) = {181, 28, -182, -20};\n\
Line Loop(1038) = {183, -29, -182, 21};\n\
Line Loop(1040) = {192, 47, -193, -39};\n\
Line Loop(1042) = {193, 48, -186, -40};\n\
Line Loop(1044) = {33, 187, -41, -186};\n\
Line Loop(1046) = {34, 188, -42, -187};\n\
Line Loop(1048) = {35, 189, -43, -188};\n\
Line Loop(1050) = {189, 44, -190, -36};\n\
Line Loop(1052) = {190, 45, -191, -37};\n\
Line Loop(1054) = {191, 46, -192, -38};\n\
Line Loop(1056) = {60, -198, -52, 197};\n\
Line Loop(1058) = {198, 61, -199, -53};\n\
Line Loop(1060) = {62, -200, -54, 199};\n\
Line Loop(1062) = {55, 201, -63, -200};\n\
Line Loop(1064) = {56, 194, -64, -201};\n\
Line Loop(1066) = {194, 57, -195, -49};\n\
Line Loop(1068) = {195, 58, -196, -50};\n\
Line Loop(1070) = {196, 59, -197, -51};\n\
Line Loop(1072) = {67, 205, -75, -204};\n\
Line Loop(1074) = {204, -74, -203, 66};\n\
Line Loop(1076) = {76, -206, -68, 205};\n\
Line Loop(1078) = {77, -207, -69, 206};\n\
Line Loop(1080) = {207, 78, -208, -70};\n\
Line Loop(1082) = {208, 79, -209, -71};\n\
Line Loop(1084) = {209, 80, -202, -72};\n\
Line Loop(1086) = {202, 73, -203, -65};\n\
Line Loop(1088) = {212, 91, -213, -83};\n\
Line Loop(1090) = {213, 92, -214, -84};\n\
Line Loop(1092) = {93, -215, -85, 214};\n\
Line Loop(1094) = {215, 94, -216, -86};\n\
Line Loop(1096) = {216, 95, -217, -87};\n\
Line Loop(1098) = {217, 96, -210, -88};\n\
Line Loop(1100) = {210, 89, -211, -81};\n\
Line Loop(1102) = {212, -90, -211, 82};\n\n\
Line Loop(1104) = {84, -159, -155, 154};\n\
Line Loop(1106) = {159, 85, -162, -161};\n\
Line Loop(1108) = {162, 86, -130, 132};\n\
Line Loop(1110) = {154, -83, -151, 153};\n\
Line Loop(1112) = {130, 87, -134, -131};\n\
Line Loop(1114) = {134, 88, -141, -135};\n\
Line Loop(1116) = {142, 144, -81, -141};\n\
Line Loop(1118) = {144, 82, -151, -145};\n\
Line Loop(1120) = {143, 19, -136, 135};\n\
Line Loop(1122) = {136, 20, -139, -137};\n\
Line Loop(1124) = {137, 140, 41, -138};\n\
Line Loop(1126) = {138, 42, -133, 131};\n\
Line Loop(1128) = {133, 43, -163, -132};\n\
Line Loop(1130) = {163, 44, -166, -165};\n\
Line Loop(1132) = {165, 167, -12, -164};\n\
Line Loop(1134) = {164, -11, -160, 161};\n\
Line Loop(1136) = {160, -10, -156, 155};\n\
Line Loop(1138) = {156, -9, -168, -157};\n\
Line Loop(1140) = {157, 169, -52, -158};\n\
Line Loop(1142) = {158, -51, -152, 153};\n\
Line Loop(1144) = {152, -50, -146, 145};\n\
Line Loop(1146) = {146, -49, -149, -147};\n\
Line Loop(1148) = {147, 150, 17, -148};\n\
Line Loop(1150) = {148, 18, -143, 142};\n\
Line Loop(1152) = {149, -56, -119, 118};\n\
Line Loop(1154) = {119, -55, -113, 112};\n\
Line Loop(1156) = {113, -54, -109, 108};\n\
Line Loop(1158) = {109, -53, -169, -110};\n\
Line Loop(1160) = {110, 168, -16, -114};\n\
Line Loop(1162) = {114, -15, -106, 104};\n\
Line Loop(1164) = {14, -106, -102, 101};\n\
Line Loop(1166) = {101, -13, -167, -100};\n\
Line Loop(1168) = {99, -45, -166, -100};\n\
Line Loop(1170) = {99, 46, -129, -98};\n\
Line Loop(1172) = {47, -125, 124, 129};\n\
Line Loop(1174) = {125, 48, -140, -126};\n\
Line Loop(1176) = {126, 139, 21, -127};\n\
Line Loop(1178) = {22, -122, 121, 127};\n\
Line Loop(1180) = {122, 23, -117, 116};\n\
Line Loop(1182) = {117, 24, -150, -118};\n\
Line Loop(1184) = {115, -73, -120, -116};\n\
Line Loop(1186) = {112, 115, 74, -111};\n\
Line Loop(1188) = {108, 111, 75, -107};\n\
Line Loop(1190) = {104, 107, 76, -103};\n\
Line Loop(1192) = {102, 103, 77, -97};\n\
Line Loop(1194) = {98, 128, -78, -97};\n\
Line Loop(1196) = {128, 79, -123, 124};\n\
Line Loop(1198) = {123, 80, -120, 121};\n\
')

# Now, define the curved surfaces (do a loop here)
for i in range(96):
    curInd = 1008+i*2
    fid.write('Ruled Surface(%d) = {%d};\n'%(curInd,curInd))

# Do extension of unit cell down to the bottom of the coating
for i in range(Nbot-1):
    fid.write('FreeSurface%d[] = Translate {0,0,-%d} {Duplicata{ Surface{'%(i+1,i+1) )
    for j in range(95):
        curInd = 1008+j*2
        fid.write('%d, '%(curInd))
    fid.write('%d}; }};\n'%(curInd+2))
for i in range(Nbot-1):
    for j in range(4):
        fid.write('P%dSurface%d[] = Translate {0,0,-%d} {Duplicata{ Surface{%d}; }};\n'%(1001+j,i+1,i+1,1001+j))

# Translate the closing surface
fid.write('BSurface[] = Translate {0,0,-%d} {Duplicata{ Surface{1005}; }};\n'%(Nbot-1))
fid.write('Delete { Surface{1005}; }\n')


# Remove the top cylinder, Surfaces, Lines, Points
for i in range(8):
    fid.write('Delete { Surface{%d}; }\n'%(1088+i*2))
for i in range(8):
    fid.write('Delete { Line{%d}; }\n'%(89+i))
for i in range(8):
    fid.write('Delete { Line{%d}; }\n'%(210+i))
fid.write('Delete { Point{13}; }\n')
for i in range(8):
    fid.write('Delete { Point{%d}; }\n'%(102+i))

# Move top corner points
fid.write('c[] = Point{131};\n')
fid.write('Translate {0,0,-c[2]} { Point{%d}; }\n'%(131))
fid.write('Translate {0,0,-c[2]} { Point{%d}; }\n'%(133))
fid.write('Translate {0,0,-c[2]} { Point{%d}; }\n'%(135))
fid.write('Translate {0,0,-c[2]} { Point{%d}; }\n'%(137))

# Now close the top fluid surface
fid.write('Line Loop(1006) = {223, 224, 221, 222, -81,-82,-83,-84,-85,-86,-87,-88};\n')
fid.write('Plane Surface(1006) = {1006};\n')

# Set finer mesh spacing at the interface
fid.write('Characteristic Length {131, 133, 135, 137} = %.10f;\n'%resM)
    
# Define the volume, add the translated surface entities
fid.write('Surface Loop(1301) = {')
for i in range(4):
    fid.write('%d, '%(1001+i))
fid.write('%d, '%(1006))
for i in range(96):
    curInd = 1008+i*2
    if (curInd>=1088 and curInd<=1102):
        continue
    fid.write('%d, '%curInd)
for i in range(Nbot-1):
    for j in range(4):
        fid.write('P%dSurface%d[], '%(1001+j,i+1))
for i in range(Nbot-1):
    fid.write('FreeSurface%d[], '%(i+1))
fid.write('BSurface[]};\n')
fid.write('Volume(1401) = {1301};\n')

# Define the physical surfaces
# Periodic boundaries, including the translated ones (need to have separate
# pysical labels, because they are disconnected)
fid.write('Physical Surface(3001) = {1001};\n')
for i in range(Nbot-1):
    fid.write('Physical Surface(3001%d) = {P1001Surface%d[]};\n'%(i+1,i+1))
fid.write('Physical Surface(3002) = {1002};\n')
for i in range(Nbot-1):
    fid.write('Physical Surface(3002%d) = {P1002Surface%d[]};\n'%(i+1,i+1))
fid.write('Physical Surface(3003) = {1003};\n')
for i in range(Nbot-1):
    fid.write('Physical Surface(3003%d) = {P1003Surface%d[]};\n'%(i+1,i+1))
fid.write('Physical Surface(3004) = {1004};\n')
for i in range(Nbot-1):
    fid.write('Physical Surface(3004%d) = {P1004Surface%d[]};\n'%(i+1,i+1))
fid.write('Physical Surface(3005) = {BSurface[]};\n' )
fid.write('Physical Surface(3006) = {1006};\n' )

# Free boundaries, including the translated ones
fid.write('Physical Surface(3007) = {')
for i in range(96):
    curInd = 1008+i*2
    if curInd>=1088 and curInd<=1102:
        continue
    fid.write('%d, '%(curInd))
for i in range(Nbot-2):
    fid.write('FreeSurface%d[], '%(i+1))
fid.write('FreeSurface%d[]};\n'%(Nbot-1))

# Define the physical volume
fid.write('Physical Volume(4001) = {1401};\n\n' )

# Define the periodic surfaces, direction vector to translate the right surface for matching
# In this case, left is top surface, right is bottom surface (must be raised by one unit)
fid.write('Periodic Surface {1002} = {1001} Translate {1,0,0};\n');
fid.write('Periodic Surface {1004} = {1003} Translate {0,1,0};\n');

# Translated periodic surfaces
for i in range(Nbot-1):
    fid.write('Periodic Surface {P1002Surface%d[]} = {P1001Surface%d[]} Translate {1,0,0};\n'%(i+1,i+1));
    fid.write('Periodic Surface {P1004Surface%d[]} = {P1003Surface%d[]} Translate {0,1,0};\n'%(i+1,i+1));

# Set the mesh resolution
fid.write('Mesh.CharacteristicLengthMax = %.4f;\n'%(resM2))

fid.close()
# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------------------------
# Construct interface cell, top part
# -----------------------------------------------------------------------------------------------------
GeoFile = '3D_MSE_bc_gen_geomT.geo'
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
fid = open(GeoFile,'w')
for i in range(npnt):
    fid.write('Point(%d) = {%.10f, %.10f, %.10f, %.10f};\n'%(i+1,Coord[i,0],Coord[i,1],Coord[i,2],Coord[i,3]))

# Define lines
fid.write('\
Circle(1) = {14, 13, 15};\n\
Circle(2) = {15, 13, 14};\n\
')
fid.write('\
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
fid.write('\
Line Loop(24) = {2, 1};\n\
Line Loop(31) = {10, 7, 8, 9, -2, -1};\n\
Plane Surface(24) = {24};\n\
')
fid.write('\
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
fid.write('Surface Loop(43) = {41, 33, 31, 37, 35, 39, 24};\n')
fid.write('Volume(43) = {43};\n')

# Define the physical surfaces
# Periodic boundaries, including the translated ones (need to have separate
# pysical labels, because they are disconnected)
fid.write('Physical Surface(3001) = {37};\n')        # (y,z)-plane, x-
fid.write('Physical Surface(3002) = {33};\n')        # (y,z)-plane, x+
fid.write('Physical Surface(3003) = {35};\n')        # (x,z)-plane, y-
fid.write('Physical Surface(3004) = {39};\n')        # (x,z)-plane, y+
fid.write('Physical Surface(3005) = {31};\n')        # Bottom
fid.write('Physical Surface(3006) = {41};\n')        # Top
fid.write('Physical Surface(3007) = {24};\n' )       # No-slip walls

# Define the physical volume
fid.write('Physical Volume(4001) = {43};\n\n' )

# Define the periodic surfaces, direction vector to translate the right surface for matching
# In this case, left is top surface, right is bottom surface (must be raised by one unit)
fid.write('Periodic Surface {33} = {37} Translate {1,0,0};\n');
fid.write('Periodic Surface {39} = {35} Translate {0,1,0};\n');

# Set the mesh resolution
fid.write('Mesh.CharacteristicLengthMax = %.4f;\n'%(resM2))

fid.close()


# -----------------------------------------------------------------------------------------------------
# Construct Lagrange multiplier space
# -----------------------------------------------------------------------------------------------------
GeoFile = '3D_MSE_bc_gen_geomL.geo'
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
fid = open(GeoFile,'w')
for i in range(npnt):
    fid.write('Point(%d) = {%.10f, %.10f, %.10f, %.10f};\n'%(i+1,Coord[i,0],Coord[i,1],Coord[i,2],Coord[i,3]))

# Define lines
fid.write('\
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
fid.write('Circle(13) = {10, 9, 11};\n')
fid.write('Circle(14) = {11, 9, 10};\n')
# Define surfaces
fid.write('\
Line Loop(19) = {13, 14};\n\
Line Loop(20) = {13, 14, 1,-2,-3,4};\n\
Plane Surface(19) = {19};\n\
')
fid.write('\
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
fid.write('Surface Loop(43) = {20,21,22,23,24,25,19};\n')
fid.write('Volume(43) = {43};\n')

# Define the physical surfaces
# Periodic boundaries, including the translated ones (need to have separate
# pysical labels, because they are disconnected)
fid.write('Physical Surface(3001) = {22};\n')        # (y,z)-plane, x-
fid.write('Physical Surface(3002) = {24};\n')        # (y,z)-plane, x+
fid.write('Physical Surface(3003) = {23};\n')        # (x,z)-plane, y-
fid.write('Physical Surface(3004) = {25};\n')        # (x,z)-plane, y+
fid.write('Physical Surface(3005) = {20};\n')        # Bottom
fid.write('Physical Surface(3006) = {21};\n')        # Top
fid.write('Physical Surface(3007) = {19};\n')        # Interior walls, not needed for Lagrange multipliers
# Define the physical volume
fid.write('Physical Volume(4001) = {43};\n\n' )
# Define the periodic surfaces, direction vector to translate the right surface for matching
# In this case, left is top surface, right is bottom surface (must be raised by one unit)
fid.write('Periodic Surface {24} = {22} Translate {1,0,0};\n');
fid.write('Periodic Surface {25} = {23} Translate {0,1,0};\n');

# Set the mesh resolution
fid.write('Mesh.CharacteristicLengthMax = %.4f;\n'%(resM2))

fid.close()
# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------
