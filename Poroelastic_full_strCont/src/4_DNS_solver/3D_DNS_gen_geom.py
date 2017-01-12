#!/usr/bin/python2.7
# -----------------------------------------------------------------------------------------------------
# Phython script to generate geo file for mesh generation.
# Boundary points are generated in order to match the required mesh.
# Note that for curved surfaces - up to 4 boundaries are allowed!!!
#
# Author:
# Ugis Lacis
# ugis@mech.kth.se
# 10.04.2015
# -----------------------------------------------------------------------------------------------------

import sys
import numpy as np
import scipy.optimize as opt
import shutil as sh

# Read the case parameters
try:
    InputParam = np.loadtxt("3D_DNS_parameters.in")
except IOError:
    sys.exit("Error: expecting parameter file <3D_DNS_parameters.in>")
nCellX = int(InputParam[0])
nCell  = int(InputParam[1])
eA     =     InputParam[2] # x axis
eB     =     InputParam[3] # y axis
eC     =     InputParam[4] # z axis
cR     =     InputParam[5]
ePhi   =     InputParam[6]
resM   =     InputParam[7]
resM2  =     InputParam[8]

# Check, how the axis are progressing? Which one is the biggest one, which one is the smallest one? Or not...
# Just check later on.

# Define epsilon number (should set such that 1/epsP - integer)
epsP = 1.0/nCellX

# Display used parameters
print "Finished reading parameters, nCellX = "+str(nCellX)+", nCell = "+str(nCell)+", resulting in eps = "+str(epsP)

# Vertical coordinate
zmax = 1.0/epsP
zmin = - nCell
# Side coordinates
xmin   =-0.5/epsP
xmax   = 0.5/epsP

# Save the obtained coordinates to a new parameter file
try:
    InputParam[11] = zmin
except IndexError:
    InputParam = np.append(InputParam, zmin)
try:
    InputParam[12] = zmax
except IndexError:
    InputParam = np.append(InputParam, zmax)
try:
    InputParam[13] = xmin
except IndexError:
    InputParam = np.append(InputParam, xmin)
try:
    InputParam[14] = xmax
except IndexError:
    InputParam = np.append(InputParam, xmax)
try:
    InputParam[15] = epsP
except IndexError:
    InputParam = np.append(InputParam, epsP)
np.savetxt("3D_DNS_parameters.txt",InputParam)
# Append the comment line at the end
fid = open("3D_DNS_parameters.txt", "a")
fid.write("# NbotX, NbotD, eA, eB, eC, cR, phi[deg], dsMin, dsMax, tolSOL, barE, zmin, zmax, xmin, xmax, epsP\n")
fid.close()

# Number of B-spline points for 1/4th of curved circle
Nbsp = 6

# Output filenames
GeoFileS = '3D_DNS_gen_geomS.geo'
GeoFileF = '3D_DNS_gen_geomF.geo'

# Rescale angle to radians, sign for correct turn direction
ePhi = -ePhi/180.0*np.pi

# Derivation of intersections with cylinders
# ----------------------------------------------------------------------------------------------------------------
# Ellipsoid equation, aligned
# x^2/a^2 + y^2/b^2 + z^2/c^2 = 1
# Ellipsoid equation, rotated, cos and sin has argument ePhi
# (x*cos - z*sin)^2/a^2 + y^2/b^2 + (x*sin + z*cos)^2/c^2 = 1
# (x^2*cos^2+z^2*sin^2-2*x*z*sin*cos)/a^2 + y^2/b^2 + (x^2*sin^2+z^2*cos^2+2*x*z*sin*cos)/c^2 = 1
# For x,y,z-cylinders, solve quadratic equation for x,y,z-coordinates
# x^2*[cos^2/a^2+sin^2/c^2] + x*[2*z*sin*cos/c^2-2*z*sin*cos/a^2] + [z^2*sin^2/a^2+y^2/b^2+z^2*cos^2/c^2-1] = 0
# y^2*[1/b^2]  + y*0  + [(x^2*cos^2+z^2*sin^2-2*x*z*sin*cos)/a^2+(x^2*sin^2+z^2*cos^2+2*x*z*sin*cos)/c^2-1] = 0
# z^2*[sin^2/a^2+cos^2/c^2] + z*[2*x*sin*cos/c^2-2*x*sin*cos/a^2] + [x^2*cos^2/a^2+y^2/b^2+x^2*sin^2/c^2-1] = 0
# z = ( -B +- sqrt(B^2-4*A*C) )/(2*A)
# A, B, C obtained using x and y coordinates from z-cylinder, parametric
# x = cR*sin(t)
# y = cR*cos(t)
# 3D curve done!!!
# Now, coefficients for all slices:
# Ax = cos^2/a^2+sin^2/c^2; Bx = 2*z*sin*cos/c^2-2*z*sin*cos/a^2; Cx = z^2*sin^2/a^2+y^2/b^2+z^2*cos^2/c^2-1; y = cR*sin(t); z = cR*cos(t)
# Az = sin^2/a^2+cos^2/c^2; Bz = 2*x*sin*cos/c^2-2*x*sin*cos/a^2; Cz = x^2*cos^2/a^2+y^2/b^2+x^2*sin^2/c^2-1; x = cR*sin(t); y = cR*cos(t)
# Ay = 1/b^2; By = 0; Cy = (x^2*cos^2+z^2*sin^2-2*x*z*sin*cos)/a^2+(x^2*sin^2+z^2*cos^2+2*x*z*sin*cos)/c^2-1; x = cR*sin(t); z = cR*cos(t)
# ----------------------------------------------------------------------------------------------------------------

# Derivation of plane curves
# (x^2*cos^2+z^2*sin^2-2*x*z*sin*cos)/a^2 + y^2/b^2 + (x^2*sin^2+z^2*cos^2+2*x*z*sin*cos)/c^2 = 1
# Slice with x = 0
# (z^2*sin^2)/a^2 + y^2/b^2 + (z^2*cos^2)/c^2 = 1 --> y^2/b^2 + z^2*(sin^2/a^2+cos^2/c^2) = 1; Still ellipse!
# Slice with y = 0
# (x^2*cos^2+z^2*sin^2-2*x*z*sin*cos)/a^2 + (x^2*sin^2+z^2*cos^2+2*x*z*sin*cos)/c^2 = 1; Still ellipse, turned with angle!
# Slice with z = 0
# Something similar, with an ellipse

# Compute the count of points {center ellipse sup} + {border cyl} + {bspline cyl} + {mid-ellipse points} + {cell-corners} + {add ell support}
npnt  = 1+3+2+2+2 + 6*9 + 6*4*Nbsp + 3*4+2*4 + 8 + 2+2+2+2
# + 6*8 + 2*8+4 + 8
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

# -----------------------------------------------------------------------------------------------------
# Additional support points for turned ellipses
# -----------------------------------------------------------------------------------------------------
nTES = npnt - 8
# For planes z = +- y
# First, rotate around y-axis with turn angle ePhi2
# Then, rotate around x-axis with turn angle +- pi/4
cs45  = np.cos(np.pi/4)
ePhi2 = ePhi
Coord[nTES+0,0] = dsub*np.cos(ePhi2); Coord[nTES+0,1] =-dsub*np.sin(ePhi2)*cs45 ; Coord[nTES+0,2] =-dsub*np.sin(ePhi2)*cs45    # hat{x} coordinates after 2nd rotation
Coord[nTES+1,0] = dsub*np.sin(ePhi2); Coord[nTES+1,1] = dsub*np.cos(ePhi2)*cs45;  Coord[nTES+1,2] = dsub*np.cos(ePhi2)*cs45    # hat{z} coordinates after 2nd rotation
Coord[nTES+2,0] = dsub*np.cos(ePhi2); Coord[nTES+2,1] = dsub*np.sin(ePhi2)*cs45 ; Coord[nTES+2,2] =-dsub*np.sin(ePhi2)*cs45    # hat{x} coordinates after 2nd rotation
Coord[nTES+3,0] = dsub*np.sin(ePhi2); Coord[nTES+3,1] =-dsub*np.cos(ePhi2)*cs45;  Coord[nTES+3,2] = dsub*np.cos(ePhi2)*cs45    # hat{z} coordinates after 2nd rotation
# For planes y = +- x
# First, rotate around y-axis with turn angle ePhi2
# Then, rotate around z-axis with turn angle +- pi/4
Coord[nTES+4,0] = dsub*np.cos(ePhi2)*cs45; Coord[nTES+4,1] = dsub*np.cos(ePhi2)*cs45 ; Coord[nTES+4,2] =-dsub*np.sin(ePhi2)    # hat{x} coordinates after 2nd rotation
Coord[nTES+5,0] = dsub*np.sin(ePhi2)*cs45; Coord[nTES+5,1] = dsub*np.sin(ePhi2)*cs45;  Coord[nTES+5,2] = dsub*np.cos(ePhi2)    # hat{z} coordinates after 2nd rotation
Coord[nTES+6,0] = dsub*np.cos(ePhi2)*cs45; Coord[nTES+6,1] =-dsub*np.cos(ePhi2)*cs45 ; Coord[nTES+6,2] =-dsub*np.sin(ePhi2)    # hat{x} coordinates after 2nd rotation
Coord[nTES+7,0] = dsub*np.sin(ePhi2)*cs45; Coord[nTES+7,1] =-dsub*np.sin(ePhi2)*cs45;  Coord[nTES+7,2] = dsub*np.cos(ePhi2)    # hat{z} coordinates after 2nd rotation
# -----------------------------------------------------------------------------------------------------

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
# (x^2*cos^2+z^2*sin^2-2*x*z*sin*cos)/a^2 + y^2/b^2 + (x^2*sin^2+z^2*cos^2+2*x*z*sin*cos)/c^2 = 1
# z = 0
# x^2/[1/(cos^2/a^2+sin^2/c^2)] + y^2/b^2 = 1 --> x^2/Ax^2 + y^2/b^2 = 1
# x = Ax*sin(t); y = b*cos(t)
Ax  = np.sqrt(1.0/(cosP2/eA**2+sinP2/eC**2))
tpm = np.arctan(eB/Ax)
Coord[npntOffs+0,0] = Ax*np.sin(tpm); Coord[npntOffs+0,1] = eB*np.cos(tpm); Coord[npntOffs+0,2] = 0
Coord[npntOffs+1,0] =-Ax*np.sin(tpm); Coord[npntOffs+1,1] = eB*np.cos(tpm); Coord[npntOffs+1,2] = 0
Coord[npntOffs+2,0] = Ax*np.sin(tpm); Coord[npntOffs+2,1] =-eB*np.cos(tpm); Coord[npntOffs+2,2] = 0
Coord[npntOffs+3,0] =-Ax*np.sin(tpm); Coord[npntOffs+3,1] =-eB*np.cos(tpm); Coord[npntOffs+3,2] = 0
XYi = Ax*np.sin(tpm)

# (x^2*cos^2+z^2*sin^2-2*x*z*sin*cos)/a^2 + y^2/b^2 + (x^2*sin^2+z^2*cos^2+2*x*z*sin*cos)/c^2 = 1
# y = 0, x = +- z
# (x*cos-z*sin)^2/a^2 + (x*sin+z*cos)^2/c^2 = 1
# For x = z, we get
# x^2*[(cos-sin)^2/a^2 + (sin+cos)^2/c^2] = 1
# x = +- sqrt(1.0 / [(cos-sin)^2/a^2 + (sin+cos)^2/c^2] )
npntOffs = npntOffs+4
xzcur = np.sqrt(1.0/( (np.cos(ePhi)-np.sin(ePhi))**2/eA**2 + (np.sin(ePhi)+np.cos(ePhi))**2/eC**2 ))
Coord[npntOffs+0,0] = xzcur; Coord[npntOffs+0,1] = 0; Coord[npntOffs+0,2] = xzcur 
Coord[npntOffs+1,0] =-xzcur; Coord[npntOffs+1,1] = 0; Coord[npntOffs+1,2] =-xzcur 
# For x =-z, we get
# x^2*[(cos+sin)^2/a^2 + (sin-cos)^2/c^2] = 1
# x = +- sqrt(1.0 / [(cos+sin)^2/a^2 + (sin-cos)^2/c^2] )
xzcur = np.sqrt(1.0/( (np.cos(ePhi)+np.sin(ePhi))**2/eA**2 + (np.sin(ePhi)-np.cos(ePhi))**2/eC**2 ))
Coord[npntOffs+2,0] = xzcur; Coord[npntOffs+2,1] = 0; Coord[npntOffs+2,2] =-xzcur 
Coord[npntOffs+3,0] =-xzcur; Coord[npntOffs+3,1] = 0; Coord[npntOffs+3,2] = xzcur 

# (x^2*cos^2+z^2*sin^2-2*x*z*sin*cos)/a^2 + y^2/b^2 + (x^2*sin^2+z^2*cos^2+2*x*z*sin*cos)/c^2 = 1
# x = 0
# z^2/[1/(sin^2/a^2+cos^2/c^2)] + y^2/b^2 = 1 --> z^2/Az^2 + y^2/b^2 = 1
# z = Az*sin(t); y = b*cos(t)
npntOffs = npntOffs+4
Az  = np.sqrt(1.0/(sinP2/eA**2+cosP2/eC**2))
tpm = np.arctan(eB/Az);
Coord[npntOffs+0,0] = 0; Coord[npntOffs+0,1] = eB*np.cos(tpm); Coord[npntOffs+0,2] = Az*np.sin(tpm)
Coord[npntOffs+1,0] = 0; Coord[npntOffs+1,1] =-eB*np.cos(tpm); Coord[npntOffs+1,2] = Az*np.sin(tpm)
Coord[npntOffs+2,0] = 0; Coord[npntOffs+2,1] = eB*np.cos(tpm); Coord[npntOffs+2,2] =-Az*np.sin(tpm)
Coord[npntOffs+3,0] = 0; Coord[npntOffs+3,1] =-eB*np.cos(tpm); Coord[npntOffs+3,2] =-Az*np.sin(tpm)
YZi = Az*np.sin(tpm)

# Mark points of the connecting 45-degree arcs, under condition x=y=z
npntOffs = npntOffs+4
# (x*cos-z*sin)^2/a^2 + y^2/b^2 + (x*sin + z*cos)^2/c^2 = 1
# Solving for x, if x = z
# x^2*[(cos-sin)^2/a^2 + 1/b^2 + (sin + cos)^2/c^2] = 1
# x = +- sqrt( 1/[(cos-sin)^2/a^2 + 1/b^2 + (sin + cos)^2/c^2] )
curxyz = np.sqrt( 1.0/( (np.cos(ePhi)-np.sin(ePhi))**2/eA**2 + 1.0/eB**2 \
                      + (np.sin(ePhi)+np.cos(ePhi))**2/eC**2 ) )

Coord[npntOffs+0,0] = curxyz; Coord[npntOffs+0,1] = curxyz; Coord[npntOffs+0,2] = curxyz
Coord[npntOffs+2,0] = curxyz; Coord[npntOffs+2,1] =-curxyz; Coord[npntOffs+2,2] = curxyz
Coord[npntOffs+5,0] =-curxyz; Coord[npntOffs+5,1] = curxyz; Coord[npntOffs+5,2] =-curxyz
Coord[npntOffs+7,0] =-curxyz; Coord[npntOffs+7,1] =-curxyz; Coord[npntOffs+7,2] =-curxyz

# Solving for x, if x = -z
# x^2*[(cos+sin)^2/a^2 + 1/b^2 + (sin - cos)^2/c^2] = 1
# x = +- sqrt( 1/[(cos+sin)^2/a^2 + 1/b^2 + (sin-cos)^2/c^2] )
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

# Shifting all the points, to have z = 0 at the tip of the solid structure, and xmin at the left side of the cell
for i in range(npnt):
    Coord[i,0] = Coord[i,0] + xmin + 0.5
    Coord[i,2] = Coord[i,2] - 0.5

# Exporting all the points
fidF = open(GeoFileF,'w')
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
Ruled Surface(230) = {230};\n\
Line Loop(232) = {22, 105, 75, -112};\n\
Ruled Surface(232) = {232};\n\
Line Loop(234) = {23, 106, 74, -105};\n\
Ruled Surface(234) = {234};\n\
Line Loop(236) = {106, -73, -107, -24};\n\
Ruled Surface(236) = {236};\n\
Line Loop(238) = {80, -107, 17, 108};\n\
Ruled Surface(238) = {238};\n\
Line Loop(240) = {108, -79, -109, -18};\n\
Ruled Surface(240) = {240};\n\
Line Loop(242) = {109, -78, -110, -19};\n\
Ruled Surface(242) = {242};\n\
Line Loop(244) = {110, -77, -111, -20};\n\
Ruled Surface(244) = {244};\n\
Line Loop(246) = {66, -99, 31, 98};\n\
Ruled Surface(246) = {246};\n\
Line Loop(248) = {32, 97, 65, -98};\n\
Ruled Surface(248) = {248};\n\
Line Loop(250) = {97, -72, -104, -25};\n\
Ruled Surface(250) = {250};\n\
Line Loop(252) = {104, -71, -103, -26};\n\
Ruled Surface(252) = {252};\n\
Line Loop(254) = {103, -70, -102, -27};\n\
Ruled Surface(254) = {254};\n\
Line Loop(256) = {28, 101, 69, -102};\n\
Ruled Surface(256) = {256};\n\
Line Loop(258) = {29, 100, 68, -101};\n\
Ruled Surface(258) = {258};\n\
Line Loop(260) = {30, 99, 67, -100};\n\
Ruled Surface(260) = {260};\n\
Line Loop(262) = {4, 127, -52, -128};\n\
Ruled Surface(262) = {262};\n\
Line Loop(264) = {5, 126, -53, -127};\n\
Ruled Surface(264) = {264};\n\
Line Loop(266) = {126, 54, -125, -6};\n\
Ruled Surface(266) = {266};\n\
Line Loop(268) = {125, 55, -124, -7};\n\
Ruled Surface(268) = {268};\n\
Line Loop(270) = {124, 56, -123, -8};\n\
Ruled Surface(270) = {270};\n\
Line Loop(272) = {1, 122, -49, -123};\n\
Ruled Surface(272) = {272};\n\
Line Loop(274) = {122, 50, -121, -2};\n\
Ruled Surface(274) = {274};\n\
Line Loop(276) = {3, 128, -51, -121};\n\
Ruled Surface(276) = {276};\n\
Line Loop(278) = {58, -113, 12, 114};\n\
Ruled Surface(278) = {278};\n\
Line Loop(280) = {59, -120, 11, 113};\n\
Ruled Surface(280) = {280};\n\
Line Loop(282) = {120, 60, 119, 10};\n\
Ruled Surface(282) = {282};\n\
Line Loop(284) = {119, -9, 118, -61};\n\
Ruled Surface(284) = {284};\n\
Line Loop(286) = {118, 62, -117, 16};\n\
Ruled Surface(286) = {286};\n\
Line Loop(288) = {15, 117, 63, -116};\n\
Ruled Surface(288) = {288};\n\
Line Loop(290) = {14, 116, 64, -115};\n\
Ruled Surface(290) = {290};\n\
Line Loop(292) = {13, 115, 57, -114};\n\
Ruled Surface(292) = {292};\n\
Line Loop(294) = {144, 84, -143, -35};\n\
Ruled Surface(294) = {294};\n\
Line Loop(296) = {143, 85, -142, -36};\n\
Ruled Surface(296) = {296};\n\
Line Loop(298) = {86, -141, -37, 142};\n\
Ruled Surface(298) = {298};\n\
Line Loop(300) = {141, 87, -140, -38};\n\
Ruled Surface(300) = {300};\n\
Line Loop(302) = {140, 88, -139, -39};\n\
Ruled Surface(302) = {302};\n\
Line Loop(304) = {139, 81, -138, -40};\n\
Ruled Surface(304) = {304};\n\
Line Loop(306) = {33, 137, -82, -138};\n\
Ruled Surface(306) = {306};\n\
Line Loop(308) = {34, 144, -83, -137};\n\
Ruled Surface(308) = {308};\n\
Line Loop(310) = {90, -129, 43, 130};\n\
Ruled Surface(310) = {310};\n\
Line Loop(312) = {130, -89, -131, -44};\n\
Ruled Surface(312) = {312};\n\
Line Loop(314) = {131, -96, -132, -45};\n\
Ruled Surface(314) = {314};\n\
Line Loop(316) = {132, -95, -133, -46};\n\
Ruled Surface(316) = {316};\n\
Line Loop(318) = {94, -133, 47, 134};\n\
Ruled Surface(318) = {318};\n\
Line Loop(320) = {134, -93, -135, -48};\n\
Ruled Surface(320) = {320};\n\
Line Loop(322) = {41, 136, 92, -135};\n\
Ruled Surface(322) = {322};\n\
Line Loop(324) = {136, -91, -129, -42};\n\
Ruled Surface(324) = {324};\n\
Line Loop(326) = {155, 179, -193, 90};\n\
Ruled Surface(326) = {326};\n\
Line Loop(328) = {193, -204, -163, 89};\n\
Ruled Surface(328) = {328};\n\
Line Loop(330) = {179, 203, -67, -156};\n\
Ruled Surface(330) = {330};\n\
Line Loop(332) = {68, -150, -194, 203};\n\
Ruled Surface(332) = {332};\n\
Line Loop(334) = {194, -149, 58, -180};\n\
Ruled Surface(334) = {334};\n\
Line Loop(336) = {180, 59, -164, 204};\n\
Ruled Surface(336) = {336};\n\
Line Loop(338) = {164, 60, -172, 202};\n\
Ruled Surface(338) = {338};\n\
Line Loop(340) = {163, -202, -185, 96};\n\
Ruled Surface(340) = {340};\n\
Line Loop(342) = {185, -171, 154, 95};\n\
Ruled Surface(342) = {342};\n\
Line Loop(344) = {153, 171, -201, 74};\n\
Ruled Surface(344) = {344};\n\
Line Loop(346) = {201, 186, -147, 73};\n\
Ruled Surface(346) = {346};\n\
Line Loop(348) = {186, 148, -61, -172};\n\
Ruled Surface(348) = {348};\n\
Line Loop(350) = {155, -178, 192, -91};\n\
Ruled Surface(350) = {350};\n\
Line Loop(352) = {156, -66, 209, 178};\n\
Ruled Surface(352) = {352};\n\
Line Loop(354) = {209, -191, -151, 65};\n\
Ruled Surface(354) = {354};\n\
Line Loop(356) = {152, -52, 177, -191};\n\
Ruled Surface(356) = {356};\n\
Line Loop(358) = {177, 210, -161, 51};\n\
Ruled Surface(358) = {358};\n\
Line Loop(360) = {210, 162, -92, -192};\n\
Ruled Surface(360) = {360};\n\
Line Loop(362) = {162, 93, -200, -211};\n\
Ruled Surface(362) = {362};\n\
Line Loop(364) = {211, -169, 50, 161};\n\
Ruled Surface(364) = {364};\n\
Line Loop(366) = {169, -199, -145, 49};\n\
Ruled Surface(366) = {366};\n\
Line Loop(368) = {170, 154, -94, -200};\n\
Ruled Surface(368) = {368};\n\
Line Loop(370) = {170, -153, 75, -212};\n\
Ruled Surface(370) = {370};\n\
Line Loop(372) = {212, 76, -146, 199};\n\
Ruled Surface(372) = {372};\n\
Line Loop(374) = {158, 83, -196, 174};\n\
Ruled Surface(374) = {374};\n\
Line Loop(376) = {174, -157, -70, -216};\n\
Ruled Surface(376) = {376};\n\
Line Loop(378) = {195, 216, -69, -150};\n\
Ruled Surface(378) = {378};\n\
Line Loop(380) = {149, 195, -173, 57};\n\
Ruled Surface(380) = {380};\n\
Line Loop(382) = {173, -215, -165, 64};\n\
Ruled Surface(382) = {382};\n\
Line Loop(384) = {215, 196, 84, -166};\n\
Ruled Surface(384) = {384};\n\
Line Loop(386) = {165, -214, -181, 63};\n\
Ruled Surface(386) = {386};\n\
Line Loop(388) = {214, 166, 85, -188};\n\
Ruled Surface(388) = {388};\n\
Line Loop(390) = {181, -187, 148, 62};\n\
Ruled Surface(390) = {390};\n\
Line Loop(392) = {187, -213, 80, 147};\n\
Ruled Surface(392) = {392};\n\
Line Loop(394) = {213, 182, 160, 79};\n\
Ruled Surface(394) = {394};\n\
Line Loop(396) = {182, -159, -86, -188};\n\
Ruled Surface(396) = {396};\n\
Line Loop(398) = {160, -78, -208, -183};\n\
Ruled Surface(398) = {398};\n\
Line Loop(400) = {183, -197, -87, 159};\n\
Ruled Surface(400) = {400};\n\
Line Loop(402) = {208, -77, -146, -198};\n\
Ruled Surface(402) = {402};\n\
Line Loop(404) = {145, -198, 184, 56};\n\
Ruled Surface(404) = {404};\n\
Line Loop(406) = {184, -55, -168, 207};\n\
Ruled Surface(406) = {406};\n\
Line Loop(408) = {207, -197, 88, 167};\n\
Ruled Surface(408) = {408};\n\
Line Loop(410) = {168, -54, -176, 206};\n\
Ruled Surface(410) = {410};\n\
Line Loop(412) = {206, -167, 81, 189};\n\
Ruled Surface(412) = {412};\n\
Line Loop(414) = {176, -53, -152, -190};\n\
Ruled Surface(414) = {414};\n\
Line Loop(416) = {151, -190, -205, 72};\n\
Ruled Surface(416) = {416};\n\
Line Loop(418) = {205, -175, -157, 71};\n\
Ruled Surface(418) = {418};\n\
Line Loop(420) = {175, -189, 82, -158};\n\
Ruled Surface(420) = {420};\n\
')
# Interior surface indices - from 230 to 420 (step 2), inclusive

# Add closing surfaces, fluid surfaces, 1006 - will be closed later
# Line Loop(1006) = {218, 219, 220, 217, -44, -43, -42, -41, -48, -47, -46, -45};\n\
fidF.write('\
Line Loop(1001) = {225, -224, -223, 219, -32, -31, -30, -29, -28, -27, -26, -25};\n\
Line Loop(1002) = {221, 227, -222, -217, -17, -24, -23, -22, -21, -20, -19, -18};\n\
Line Loop(1003) = {226, -221, -220, 225, -8, -7, -6, -5, -4, -3, -2, -1};\n\
Line Loop(1004) = {218, 223, -228, -222, -16, -15, -14, -13, -12, -11, -10, -9};\n\
Line Loop(1005) = {227, 228, 224, 226, -36, -35, -34, -33, -40, -39, -38, -37};\n\
')
# Define plane surfaces
for i in range(5):
    fidF.write('Plane Surface(%d) = {%d};\n'%(1001+i,1001+i))

# Add closing surfaces, solid surfaces, 11006 - will be closed later
# Line Loop(11006) = {46, 47, 48, 41, 42, 43, 44, 45};\n\
fidF.write('\
Line Loop(11001) = {32, 25, 26, 27, 28, 29, 30, 31};\n\
Line Loop(11002) = {23, 24, 17, 18, 19, 20, 21, 22};\n\
Line Loop(11003) = {1, 2, 3, 4, 5, 6, 7, 8};\n\
Line Loop(11004) = {9, 10, 11, 12, 13, 14, 15, 16};\n\
Line Loop(11005) = {38, 39, 40, 33, 34, 35, 36, 37};\n\
')
# Define plane surfaces
for i in range(5):
    fidF.write('Plane Surface(%d) = {%d};\n'%(11001+i,11001+i))


# Do extension of unit cell down to the bottom of the coating
for i in range(nCell-1):
    fidF.write('FreeSurface%d[] = Translate {0,0,-%d} {Duplicata{ Surface{'%(i+1,i+1) )
    for j in range(95):
        curInd = 230+j*2
        fidF.write('%d, '%(curInd))
    fidF.write('%d}; }};\n'%(curInd+2))
for i in range(nCell-1):
    for j in range(4):
        fidF.write('P%dSurface%d[]  = Translate {0,0,-%d} {Duplicata{ Surface{%d}; }};\n'%(1001+j,i+1,i+1, 1001+j))
        fidF.write('P%dSurface%dS[] = Translate {0,0,-%d} {Duplicata{ Surface{%d}; }};\n'%(1001+j,i+1,i+1,11001+j))

# Translate the closing surface
fidF.write('BSurface[]  = Translate {0,0,-%d} {Duplicata{ Surface{ 1005}; }};\n'%(nCell-1))
fidF.write('Delete { Surface{ 1005}; }\n')
fidF.write('BSurfaceS[] = Translate {0,0,-%d} {Duplicata{ Surface{11005}; }};\n'%(nCell-1))
fidF.write('Delete { Surface{11005}; }\n')


# Close the surface, where the cylinder was attached
fidF.write('Line Loop(10061) = {')
for j in range(7):
    fidF.write('%d, '%(41+j))
fidF.write('48};\n')
fidF.write('Plane Surface(10061) = {10061};\n')

# Now close the top fluid surface
fidF.write('Line Loop(1006) = {217, 218, 219, 220};\n')
fidF.write('Plane Surface(1006) = {1006};\n')

# And translate the top fluid surface
fidF.write('Translate {0,0,%.4f} { Surface{1006}; }\n'%(zmax))

# ------------------------------------------------------------------------------------------------
# Extending the domain to describe full cavity
# ------------------------------------------------------------------------------------------------
# Set the variable, which contains all the surfaces of the top solid structure
fidF.write('FreeSurfaceTOP = {')
for i in range(96):
    curInd = 230+i*2
    fidF.write('%d, '%(curInd))
fidF.write('10061};\n')

# Copying the solid structure in the stream-wise direction
for j in range(nCellX-1):
    fidF.write('FreeSurfaceTOP%d[] = Translate {%d,0,0} {Duplicata{ Surface{FreeSurfaceTOP[]}; }};\n'%(j+1,j+1) )
    for i in range(nCell-1):
        fidF.write('FreeSurface%d%d[] = Translate {%d,0,0} {Duplicata{ Surface{FreeSurface%d[]}; }};\n'%(i+1,j+1,j+1,i+1) )
# Copying the side, bottom and top walls
for j in range(nCellX-1):
    fidF.write('S1006%d[] = Translate {%d,0,0} {Duplicata{ Surface{1006}; }};\n'%(j+1,j+1))
    fidF.write('BSurface%d[]  = Translate {%d,0,0} {Duplicata{ Surface{BSurface []}; }};\n'%(j+1,j+1))
    fidF.write('BSurface%dS[] = Translate {%d,0,0} {Duplicata{ Surface{BSurfaceS[]}; }};\n'%(j+1,j+1))
    fidF.write('P1003SurfaceTOP%d[]  = Translate {%d,0,0} {Duplicata{ Surface{ 1003}; }};\n'%(j+1,j+1))
    fidF.write('P1004SurfaceTOP%d[]  = Translate {%d,0,0} {Duplicata{ Surface{ 1004}; }};\n'%(j+1,j+1))
    fidF.write('P1003SurfaceTOP%dS[] = Translate {%d,0,0} {Duplicata{ Surface{11003}; }};\n'%(j+1,j+1))
    fidF.write('P1004SurfaceTOP%dS[] = Translate {%d,0,0} {Duplicata{ Surface{11004}; }};\n'%(j+1,j+1))
    for i in range(nCell-1):    
        fidF.write('P%dSurface%d%d[]  = Translate {%d,0,0} {Duplicata{ Surface{P%dSurface%d[]};  }};\n'%(1003,i+1,j+1,j+1,1003,i+1))
        fidF.write('P%dSurface%d%d[]  = Translate {%d,0,0} {Duplicata{ Surface{P%dSurface%d[]};  }};\n'%(1004,i+1,j+1,j+1,1004,i+1))
        fidF.write('P%dSurface%d%dS[] = Translate {%d,0,0} {Duplicata{ Surface{P%dSurface%dS[]}; }};\n'%(1003,i+1,j+1,j+1,1003,i+1))
        fidF.write('P%dSurface%d%dS[] = Translate {%d,0,0} {Duplicata{ Surface{P%dSurface%dS[]}; }};\n'%(1004,i+1,j+1,j+1,1004,i+1))
# Copy and move the closing surfaces
fidF.write('S1002n[]  = Translate {%d,0,0} {Duplicata{ Surface{ 1002}; }};\n'%(nCellX-1) )
fidF.write('S1002nS[] = Translate {%d,0,0} {Duplicata{ Surface{11002}; }};\n'%(nCellX-1) )
for i in range(nCell-1):
    fidF.write('P%dnSurface%d[]  = Translate {%d,0,0} {Duplicata{ Surface{P%dSurface%d[]};  }};\n'%(1002,i+1,nCellX-1,1002,i+1) )
    fidF.write('P%dnSurface%dS[] = Translate {%d,0,0} {Duplicata{ Surface{P%dSurface%dS[]}; }};\n'%(1002,i+1,nCellX-1,1002,i+1) )
# Delete the remaining interior surfaces
fidF.write('Delete { Surface{ 1002}; }\n')
fidF.write('Delete { Surface{11002}; }\n')
for i in range(nCell-1):
    fidF.write('Delete { Surface{P%dSurface%d[]};  }\n'%(1002,i+1))
    fidF.write('Delete { Surface{P%dSurface%dS[]}; }\n'%(1002,i+1))


# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------

# Define the volume, add the translated surface entities
# Fluid volume
fidF.write('Surface Loop(1301) = {')
fidF.write('%d, '%(1001))
fidF.write('%d, '%(1003))
fidF.write('%d, '%(1004))
fidF.write('%d, '%(1006))
fidF.write('%d, '%(10061))
for i in range(96):
    curInd = 230+i*2
    fidF.write('%d, '%curInd)
for i in range(nCell-1):
    fidF.write('P%dSurface%d[], '%(1001,i+1))
    fidF.write('P%dSurface%d[], '%(1003,i+1))
    fidF.write('P%dSurface%d[], '%(1004,i+1))
for i in range(nCell-1):
    fidF.write('FreeSurface%d[], '%(i+1))
fidF.write('BSurface[], ')
for j in range(nCellX-1):
    fidF.write('S1006%d[], '%(j+1))
    fidF.write('P1003SurfaceTOP%d[], '%(j+1))
    fidF.write('P1004SurfaceTOP%d[], '%(j+1))
    fidF.write('BSurface%d[], '%(j+1))
    fidF.write('FreeSurfaceTOP%d[], '%(j+1) )
    for i in range(nCell-1):
        fidF.write('FreeSurface%d%d[], '%(i+1,j+1) )
        fidF.write('P%dSurface%d%d[], '%(1003,i+1,j+1))
        fidF.write('P%dSurface%d%d[], '%(1004,i+1,j+1))
fidF.write('S1002n[], ')
for i in range(nCell-2):
    fidF.write('P%dnSurface%d[], '%(1002,i+1) )
fidF.write('P%dnSurface%d[]};\n'%(1002,i+2) )
fidF.write('Volume(1401) = {1301};\n')

# Solid volume
fidF.write('Surface Loop(1302) = {')
fidF.write('%d, '%(11001))
fidF.write('%d, '%(11003))
fidF.write('%d, '%(11004))
fidF.write('%d, '%(10061))
for i in range(96):
    curInd = 230+i*2
    fidF.write('%d, '%curInd)
for i in range(nCell-1):
    fidF.write('P%dSurface%dS[], '%(1001,i+1))
    fidF.write('P%dSurface%dS[], '%(1003,i+1))
    fidF.write('P%dSurface%dS[], '%(1004,i+1))
for i in range(nCell-1):
    fidF.write('FreeSurface%d[], '%(i+1))
fidF.write('BSurfaceS[], ')
for j in range(nCellX-1):
    fidF.write('P1003SurfaceTOP%dS[], '%(j+1))
    fidF.write('P1004SurfaceTOP%dS[], '%(j+1))
    fidF.write('BSurface%dS[], '%(j+1))
    fidF.write('FreeSurfaceTOP%d[], '%(j+1) )
    for i in range(nCell-1):
        fidF.write('FreeSurface%d%d[], '%(i+1,j+1) )
        fidF.write('P%dSurface%d%dS[], '%(1003,i+1,j+1))
        fidF.write('P%dSurface%d%dS[], '%(1004,i+1,j+1))
fidF.write('S1002nS[], ')
for i in range(nCell-2):
    fidF.write('P%dnSurface%dS[], '%(1002,i+1) )
fidF.write('P%dnSurface%dS[]};\n'%(1002,i+2) )
fidF.write('Volume(1402) = {1302};\n')

# Define the periodic surfaces, direction vector to translate the right surface for matching
# In this case, left is top surface, right is bottom surface (must be raised by one unit)
fidF.write('Periodic Surface {S1002n[]}  = { 1001} Translate {%d,0,0};\n'%nCellX);
fidF.write('Periodic Surface {S1002nS[]} = {11001} Translate {%d,0,0};\n'%nCellX);
fidF.write('Periodic Surface { 1004}     = { 1003} Translate {0 ,1,0};\n');
fidF.write('Periodic Surface {11004}     = {11003} Translate {0 ,1,0};\n');
# Translated periodic surfaces
for i in range(nCell-1):
    fidF.write('Periodic Surface {P1002nSurface%d[]}  = {P1001Surface%d[]}  Translate {%d,0,0};\n'%(i+1,i+1,nCellX));
    fidF.write('Periodic Surface {P1004Surface%d[]}   = {P1003Surface%d[]}  Translate {0,1,0};\n'%(i+1,i+1));    
    fidF.write('Periodic Surface {P1002nSurface%dS[]} = {P1001Surface%dS[]} Translate {%d,0,0};\n'%(i+1,i+1,nCellX));
    fidF.write('Periodic Surface {P1004Surface%dS[]}  = {P1003Surface%dS[]} Translate {0,1,0};\n'%(i+1,i+1));    
# Periodic surfaces in the remaining direction
for j in range(nCellX-1):
    fidF.write('Periodic Surface {P1004SurfaceTOP%d[]}  = {P1003SurfaceTOP%d[]}  Translate {0,1,0};\n'%(j+1,j+1))
    fidF.write('Periodic Surface {P1004SurfaceTOP%dS[]} = {P1003SurfaceTOP%dS[]} Translate {0,1,0};\n'%(j+1,j+1))
    for i in range(nCell-1):
        fidF.write('Periodic Surface {P1004Surface%d%d[]}  = {P1003Surface%d%d[]}  Translate {0,1,0};\n'%(i+1,j+1,i+1,j+1));
        fidF.write('Periodic Surface {P1004Surface%d%dS[]} = {P1003Surface%d%dS[]} Translate {0,1,0};\n'%(i+1,j+1,i+1,j+1));

# Make sure that we output only entities denoted using physical labels
fidF.write('Mesh.SaveAll = 0;\n')
# Set the mesh resolution
fidF.write('Mesh.CharacteristicLengthMax = %.4f;\n'%(resM2))

fidF.close()

# Copy the fluid file to the solid file
sh.copyfile(GeoFileF,GeoFileS)

# Append physical definitions
fidF = open(GeoFileF,'a')
# Define the physical surfaces, FLUID
# Periodic boundaries, including the translated ones (need to have separate
# pysical labels, because they are disconnected)
# LEFT BOUNDARY
fidF.write('Physical Surface(3001) = {1001, ')
for i in range(nCell-2):
    fidF.write('P1001Surface%d[], '%(i+1))
fidF.write('P1001Surface%d[]};\n'%(i+2))
# RIGHT BOUNDARY
fidF.write('Physical Surface(3002) = {S1002n[], ')
for i in range(nCell-2):
    fidF.write('P1002nSurface%d[], '%(i+1))
fidF.write('P1002nSurface%d[]};\n'%(i+2))
# FRONT BOUNDARY
fidF.write('Physical Surface(3003) = {1003, ')
for i in range(nCell-2):
    fidF.write('P1003Surface%d[], '%(i+1))
for j in range(nCellX-1):
    fidF.write('P1003SurfaceTOP%d[], '%(j+1))    
    for i in range(nCell-1):
        fidF.write('P1003Surface%d%d[], '%(i+1,j+1))
fidF.write('P1003Surface%d[]};\n'%(nCell-1))
# BACK BOUNDARY
fidF.write('Physical Surface(3004) = {1004, ')
for i in range(nCell-2):
    fidF.write('P1004Surface%d[], '%(i+1))
for j in range(nCellX-1):
    fidF.write('P1004SurfaceTOP%d[],'%(j+1))    
    for i in range(nCell-1):
        fidF.write('P1004Surface%d%d[], '%(i+1,j+1))
fidF.write('P1004Surface%d[]};\n'%(nCell-1))
# BOTTOM BOUNDARY
fidF.write('Physical Surface(3005) = {BSurface[], ' )
for j in range(nCellX-2):
    fidF.write('BSurface%d[], '%(j+1))    
fidF.write('BSurface%d[]};\n'%(j+2))    
# TOP BOUNDARY
fidF.write('Physical Surface(3006) = {1006, ' )
for j in range(nCellX-2):
    fidF.write('S1006%d[], '%(j+1))    
fidF.write('S1006%d[]};\n'%(j+2))    
# SOLID BOUNDARY
fidF.write('Physical Surface(3007) = {')
for i in range(96):
    curInd = 230+i*2
    fidF.write('%d, '%(curInd))
for i in range(nCell-1):
    fidF.write('FreeSurface%d[], '%(i+1))
for j in range(nCellX-1):
    fidF.write('FreeSurfaceTOP%d[], '%(j+1) )    
    for i in range(nCell-1):
        fidF.write('FreeSurface%d%d[], '%(i+1,j+1) )        
fidF.write('10061};\n')
# Define the physical volume, FLUID
fidF.write('Physical Volume(4001) = {1401};\n\n' )
# Close the file
fidF.close()


# Append physical definitions
fidS = open(GeoFileS,'a')
# Define the physical surfaces, SOLID
# Periodic boundaries, including the translated ones (need to have separate
# pysical labels, because they are disconnected)
# LEFT BOUNDARY
fidS.write('Physical Surface(3001) = {11001, ')
for i in range(nCell-2):
    fidS.write('P1001Surface%dS[], '%(i+1))
fidS.write('P1001Surface%dS[]};\n'%(i+2))
# RIGHT BOUNDARY
fidS.write('Physical Surface(3002) = {S1002nS[], ')
for i in range(nCell-2):
    fidS.write('P1002nSurface%dS[], '%(i+1))
fidS.write('P1002nSurface%dS[]};\n'%(i+2))
# FRONT BOUNDARY
fidS.write('Physical Surface(3003) = {11003, ')
for i in range(nCell-2):
    fidS.write('P1003Surface%dS[], '%(i+1))
for j in range(nCellX-1):
    fidS.write('P1003SurfaceTOP%dS[], '%(j+1))    
    for i in range(nCell-1):
        fidS.write('P1003Surface%d%dS[], '%(i+1,j+1))
fidS.write('P1003Surface%dS[]};\n'%(nCell-1))
# BACK BOUNDARY
fidS.write('Physical Surface(3004) = {11004, ')
for i in range(nCell-2):
    fidS.write('P1004Surface%dS[], '%(i+1))
for j in range(nCellX-1):
    fidS.write('P1004SurfaceTOP%dS[],'%(j+1))    
    for i in range(nCell-1):
        fidS.write('P1004Surface%d%dS[], '%(i+1,j+1))
fidS.write('P1004Surface%dS[]};\n'%(nCell-1))
# BOTTOM BOUNDARY
fidS.write('Physical Surface(3005) = {BSurfaceS[], ' )
for j in range(nCellX-2):
    fidS.write('BSurface%dS[], '%(j+1))    
fidS.write('BSurface%dS[]};\n'%(j+2))    
# SOLID BOUNDARY
fidS.write('Physical Surface(3007) = {')
for i in range(96):
    curInd = 230+i*2
    fidS.write('%d, '%(curInd))
for i in range(nCell-1):
    fidS.write('FreeSurface%d[], '%(i+1))
for j in range(nCellX-1):
    fidS.write('FreeSurfaceTOP%d[], '%(j+1) )    
    for i in range(nCell-1):
        fidS.write('FreeSurface%d%d[], '%(i+1,j+1) )        
fidS.write('10061};\n')
# Define the physical volume, SOLID
fidS.write('Physical Volume(4001) = {1402};\n\n' )
# Close the file
fidS.close()
