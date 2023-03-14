#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------
# Created By  : Sebastian Widmann
# Institution : TU Munich, Department of Aerospace and Geodesy
# Created Date: September 28, 2022
# version ='1.0'
# ---------------------------------------------------------------------------
"""
Class implementation to generate mesh using cfMesh. meshDict is generated by
specifying aerofoil and angle of attack and dictionary will be written into
the "/system" directory..
"""
# ---------------------------------------------------------------------------

import argparse
import subprocess
import numpy as np
import generateStlFile
from flowProperties import *


class generateMeshDict(object):
    def __init__(self, maxCellSize):
        self.airfoil = args.airfoil
        self.angle = args.alpha
        self.mach = args.mach
        self.maxCellSize = maxCellSize

        self.writeToFile()

    def calculateFirstLayerThickness(self):
        """
        Returns
        -------
        Minimum cell height in boundary layer
        """
        yplus = 30.0

        p = calculateStaticPressure(self.mach)
        T = calculateStaticTemperature(self.mach)
        rho = calculateStaticDensity(p, T)
        u = self.mach * calculateSpeedofSound(T)
        Re = rho * u / mu  # [-] Freestream Reynolds number
        cf = np.power(2 * np.log10(Re) - 0.65, -2.3)  # [-] Skin friction coefficient based on Schlichting
        Tau_w = 0.5 * cf * rho * u ** 2  # [Pa] Wall shear stress
        u_star = np.sqrt(Tau_w / rho)  # [m*s^-1] Friction velocity
        return yplus * mu / (rho * u_star)

    def writeToFile(self):
        f = open('system/meshDict', 'w+')

        f.write('/*--------------------------------*- C++ -*----------------------------------*\\   \n')
        f.write('| =========                 |                                                 |    \n')
        f.write('| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |  \n')
        f.write('|  \\\\    /   O peration     | Version:  v2206                                 |  \n')
        f.write('|   \\\\  /    A nd           | Website:  www.openfoam.com                      |  \n')
        f.write('|    \\\\/     M anipulation  |                                                 |  \n')
        f.write('\\*---------------------------------------------------------------------------*/   \n')
        f.write('FoamFile                                                                           \n')
        f.write('{                                                                                  \n')
        f.write('    version     2.0;                                                               \n')
        f.write('    format      ascii;                                                             \n')
        f.write('    class       dictionary;                                                        \n')
        f.write('    object      meshDict;                                                 \n')
        f.write('}                                                                                  \n')
        f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //    \n')
        f.write('surfaceFile    "constant/triSurface/naca{}_{}deg.stl"; \n'.format(self.airfoil, self.angle))
        f.write('                                                                                   \n')
        f.write('maxCellSize    {}; \n'.format(self.maxCellSize))
        f.write('                                                                                   \n')
        f.write('objectRefinements                                                                  \n')
        f.write('{                                                                                  \n')
        f.write('   refinementBox1                                                                  \n')
        f.write('   {                                                                               \n')
        f.write('       cellSize    {}; \n'.format(self.maxCellSize/(2**1)))
        f.write('       type        box;                                                            \n')
        f.write('       centre      (7.5 0 0);                                                      \n')
        f.write('       lengthX     25;                                                             \n')
        f.write('       lengthY     8;                                                              \n')
        f.write('       lengthZ     1;                                                              \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   refinementBox2                                                                  \n')
        f.write('   {                                                                               \n')
        f.write('       cellSize    {}; \n'.format(self.maxCellSize/(2**2)))
        f.write('       type        box;                                                            \n')
        f.write('       centre      (3.75 0 0);                                                     \n')
        f.write('       lengthX     12.5;                                                           \n')
        f.write('       lengthY     4;                                                              \n')
        f.write('       lengthZ     1;                                                              \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   refinementBox3                                                                  \n')
        f.write('   {                                                                               \n')
        f.write('       cellSize    {}; \n'.format(self.maxCellSize/(2**3)))
        f.write('       type        box;                                                            \n')
        f.write('       centre      (4 0 0);                                                        \n')
        f.write('       lengthX     4;                                                              \n')
        f.write('       lengthY     2;                                                              \n')
        f.write('       lengthZ     1;                                                              \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   refinementBox4                                                                  \n')
        f.write('   {                                                                               \n')
        f.write('       cellSize    {}; \n'.format(self.maxCellSize/(2**4)))
        f.write('       type        box;                                                            \n')
        f.write('       centre      (0.5 0 0);                                                      \n')
        f.write('       lengthX     2;                                                              \n')
        f.write('       lengthY     1;                                                              \n')
        f.write('       lengthZ     1;                                                              \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   refinementBox5                                                                  \n')
        f.write('   {                                                                               \n')
        f.write('       cellSize    {}; \n'.format(self.maxCellSize/(2**5)))
        f.write('       type        box;                                                            \n')
        f.write('       centre      (0.25 0 0);                                                     \n')
        f.write('       lengthX     1.2;                                                            \n')
        f.write('       lengthY     0.5;                                                            \n')
        f.write('       lengthZ     1;                                                              \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('}                                                                                  \n')
        f.write('                                                                                   \n')
        f.write('localRefinement                                                                    \n')
        f.write('{                                                                                  \n')
        f.write('   patch0                                                                          \n')
        f.write('   {                                                                               \n')
        f.write('       additionalRefinementLevels  6;                                              \n')
        # f.write('       refinementThickness         0.01;                                           \n')
        f.write('   }                                                                               \n')
        f.write('}                                                                                  \n')
        f.write('                                                                                   \n')
        f.write('boundaryLayers                                                                     \n')
        f.write('{                                                                                  \n')
        f.write('   patchBoundaryLayers                                                             \n')
        f.write('   {                                                                               \n')
        f.write('       patch0                                                                      \n')
        f.write('       {                                                                           \n')
        f.write('           nLayers                 10;                                             \n')
        f.write('           thicknessRatio          1.2;                                            \n')
        f.write('           maxFirstLayerThickness  {}; \n'.format(self.calculateFirstLayerThickness()))
        f.write('           allowDiscontinuity      1;                                              \n')
        f.write('                                                                                   \n')
        f.write('           optimiseLayer           1;                                              \n')
        f.write('       }                                                                           \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   optimisationParameters                                                          \n')
        f.write('   {                                                                               \n')
        f.write('       nSmoothNormals              10;                                             \n')
        f.write('       maxNumInterations           10;                                             \n')
        f.write('       featureSizeFactor           0.1;                                            \n')
        f.write('       relThicknessTol             0.1;                                            \n')
        f.write('   }                                                                               \n')
        f.write('}                                                                                  \n')
        f.write('                                                                                   \n')
        f.write('renameBoundary                                                                     \n')
        f.write('{                                                                                  \n')
        f.write('   newPatchNames                                                                   \n')
        f.write('   {                                                                               \n')
        f.write('       "xMin"                                                                      \n')
        f.write('       {                                                                           \n')
        f.write('           newName inlet;                                                          \n')
        f.write('           type    patch;                                                          \n')
        f.write('       }                                                                           \n')
        f.write('                                                                                   \n')
        f.write('       "xMax"                                                                      \n')
        f.write('       {                                                                           \n')
        f.write('           newName outlet;                                                         \n')
        f.write('           type    patch;                                                          \n')
        f.write('       }                                                                           \n')
        f.write('                                                                                   \n')
        f.write('       "yMin"                                                                      \n')
        f.write('       {                                                                           \n')
        f.write('           newName bottom;                                                         \n')
        f.write('           type    patch;                                                          \n')
        f.write('       }                                                                           \n')
        f.write('                                                                                   \n')
        f.write('       "yMax"                                                                      \n')
        f.write('       {                                                                           \n')
        f.write('           newName top;                                                            \n')
        f.write('           type    patch;                                                          \n')
        f.write('       }                                                                           \n')
        # f.write('                                                                                   \n')
        # f.write('       "zMin"                                                                      \n')
        # f.write('       {                                                                           \n')
        # f.write('           newName back;                                                           \n')
        # f.write('           type    empty;                                                          \n')
        # f.write('       }                                                                           \n')
        # f.write('                                                                                   \n')
        # f.write('       "zMax"                                                                      \n')
        # f.write('       {                                                                           \n')
        # f.write('           newName front;                                                          \n')
        # f.write('           type    empty;                                                          \n')
        # f.write('       }                                                                           \n')
        f.write('                                                                                   \n')
        f.write('       "patch0"                                                                    \n')
        f.write('       {                                                                           \n')
        f.write('           newName wing;                                                           \n')
        f.write('           type    wall;                                                           \n')
        f.write('       }                                                                           \n')
        f.write('   }                                                                               \n')
        f.write('}                                                                                  \n')
        f.write('                                                                                   \n')
        f.write('// ************************************************************************* //    \n')
        f.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate snappyHexMeshDict file and save into "system/blockMeshFile"')
    parser.add_argument('airfoil', type=str, help='NACA airfoil digits')
    parser.add_argument('alpha', type=float, help='Angle of attack [deg]')
    parser.add_argument('mach', type=float, help='Freestream Mach number [-]')
    parser.add_argument('maxCellSize', type=float, help='Base cell size of domain [m]')
    args = parser.parse_args()

    # Generate required NACA airfoil
    subprocess.call("python scripts/generateStlFile.py {} {}".format(args.airfoil, args.alpha), shell=True)

    # Generate meshDict into system/ directory
    generateMeshDict(args.maxCellSize)
