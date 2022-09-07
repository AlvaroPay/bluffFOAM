#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------
# Created By  : Sebastian Widmann
# Institution : TU Munich, Department of Aerospace and Geodesy
# Created Date: June 16, 2022
# version ='1.0'
# ---------------------------------------------------------------------------
"""
Class implementation to generate background mesh with blockMesh. Background
mesh is specified by base cell size and blockMeshDict will be written into
the "/system" directory.
"""
# ---------------------------------------------------------------------------

import argparse
import numpy as np
from scipy.optimize import newton

class generateBlockMeshDict(object):
    def __init__(self, baseCellSize):
        # Domain size
        self.xMin = self.xMax = None
        self.zMin = self.zMax = None
        self.yMin = self.yMax = None

        # Mesh parameter
        self.xCells = self.yCells = None
        self.baseCellSize = baseCellSize

        self.setDomainSize()
        self.writeToFile()

    def setDomainSize(self, xMin=-10, xMax=30, yMin=-10, yMax=10):
        # Domain size parameter
        self.xMin = xMin
        self.xMax = xMax
        self.yMin = yMin
        self.yMax = yMax
        self.zMin = -0.1
        self.zMax = 0.1

        self.xCells = self.setNumberOfCells(1, abs(self.xMax) + abs(self.xMin), self.baseCellSize)
        self.yCells = self.setNumberOfCells(1, abs(self.yMax) + abs(self.yMin), self.baseCellSize)

    def setNumberOfCells(self, ratio, length, deltaS):
        """
        Parameters
        ----------
        ratio: total expansion ratio
        length: length of block edge
        deltaS: width / height of start cell

        Returns
        -------
        Number of cells as integer value
        """

        if ratio > 1:
            cMin = deltaS
        else:
            cMin = deltaS * ratio

        if abs(ratio - 1) < 1e-5:
            return round(length / cMin)
        else:
            return round(newton(
                lambda n: (1 - np.power(ratio, (n / (n - 1)))) / (1 - np.power(ratio, (1 / (n - 1)))) - length / deltaS,
                length / cMin))

    def writeToFile(self):
        f = open('system/blockMeshDict', 'w+')

        f.write('/*--------------------------------*- C++ -*----------------------------------*\\   \n')
        f.write('| =========                 |                                                 |    \n')
        f.write('| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |  \n')
        f.write('|  \\\\    /   O peration     | Version:  v2112                                 |  \n')
        f.write('|   \\\\  /    A nd           | Website:  www.openfoam.com                      |  \n')
        f.write('|    \\\\/     M anipulation  |                                                 |  \n')
        f.write('\\*---------------------------------------------------------------------------*/   \n')
        f.write('FoamFile                                                                           \n')
        f.write('{                                                                                  \n')
        f.write('    version     2.0;                                                               \n')
        f.write('    format      ascii;                                                             \n')
        f.write('    class       dictionary;                                                        \n')
        f.write('    object      blockMeshDict;                                                     \n')
        f.write('}                                                                                  \n')
        f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //    \n')
        f.write('                                                                                   \n')
        f.write('domain                                                                             \n')
        f.write('{                                                                                  \n')
        f.write('   xMin        {}; \n'.format(self.xMin))
        f.write('   xMax        {}; \n'.format(self.xMax))
        f.write('   yMin        {}; \n'.format(self.yMin))
        f.write('   yMax        {}; \n'.format(self.yMax))
        f.write('   zMin        {}; \n'.format(self.zMin))
        f.write('   zMax        {}; \n'.format(self.zMax))
        f.write('}                                                                                  \n')
        f.write('                                                                                   \n')
        f.write('vertices                                                                           \n')
        f.write('(                                                                                  \n')
        f.write('   ($domain.xMin   $domain.yMin $domain.zMin)                                      \n')
        f.write('   ($domain.xMax   $domain.yMin $domain.zMin)                                      \n')
        f.write('   ($domain.xMax   $domain.yMax $domain.zMin)                                      \n')
        f.write('   ($domain.xMin   $domain.yMax $domain.zMin)                                      \n')
        f.write('                                                                                   \n')
        f.write('   ($domain.xMin   $domain.yMin $domain.zMax)                                      \n')
        f.write('   ($domain.xMax   $domain.yMin $domain.zMax)                                      \n')
        f.write('   ($domain.xMax   $domain.yMax $domain.zMax)                                      \n')
        f.write('   ($domain.xMin   $domain.yMax $domain.zMax)                                      \n')
        f.write(');                                                                                 \n')
        f.write('                                                                                   \n')
        f.write('blocks                                                                             \n')
        f.write('(                                                                                  \n')
        f.write('   hex (0 1 2 3 4 5 6 7)                                                           \n')
        f.write('   ({} {} 1) \n'.format(self.xCells, self.yCells))
        f.write('   simpleGrading (1 1 1)                                                           \n')
        f.write(');                                                                                 \n')
        f.write('                                                                                   \n')
        f.write('boundary                                                                           \n')
        f.write('(                                                                                  \n')
        f.write('   inlet                                                                           \n')
        f.write('   {                                                                               \n')
        f.write('       type patch;                                                                 \n')
        f.write('       faces                                                                       \n')
        f.write('       (                                                                           \n')
        f.write('           (0 4 7 3)                                                               \n')
        f.write('       );                                                                          \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   outlet                                                                          \n')
        f.write('   {                                                                               \n')
        f.write('       type patch;                                                                 \n')
        f.write('       faces                                                                       \n')
        f.write('       (                                                                           \n')
        f.write('           (1 2 6 5)                                                               \n')
        f.write('       );                                                                          \n')
        f.write('                                                                                   \n')
        f.write('   }                                                                               \n')
        f.write('   top                                                                             \n')
        f.write('                                                                                   \n')
        f.write('   {                                                                               \n')
        f.write('       type patch;                                                                 \n')
        f.write('       faces                                                                       \n')
        f.write('       (                                                                           \n')
        f.write('           (2 3 7 6)                                                               \n')
        f.write('       );                                                                          \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   bottom                                                                          \n')
        f.write('   {                                                                               \n')
        f.write('       type patch;                                                                 \n')
        f.write('       faces                                                                       \n')
        f.write('       (                                                                           \n')
        f.write('           (0 1 5 4)                                                               \n')
        f.write('       );                                                                          \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   back                                                                            \n')
        f.write('   {                                                                               \n')
        f.write('       type symmetryPlane;                                                         \n')
        f.write('       faces                                                                       \n')
        f.write('       (                                                                           \n')
        f.write('           (0 1 2 3)                                                               \n')
        f.write('       );                                                                          \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   front                                                                           \n')
        f.write('   {                                                                               \n')
        f.write('       type symmetryPlane;                                                         \n')
        f.write('       faces                                                                       \n')
        f.write('       (                                                                           \n')
        f.write('           (4 5 6 7)                                                               \n')
        f.write('       );                                                                          \n')
        f.write('   }                                                                               \n')
        f.write(');                                                                                 \n')
        f.write('                                                                                   \n')
        f.write('// ************************************************************************* //    \n')
        f.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate blockMeshDict File')
    parser.add_argument('cellSize', type=float, help='Background mesh cell size')
    args = parser.parse_args()

    blockMesh = generateBlockMeshDict(baseCellSize=args.cellSize)
