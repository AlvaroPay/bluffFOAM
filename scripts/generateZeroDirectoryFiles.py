# ----------------------------------------------------------------------------
# Created By  : Sebastian Widmann
# Institution : TU Munich, Department of Aerospace and Geodesy
# Created Date: June 18, 2022
# version ='1.0'
# ---------------------------------------------------------------------------
"""
Class implementation to generate the initial condition files in the zero
directory. Initial values generated for kinematic eddy viscosity nut,
initial pressure p, specific dissipation rate omega, initial temperature T,
turbulent kinetic energy k, turbulent thermal diffusivity alpha_t and initial
velocity u. Initial condition files written into the "/0.orig" directory.
"""
# ---------------------------------------------------------------------------

import numpy as np
import os
from math import ceil, floor
from scripts.flowProperties import *


class generateInitialConditionsFiles(object):
    def __init__(self, caseDir, mach):
        self.caseDir = caseDir
        self.mach = mach

        self.p = None  # [Pa] Static pressure
        self.T = None  # [K] Static temperature
        self.mDot = None # [kg * s^-1] Mass flow rate
        self.u = None  # [m * s^-1] Freestream velocity vector
        self.k_inf = None  # [m^2 * s^-2] Freestream turbulent kinetic energy
        self.omega_inf = None  # [s^-1] Freestream turbulence specific dissipation rate
        self.omega_wall = None  # [s^-1] Wall turbulence specific dissipation rate
        self.rho = 1.223
        
        # Constants
        self.I = 0.05  # [%] Turbulence intensity
        self.C_mu = 0.09  # [-] Constant related to k-omega SST model
        self.L = 1  # [m] Reference length

        self.calculatePressure()
        self.calculateTemperature()
        self.calculateVelocity()
        self.calculateMassFlowRate()
        self.calculateTurbulentKineticEnergy()
        self.calculateSpecificDissipationRate()
        self.writeToFile()

    def writeToFile(self):
        self.writeToFile_KinematicEddyViscosity()
        self.writeToFile_Pressure()
        self.writeToFile_SpecificDissipationRate()
        self.writeToFile_TurbulentKineticEnergy()
        self.writeToFile_Velocity()

    def calculatePressure(self):
        self.p = calculateStaticPressure(self.mach)

    def calculateTemperature(self):
        self.T = calculateStaticTemperature(self.mach)

    def calculateVelocity(self):
        u = self.mach * calculateSpeedofSound(self.T)  # [m*s^-1] Freestream velocity
        self.u = np.array([u, 0, 0])

    def calculateMassFlowRate(self):
        self.mDot = calculateStaticDensity(self.p, self.T) * np.linalg.norm(self.u) * (25 * 1.0)

    def calculateTurbulentKineticEnergy(self):
        kLower = 1e-5 * np.linalg.norm(self.u)**2 / calculateReynoldsNumber(calculateStaticDensity(self.p, self.T), np.linalg.norm(self.u), 40, mu)
        kUpper = 1e-1 * np.linalg.norm(self.u)**2 / calculateReynoldsNumber(calculateStaticDensity(self.p, self.T), np.linalg.norm(self.u), 40, mu)

        self.k_inf = kLower

    def calculateSpecificDissipationRate(self):
        # TODO: implement inheritance for minimum cell layer thickness and domain size (replace 40 by windtunnel length)
        omegaLower = ceil(np.linalg.norm(self.u) / 40)
        omegaUpper = floor(10 * np.linalg.norm(self.u) / 40)

        self.omega_inf = omegaLower
        self.omega_wall = ceil(60 * (mu / calculateStaticDensity(self.p, self.T)) / (0.075 * calculateFirstLayerThickness(self.mach, 30)**2))

    def writeToFile_KinematicEddyViscosity(self):
        saveDir = os.path.join(self.caseDir, '0.org/nut')
        f = open(saveDir, 'w+')

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
        f.write('    class       volScalarField;                                                    \n')
        f.write('    object      nut;                                                               \n')
        f.write('}                                                                                  \n')
        f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //    \n')
        f.write('                                                                                   \n')
        f.write('dimensions             [0 2 -1 0 0 0 0];                                           \n')
        f.write('                                                                                   \n')
        f.write('internalField          uniform 0;                                                  \n')
        f.write('                                                                                   \n')
        f.write('boundaryField                                                                      \n')
        f.write('{                                                                                  \n')
        f.write('   inlet                                                                           \n')
        f.write('   {                                                                               \n')
        f.write('       type            calculated;                                                 \n')
        f.write('       value           $internalField;                                             \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   outlet                                                                          \n')
        f.write('   {                                                                               \n')
        f.write('       type            zeroGradient;                                               \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   wall                                                                            \n')
        f.write('   {                                                                               \n')
        f.write('       type            nutkWallFunction;                                           \n')
        f.write('       value           $internalField;                                             \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   "(top|bottom)"                                                                  \n')
        f.write('   {                                                                               \n')
        f.write('       type            slip;                                                       \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   #includeEtc "caseDicts/setConstraintTypes"                                      \n')
        f.write('}                                                                                  \n')
        f.write('                                                                                   \n')
        f.write('// ************************************************************************* //    \n')
        f.close()

    def writeToFile_Pressure(self):
        saveDir = os.path.join(self.caseDir, '0.org/p')
        f = open(saveDir, 'w+')

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
        f.write('    class       volScalarField;                                                    \n')
        f.write('    object      p;                                                                 \n')
        f.write('}                                                                                  \n')
        f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //    \n')
        f.write('                                                                                   \n')
        f.write('dimensions             [0 2 -2 0 0 0 0];                                           \n')
        f.write('                                                                                   \n')
        f.write('internalField          uniform {}; \n'.format(self.p))
        f.write('                                                                                   \n')
        f.write('boundaryField                                                                      \n')
        f.write('{                                                                                  \n')
        f.write('   inlet                                                                           \n')
        f.write('   {                                                                               \n')
        f.write('       type            zeroGradient;                                               \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   outlet                                                                          \n')
        f.write('   {                                                                               \n')
        f.write('       type            fixedValue;                                                 \n')
        f.write('       value           $internalField;                                             \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   wall                                                                            \n')
        f.write('   {                                                                               \n')
        f.write('       type            zeroGradient;                                               \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   "(top|bottom)"                                                                  \n')
        f.write('   {                                                                               \n')
        f.write('       type            slip;                                                       \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   #includeEtc "caseDicts/setConstraintTypes"                                      \n')
        f.write('}                                                                                  \n')
        f.write('                                                                                   \n')
        f.write('// ************************************************************************* //    \n')
        f.close()

    def writeToFile_SpecificDissipationRate(self):
        saveDir = os.path.join(self.caseDir, '0.org/omega')
        f = open(saveDir, 'w+')

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
        f.write('    class       volScalarField;                                                    \n')
        f.write('    object      omega;                                                             \n')
        f.write('}                                                                                  \n')
        f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //    \n')
        f.write('                                                                                   \n')
        f.write('dimensions             [0 0 -1 0 0 0 0];                                           \n')
        f.write('                                                                                   \n')
        f.write('internalField          uniform {}; \n'.format(self.omega_inf))
        f.write('                                                                                   \n')
        f.write('boundaryField                                                                      \n')
        f.write('{                                                                                  \n')
        f.write('   inlet                                                                           \n')
        f.write('   {                                                                               \n')
        f.write('       type            fixedValue;                                                 \n')
        f.write('       value           $internalField;                                             \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   outlet                                                                          \n')
        f.write('   {                                                                               \n')
        f.write('       type            zeroGradient;                                               \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   wall                                                                            \n')
        f.write('   {                                                                               \n')
        f.write('       type            omegaWallFunction;                                          \n')
        f.write('       value           uniform {}; \n'.format(self.omega_wall))
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   "(top|bottom)"                                                                  \n')
        f.write('   {                                                                               \n')
        f.write('       type            slip;                                                       \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   #includeEtc "caseDicts/setConstraintTypes"                                      \n')
        f.write('}                                                                                  \n')
        f.write('                                                                                   \n')
        f.write('// ************************************************************************* //    \n')
        f.close()

    def writeToFile_TurbulentKineticEnergy(self):
        saveDir = os.path.join(self.caseDir, '0.org/k')
        f = open(saveDir, 'w+')

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
        f.write('    class       volScalarField;                                                    \n')
        f.write('    object      k;                                                                 \n')
        f.write('}                                                                                  \n')
        f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //    \n')
        f.write('                                                                                   \n')
        f.write('dimensions             [0 2 -2 0 0 0 0];                                           \n')
        f.write('                                                                                   \n')
        f.write('internalField          uniform {}; \n'.format(self.k_inf))
        f.write('                                                                                   \n')
        f.write('boundaryField                                                                      \n')
        f.write('{                                                                                  \n')
        f.write('   inlet                                                                           \n')
        f.write('   {                                                                               \n')
        f.write('       type            fixedValue;                                                 \n')
        f.write('       value           $internalField;                                             \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   outlet                                                                          \n')
        f.write('   {                                                                               \n')
        f.write('       type            zeroGradient;                                               \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   wall                                                                            \n')
        f.write('   {                                                                               \n')
        f.write('       type            kqRWallFunction;                                            \n')
        f.write('       value           uniform 1e-20;                                              \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   "(top|bottom)"                                                                  \n')
        f.write('   {                                                                               \n')
        f.write('       type            slip;                                                       \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   #includeEtc "caseDicts/setConstraintTypes"                                      \n')
        f.write('}                                                                                  \n')
        f.write('                                                                                   \n')
        f.write('// ************************************************************************* //    \n')
        f.close()

    def writeToFile_Velocity(self):
        saveDir = os.path.join(self.caseDir, '0.org/U')
        f = open(saveDir, 'w+')

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
        f.write('    class       volVectorField;                                                    \n')
        f.write('    object      U;                                                                 \n')
        f.write('}                                                                                  \n')
        f.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //    \n')
        f.write('                                                                                   \n')
        f.write('dimensions             [0 1 -1 0 0 0 0];                                           \n')
        f.write('                                                                                   \n')
        f.write('internalField          uniform ({} {} {}); \n'.format(*self.u))
        f.write('                                                                                   \n')
        f.write('boundaryField                                                                      \n')
        f.write('{                                                                                  \n')
        f.write('   inlet                                                                           \n')
        f.write('   {                                                                               \n')
        f.write('       type            fixedValue;                                                 \n')
        f.write('       value           $internalField;                                             \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   outlet                                                                          \n')
        f.write('   {                                                                               \n')
        f.write('       type            inletOutlet;                                                \n')
        f.write('       inletValue      uniform (0 0 0);                                            \n')
        f.write('       value           $internalField;                                             \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   wall                                                                            \n')
        f.write('   {                                                                               \n')
        f.write('       type            noSlip;                                                     \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   "(top|bottom)"                                                                  \n')
        f.write('   {                                                                               \n')
        f.write('       type            slip;                                                       \n')
        f.write('   }                                                                               \n')
        f.write('                                                                                   \n')
        f.write('   #includeEtc "caseDicts/setConstraintTypes"                                      \n')
        f.write('}                                                                                  \n')
        f.write('                                                                                   \n')
        f.write('// ************************************************************************* //    \n')
        f.close()
