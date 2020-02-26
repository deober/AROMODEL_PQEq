import Atom
import Molecule
import OPLS
import System
import sys
import DA_Polymer
import math

#Modify and run this script to generate a lammps data file and control file
#Molecule class accepts .xyz or Lammps single molecule data files
#System combines molecule data files into a single data file describing the full system of solvents and solutes


def main():

    #thisScript is junk, just necessary for unpacking the system argument vector (sys.argv)
    #outputFilenames is a general simulation name that will be appended as in.<outputFilename> and dat.<outputFilename>

    soluteFilename = '/home/derick/SimSandbox/Projects/PQEq/PDTSTPD/PDTSTPD_dimer.data'
    nSoluteMolecules = 1
    solventFilename = '/home/derick/SimSandbox/Projects/PQEq/PDTSTPD/noSolvent.dat'
    nSolventMolecules = 0   

    systemPrefix = 'PDTSTPD-PQEq-test'
    systemName = '%s_%s%d_%s%d' % (systemPrefix,soluteFilename.split('/')[-1].split('.')[0],int(nSoluteMolecules),solventFilename.split('/')[-1].split('.')[0],int(nSolventMolecules)) 


    nSoluteMolecules = int(nSoluteMolecules)
    nSolventMolecules = int(nSolventMolecules)
    dataName = '%s.data' % (systemName)
    print(systemName)
    inputName = 'in.%s' % (systemName)
    print('file names will be %s and %s' % (inputName,dataName))


    solute = Molecule.Molecule(soluteFilename)
    solute.UnConverged = True
    solute.Set_Up_FF(run_orca=True,local=True)

    fullSystem = System.System([solute],[nSoluteMolecules],Box_Size=20,Name=systemPrefix)
    fullSystem.Assign_PQEq(atomParameters=fullSystem.Atom_Params)
    fullSystem.Write_LAMMPS_PQEq(inputName,dataName,systemPrefix,300,fullSystem,)

if __name__=='__main__': main()
