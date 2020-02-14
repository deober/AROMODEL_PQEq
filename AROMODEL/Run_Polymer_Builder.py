import Molecule
import sys
import OPLS
import System

Script, Name = sys.argv

Polymer = Molecule.Molecule('%s.xyz' % Name)
Polymer.Set_Up_FF(local = False)
OPLS.Assign_OPLS(Polymer, ChelpG = False)
Initial_System = System.System([Polymer], [1], 100.0, Name)
Initial_System.Gen_Rand()
Initial_System.Write_LAMMPS_Data()
