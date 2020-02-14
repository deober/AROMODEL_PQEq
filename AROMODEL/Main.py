import Atom
import Molecule
import OPLS
import System
import sys
import DA_Polymer
import math

def main():
    """Script, File_Name, N = sys.argv
    N = int(N)
    print File_Name
    Name = File_Name.split('.')[0] + "_%s" % N"""
    Solvent = Molecule.Molecule(File_Name)
    Solvent.Scan_Dihedrals([4,5,10,11],local = True)
    Solvent.UnConverged = True
    Solvent.Set_Up_FF(run_orca=True, local = False)
    OPLS.Assign_OPLS(Solvent, ChelpG = False)
    Solvent_System = System.System([Solvent], [N], 200.0, Name)
    Solvent_System.Gen_Rand()
    Solvent_System.Write_LAMMPS_Data()
    """Solvent_System.Run_Lammps_Init()
    Solvent_System.Run_Lammps_NPT""
    

if __name__=='__main__': main()


