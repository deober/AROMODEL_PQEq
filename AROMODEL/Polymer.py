#! usr/bin/python

# Header file containing functions specific to D14

Mass = { 1:32.060, 2:12.011, 3:12.011, 4:12.011, 5:12.011, 6:12.011, 7:1.008, 8:15.99, 9:12.011, 10: 1.008, 11:12.011, 12:1.008, 13:12.011, 14:12.011, 15:32.060, 16:12.011, 17:12.011, 18:32.060, 19:12.011, 20:12.011, 21:12.011, 22:12.011, 23:18.998, 24:15.999, 25:15.999, 26:12.011, 27:12.011, 28:1.008, 29:12.011, 30:1.008, 31:1.008 }

import Analyze_Traj

Num_CG = 7920 # Number of CG Beads in simulation
DP = 12 # Degree of polymerization
Atoms_Per_Monomer = 103
Donor_Index = 0
Acceptor_Index = 7
Num_Beads = 11 # Number of CG_Beads per monomer
Bond_Map = [[1,4,7,-4],[-1,1],[-1,1],[-1],[-4,1],[-1,1],[-1],[-7,1,4],[-1,1],[-1,1],[-1]] #relative positions of bonded atoms
Bond_Groups = {}

def Map2CG(SS):
    k=1
    Mol_Num = 1
    for Mol in SS.Mol_List:
        for i in range(0,DP):
            if i == 0 or i == DP:
                j = i*Atoms_Per_Monomer
                SS.CG_List.append(Analyze_Traj.CG_Bead(k, "D", 1, Mol.Atom_List[j:j+13],Mol_Num))
                SS.CG_List[-1].End_Group = True
                SS.CG_List[-1].Assign_Basis([1,7], [3,9])
                j += 14
                SS.CG_List.append(Analyze_Traj.CG_Bead(k+1, "SO", 2, Mol.Atom_List[j:j+3],Mol_Num))
                SS.CG_List[-1].End_Group = True
                j += 4
                SS.CG_List.append(Analyze_Traj.CG_Bead(k+2, "SA1", 3, Mol.Atom_List[j:j+8],Mol_Num))
                SS.CG_List[-1].End_Group = True
                j += 9
                SS.CG_List.append(Analyze_Traj.CG_Bead(k+3, "SA2", 4, Mol.Atom_List[j:j+12],Mol_Num))
                SS.CG_List[-1].End_Group = True
                j += 13
                SS.CG_List.append(Analyze_Traj.CG_Bead(k+4, "SO", 2, Mol.Atom_List[j:j+3],Mol_Num))
                SS.CG_List[-1].End_Group = True
                j += 4
                SS.CG_List.append(Analyze_Traj.CG_Bead(k+5, "SA1", 3, Mol.Atom_List[j:j+8],Mol_Num))
                SS.CG_List[-1].End_Group = True
                j += 9
                SS.CG_List.append(Analyze_Traj.CG_Bead(k+6, "SA2", 4, Mol.Atom_List[j:j+12],Mol_Num))
                SS.CG_List[-1].End_Group = True
                j += 13
                SS.CG_List.append(Analyze_Traj.CG_Bead(k+7, "A", 5, Mol.Atom_List[j:j+9],Mol_Num))
                SS.CG_List[-1].End_Group = True
                SS.CG_List[-1].Assign_Basis([1,7], [0,4])
                j += 10
                SS.CG_List.append(Analyze_Traj.CG_Bead(k+8, "SCO", 7, Mol.Atom_List[j:j+2],Mol_Num))
                SS.CG_List[-1].End_Group = True
                j += 3
                SS.CG_List.append(Analyze_Traj.CG_Bead(k+9, "SA3", 8, Mol.Atom_List[j:j+11],Mol_Num))
                SS.CG_List[-1].End_Group = True
                j += 12
                SS.CG_List.append(Analyze_Traj.CG_Bead(k+10, "SA4", 9, Mol.Atom_List[j:j+12],Mol_Num))
                SS.CG_List[-1].End_Group = True
            else:
                j = i*Atoms_Per_Monomer
                SS.CG_List.append(Analyze_Traj.CG_Bead(k, "D", 1, Mol.Atom_List[j:j+13],Mol_Num))
                SS.CG_List[-1].Assign_Basis([1,7], [3,9])
                j += 14
                SS.CG_List.append(Analyze_Traj.CG_Bead(k+1, "SO", 2, Mol.Atom_List[j:j+3],Mol_Num))
                j += 4
                SS.CG_List.append(Analyze_Traj.CG_Bead(k+2, "SA1", 3, Mol.Atom_List[j:j+8],Mol_Num))
                j += 9
                SS.CG_List.append(Analyze_Traj.CG_Bead(k+3, "SA2", 4, Mol.Atom_List[j:j+12],Mol_Num))
                j += 13
                SS.CG_List.append(Analyze_Traj.CG_Bead(k+4, "SO", 2, Mol.Atom_List[j:j+3],Mol_Num))
                j += 4
                SS.CG_List.append(Analyze_Traj.CG_Bead(k+5, "SA1", 3, Mol.Atom_List[j:j+8],Mol_Num))
                j += 9
                SS.CG_List.append(Analyze_Traj.CG_Bead(k+6, "SA2", 4, Mol.Atom_List[j:j+12],Mol_Num))
                j += 13
                SS.CG_List.append(Analyze_Traj.CG_Bead(k+7, "A", 5, Mol.Atom_List[j:j+9],Mol_Num))
                SS.CG_List[-1].Assign_Basis([1,7], [0,4])
                j += 10
                SS.CG_List.append(Analyze_Traj.CG_Bead(k+8, "SCO", 7, Mol.Atom_List[j:j+2],Mol_Num))
                j += 3
                SS.CG_List.append(Analyze_Traj.CG_Bead(k+9, "SA3", 8, Mol.Atom_List[j:j+11],Mol_Num))
                j += 12
                SS.CG_List.append(Analyze_Traj.CG_Bead(k+10, "SA4", 9, Mol.Atom_List[j:j+12],Mol_Num))
            k += Num_Beads
        Mol_Num += 1
    for Bead in SS.CG_List:
         for bond in Bond_Map[(Bead.CG_ID - 1) % Num_Beads ]:
             if (Bead.CG_ID+bond) >= 1 and (Bead.CG_ID+bond) < len(SS.CG_List) and Bead.Mol_Num == SS.Find_Bead(Bead.CG_ID+bond).Mol_Num:
                Bead.Bond_List.append(Bead.CG_ID+bond)
                if Bead.CG_ID < Bead.CG_ID+bond:
                    try:
                        SS.CG_Bond_List.append(Analyze_Traj.CG_Bond((Bond_Groups[Bead.CG_Type + SS.Find_Bead(Bead.CG_ID+bond).CG_Type]),Bead,SS.Find_Bead(Bead.CG_ID+bond)))
                    except:
                        Bond_Groups[Bead.CG_Type + SS.Find_Bead(Bead.CG_ID+bond).CG_Type] = len(Bond_Groups)/2+1
                        Bond_Groups[SS.Find_Bead(Bead.CG_ID+bond).CG_Type + Bead.CG_Type] = len(Bond_Groups)/2+1
                        SS.CG_Bond_List.append(Analyze_Traj.CG_Bond((Bond_Groups[Bead.CG_Type + SS.Find_Bead(Bead.CG_ID+bond).CG_Type]),Bead,SS.Find_Bead(Bead.CG_ID+bond)))
    SS.Bond_Groups = Bond_Groups
#print SS.CG_Bond_List
        # Re-wrap coordinates corresponding to periodic boundary conditions
    for CG_Mol_Obj in SS.CG_Mol_List:
        for CG_Bead_Obj in CG_Mol_Obj.CG_Bead_List:
            for i in range(3):
                if CG_Bead_Obj.COM[i] > SS.Box_Dim[i]:
                    CG_Bead_Obj.COM_PBC[i] = CG_Bead_Obj.COM[i] - SS.Box_Dim[i]
                    CG_Bead_Obj.IF[i] += 1
                elif CG_Bead_Obj.COM[i] < 0.0:
                    CG_Bead_Obj.COM_PBC[i] = CG_Bead_Obj.COM[i] + SS.Box_Dim[i]
                    CG_Bead_Obj.IF[i] -= 1
                else:
                    CG_Bead_Obj.COM_PBC[i] = CG_Bead_Obj.COM[i]
                    CG_Bead_Obj.IF[i] = 0


    return
