#! usr/bin/python
import Molecule
import copy
import numpy as np
import math
import OPLS
import Bond
import Angle
import Dihedral
import Improper

class Conjugated_Polymer(Molecule.Molecule):
    def __init__(self,Ring_List,Scheduler = "Torque",Cluster_Location='/oasis/tscc/scratch/andrewk/Optimized_Monomers'):
        self.Ring_List = copy.deepcopy(Ring_List)
        self.Name = ""
        self.OOP_Rotation = 0
        self.Dih_Rotation = 0
        self.Atom_List = []
        self.Bond_List = []
        self.Angle_List = []
        self.Dihedral_List = []
        self.Improper_List = []
        self.Interring_Bond_List = []
        self.UnConverged = True
        for ring in self.Ring_List:
            self.Name = self.Name + ring.Name + "_"
        info = []
        for i in range(len(self.Ring_List)-1):
            print("Linking Now")
            center_position,bond_vector,new_ring_bonded_atom = self.Ring_List[i].Link_Rings(self.Ring_List[i+1])
            info.append((center_position,bond_vector,new_ring_bonded_atom))
            self.Ring_List[i].Link_Bonded_Atoms()
            self.Ring_List[i+1].Link_Bonded_Atoms()
        """for ring in self.Ring_List:
            center_position = info[i][0]
            bond_vector = info[i][1]
            new_ring_bonded_atom = info[i][2]
            ring.Update_Normal_Vector()
            New_X = self.Ring_List[i].Normal_Vector
            New_Y = bond_vector/np.linalg.norm(bond_vector)
            New_X = (New_X - np.dot(New_X,New_Y)*New_Y)/np.linalg.norm(New_X - np.dot(New_X,New_Y)*New_Y)
            New_Z = np.cross(New_X,New_Y)/np.linalg.norm(np.cross(New_X,New_Y))
            Rotation_Matrix = [New_X,New_Y,New_Z]
            ring.Rotate_Ring(Rotation_Matrix)"""

        for i in range(len(self.Ring_List)-1):
            """for j in range(i+1):
                self.Ring_List[j].Translate_Ring(center_position)"""
            for b_atom in self.Ring_List[i+1].Bonded_Atoms:
                if b_atom.Is_Linked and b_atom.Bonded_Ring == self.Ring_List[i]:
                    center = b_atom.Central_Atom.Position
                    center_atom = b_atom
            self.Ring_List[i+1].Translate_Ring(center)
            self.Ring_List[i+1].Translate_Ring(-1*center_atom.Interring_Bond_Atom.Central_Atom.Position)
            self.Ring_List[i+1].Translate_Ring(center_atom.Bonded_Vector)
            """deviation_angle = new_ring_bonded_atom.Check_Alignment()
            self.Ring_List[i+1].Rotate_Ring_Not_Interring([[1,0,0],[0,math.cos(deviation_angle),-math.sin(deviation_angle)],[0,math.sin(deviation_angle),math.cos(deviation_angle)]])"""
            #self.Ring_List[i+1].Translate_Ring(new_ring_bonded_atom.Bonded_Vector)
        self.N = 0
        for ring in self.Ring_List:
            self.N += len(ring.Atom_List)
            for b_atom in ring.Bonded_Atoms:
                if not b_atom.Is_Linked:
                    self.N += 1
        atom_id = 1
        for ring in self.Ring_List:
            for atom in ring.Atom_List:
                atom.Atom_ID = atom_id
                atom_id += 1
            for b_atom in ring.Bonded_Atoms:
                if not b_atom.Is_Linked:
                    b_atom.H_Atom.Atom_ID = atom_id
                    atom_id += 1
        for ring in self.Ring_List:
            ring.Add_Bond_List()
        for ring in self.Ring_List:
            for atom in ring.Atom_List:
                if atom.Element == "C":
                    atom.Assign_OPLS_ID()
                    self.Atom_List.append(atom)
            for atom in ring.Atom_List:
                if atom.Element != "C":
                    atom.Assign_OPLS_ID()
                    if atom.Element == "H":
                        print(atom.OPLS_Class)
                        print(atom.OPLS_Type)
                    self.Atom_List.append(atom)
                    if atom.Element == "H":
                        print(self.Atom_List[-1].OPLS_Class)
                        print(self.Atom_List[-1].OPLS_Type)
            for b_atom in ring.Bonded_Atoms:
                if not b_atom.Is_Linked:
                    b_atom.H_Atom.Assign_OPLS_ID()
                    self.Atom_List.append(b_atom.H_Atom)
        for ring in self.Ring_List:
            for atom in ring.Atom_List:
                if atom.Element == "H":
                    print("Atom ID: %d Bonded OPLS Class: %d" % (atom.Atom_ID,atom.Bond_List[0].OPLS_Class))
        for i in range(len(self.Ring_List) - 1):
            for b_atom in self.Ring_List[i].Bonded_Atoms:
                if b_atom.Is_Linked and b_atom.Bonded_Ring == self.Ring_List[i+1]:
                    bond_master = b_atom.Central_Atom
            for b_atom in self.Ring_List[i+1].Bonded_Atoms:
                if b_atom.Is_Linked and b_atom.Bonded_Ring == self.Ring_List[i]:
                    bond_slave = b_atom.Central_Atom
            self.Interring_Bond_List.append(Bond.Bond(bond_master,bond_slave,bond_master.Position - bond_slave.Position))
        for ring in self.Ring_List:
            for bond in ring.Bond_List:
                self.Bond_List.append(bond)
            for b_atom in ring.Bonded_Atoms:
                if not b_atom.Is_Linked:
                    self.Bond_List.append(Bond.Bond(b_atom.Central_Atom,b_atom.H_Atom,np.linalg.norm(b_atom.H_Bond_Vector)))
        for bond in self.Interring_Bond_List:
            self.Bond_List.append(bond)
        self.Map_From_Bonds()
        OPLS.Assign_OPLS(self)
        for atom in self.Atom_List:
            print("Atom ID: %d OPLS Type: %d OPLS Class: %d" % (atom.Atom_ID,atom.OPLS_Type,atom.OPLS_Class))
        print("END OF CP INITIALIZATION")



    def Write_XYZ(self):
        f = open("%sPhi_%d_Theta_%d.xyz" % (self.Name,self.OOP_Rotation,self.Dih_Rotation ),'w')
        f.write("%d\n\n" % (self.N))
        i = 1
        for ring in self.Ring_List:
            for atom in ring.Atom_List:
                f.write("%s\t%f\t%f\t%f\n" % (atom.Element,atom.Position[0],atom.Position[1],atom.Position[2]))
            for b_atom in ring.Bonded_Atoms:
                if not b_atom.Is_Linked:
                    f.write("%s\t%f\t%f\t%f\n" % ("H",b_atom.H_Atom.Position[0],b_atom.H_Atom.Position[1],b_atom.H_Atom.Position[2]))
        f.close()

    def Map_From_Bonds(self):
        for atom1 in self.Atom_List:
            for atom2 in atom1.Bond_List:
                for atom3 in atom2.Bond_List:
                    if atom1 != atom3:
                        Identical_Angle = False
                        for angle in self.Angle_List:
                            if angle.Compare_Angles(atom1,atom2,atom3):
                                Identical_Angle = True
                                break
                        if not Identical_Angle:
                            self.Angle_List.append(Angle.Angle(atom2,atom1,atom3,0.0))
                        for atom4 in atom3.Bond_List:
                            if atom1 != atom4 and atom2 != atom4:
                                Identical_Dihedral = False
                                for dih in self.Dihedral_List:
                                    if dih.Compare_Dihedrals(atom1,atom2,atom3,atom4):
                                        Identical_Dihedral = True
                                        break
                                if not Identical_Dihedral:
                                    self.Dihedral_List.append(Dihedral.Dihedral(atom2,atom3,atom1,atom4,0.0))
        for ring in self.Ring_List:
            for atom1 in ring.Core_Atom_List:
                for atom2 in atom1.Bond_List:
                    if atom2 in ring.Core_Atom_List:
                        for atom3 in atom2.Bond_List:
                            if atom1 != atom3 and atom3 in ring.Core_Atom_List:
                                for atom4 in atom2.Bond_List:
                                    if atom4 != atom1 and atom4 != atom3 and atom4 not in ring.Core_Atom_List and atom4 in ring.Atom_List:
                                        Identical_Improper = False
                                        for imp in self.Improper_List:
                                            if imp.Compare_Impropers(atom1,atom2,atom3,atom4):
                                                Identical_Improper = True
                                                break
                                        if not Identical_Improper:
                                            self.Improper_List.append(Improper.Improper(atom2,atom3,atom1,atom4,3.5,2,len(self.Improper_List)+1))









