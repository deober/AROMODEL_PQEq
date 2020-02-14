#! usr/bin/python

# Open relevant modules

import Bond
import numpy as np
import copy
import Molecule
import os
import subprocess
import Write_Inputs
import Write_Submit_Script
import Bonded_Atom
import Cluster_IO

class Ring(Molecule.Molecule):
   
    def __init__(self,Supermolecule,Atom_Numbers,Core_Atom_List,Bond_Atoms,Bond_Atom_Vectors,Name,ID,Symmetric = True,Scheduler = "Torque",Cluster_Location='/oasis/tscc/scratch/andrewk/Optimized_Monomers'):
        self.Atom_List = []
        self.Core_Atom_List = []
        self.Name = Name
        self.Ring_ID = ID
        self.Symmetric = Symmetric
        for atom_num in Atom_Numbers:
            self.Atom_List.append(copy.deepcopy(Supermolecule.Get_Atom(atom_num)))
            if atom_num in Core_Atom_List:
                self.Core_Atom_List.append(self.Atom_List[-1])

        #Calculate Normal Vector using non-hydrogen atoms

        self.Normal_Vector = np.array([0.0,0.0,0.0])
        for cut1,atom1 in enumerate(self.Core_Atom_List):
            for cut2,atom2 in enumerate(self.Core_Atom_List[cut1+1:]):
                for atom3 in self.Core_Atom_List[cut1+cut2+2:]:
                    vec1 = atom1.Position - atom2.Position
                    vec2 = atom2.Position - atom3.Position
                    self.Normal_Vector = self.Normal_Vector + (np.cross(vec1,vec2)/np.linalg.norm(np.cross(vec1,vec2)))
        self.Normal_Vector = self.Normal_Vector/np.linalg.norm(self.Normal_Vector)

        #Link hydrogens to Bond_Atoms
        self.Bonded_Atoms = np.empty(len(Bond_Atoms),dtype = object)
        max_atom_id = len(Supermolecule.Atom_List)
        for i,atoms,bond_atom_vec in zip(range(len(Bond_Atoms)),Bond_Atoms,Bond_Atom_Vectors):
            Same_Ring_Atom_List = atoms[2:]
            self.Bonded_Atoms[i] = Bonded_Atom.Bonded_Atom(self.Get_Atom(atoms[0]),bond_atom_vec,atoms[1],(max_atom_id + i + 1),self,Same_Ring_Atom_List)

        new_atom_id = 1
        for atom in self.Atom_List:
            atom.Atom_ID = new_atom_id 
            new_atom_id += 1
        for b_atom in self.Bonded_Atoms:
            b_atom.H_Atom.Atom_ID = new_atom_id
            new_atom_id += 1

        #Optimize geometry
        self.Optimize_H_Positions(Cluster_Location,Shared_File_Location = "/Users/andrewkleinschmidt/Shared_Files_Dihedral_Parameterization/")

    def Find_Coul_Cutoff(self):
        coul_cutoff = 0
        for atom1 in self.Atom_List:
            for atom2 in self.Atom_List:
                if atom1 != atom2 and np.linalg.norm(atom1.Position - atom2.Position) > coul_cutoff:
                    coul_cutoff = np.linalg.norm(atom1.Position - atom2.Position)
        coul_cutoff += 2.5
        return coul_cutoff

    def Update_Positions_Data_File(self,Updated_Data_File):
        f = open(Updated_Data_File,'r')
        lines = f.readlines()
        bond_atom_cutoff = (-1*len(self.Bonded_Atoms))
        regular_atom_lines = lines[2:bond_atom_cutoff]
        bonded_atom_lines = lines[bond_atom_cutoff:]
        if len(self.Atom_List) != len(regular_atom_lines):
            raise Exception("Length mismatch in regular atoms and xyz file")
        for atom,line in zip(self.Atom_List,regular_atom_lines):
            atom.Position = np.array([float(line.strip().split()[1]),float(line.strip().split()[2]),float(line.strip().split()[3])])
        if len(self.Bonded_Atoms) != len(bonded_atom_lines):
            raise Exception("Length mismatch in bonded atoms and xyz file")
        for b_atom,line in zip(self.Bonded_Atoms,bonded_atom_lines):
            b_atom.H_Atom.Position = np.array([float(line.strip().split()[1]),float(line.strip().split()[2]),float(line.strip().split()[3])])
            b_atom.H_Bond_Vector = b_atom.H_Atom.Position - b_atom.Central_Atom.Position

    def Update_Positions_Orca_Output(self,Orca_Output_File):
        f = open(Orca_Output_File,'r')
        lines = f.readlines()
        read_atoms = False
        off_count = 0
        for line in lines:
            if read_atoms and len(line.strip().split()) == 4:
                atom_list.append(line)
            elif read_atoms:
                off_count +=1
            if off_count >= 2:
                read_atoms = False
            if len(line.strip().split()) >= 3 and line.strip().split()[0] == "CARTESIAN" and line.strip().split()[1] == "COORDINATES" and line.strip().split()[2] == "(ANGSTROEM)":
                read_atoms = True
                off_count = 0
                atom_list = []
        bond_atom_cutoff = (-1*len(self.Bonded_Atoms))
        regular_atom_lines = atom_list[:bond_atom_cutoff]
        bonded_atom_lines = atom_list[bond_atom_cutoff:]
        if len(self.Atom_List) != len(regular_atom_lines):
            raise Exception("Length mismatch in regular atoms and xyz file")
        for atom,line in zip(self.Atom_List,regular_atom_lines):
            atom.Position = np.array([float(line.strip().split()[1]),float(line.strip().split()[2]),float(line.strip().split()[3])])
        if len(self.Bonded_Atoms) != len(bonded_atom_lines):
            raise Exception("Length mismatch in bonded atoms and xyz file")
        for b_atom,line in zip(self.Bonded_Atoms,bonded_atom_lines):
            b_atom.H_Atom.Position = np.array([float(line.strip().split()[1]),float(line.strip().split()[2]),float(line.strip().split()[3])])
            b_atom.H_Bond_Vector = b_atom.H_Atom.Position - b_atom.Central_Atom.Position
        self.Update_Normal_Vector()

    def Update_Positions_Data_File(self,Updated_Data_File):
        f = open(Updated_Data_File,'r')
        lines = f.readlines()
        bond_atom_cutoff = (-1*len(self.Bonded_Atoms))
        regular_atom_lines = lines[2:bond_atom_cutoff]
        bonded_atom_lines = lines[bond_atom_cutoff:]
        if len(self.Atom_List) != len(regular_atom_lines):
            raise Exception("Length mismatch in regular atoms and xyz file")
        for atom,line in zip(self.Atom_List,regular_atom_lines):
            atom.Position = np.array([float(line.strip().split()[1]),float(line.strip().split()[2]),float(line.strip().split()[3])])
        if len(self.Bonded_Atoms) != len(bonded_atom_lines):
            raise Exception("Length mismatch in bonded atoms and xyz file")
        for b_atom,line in zip(self.Bonded_Atoms,bonded_atom_lines):
            b_atom.H_Atom.Position = np.array([float(line.strip().split()[1]),float(line.strip().split()[2]),float(line.strip().split()[3])])
            b_atom.H_Bond_Vector = b_atom.H_Atom.Position - b_atom.Central_Atom.Position
        self.Update_Normal_Vector()

    def Set_Up_Bonds(self,Orca_Output_File):
        f = open(Orca_Output_File,'r')
        lines = f.readlines()
        read_bonds = False
        for line in lines:
            if read_bonds and len(line.strip().split()) >= 1:
                for bond in line.strip().split("B("):
                    if len(bond.strip().split()) >= 3 and float(bond.strip().split()[-1]) > 0.5:
                        Master_ID = int(bond.strip().split()[0].split('-')[0]) + 1
                        if any(atom.Atom_ID == Master_ID for atom in self.Atom_List):
                            Master_Atom = self.Get_Atom(Master_ID)
                            Slave_ID = int(bond.strip().split(',')[1].split()[0].split('-')[0].strip(")")) + 1
                            if any(atom.Atom_ID == Slave_ID for atom in self.Atom_List):
                                Slave_Atom = self.Get_Atom(Slave_ID)
                                req = np.linalg.norm(Master_Atom.Position - Slave_Atom.Position)
                                self.Bond_List.append(Bond.Bond(Master_Atom,Slave_Atom,req))
            elif read_bonds:
                read_bonds = False
            if len(line.strip().split()) >= 3 and line.strip().split()[0] == "Mayer" and line.strip().split()[1] == "bond" and line.strip().split()[2] == "orders":
                read_bonds = True
                self.Bond_List = []
    def Optimize_H_Positions(self,Cluster_Location,Shared_File_Location = ""):
        f = open("%s.xyz" % self.Name,'w')
        num_atoms = len(self.Atom_List) + len(self.Bonded_Atoms)
        f.write("%d\n\n" % (num_atoms))
        for atom in self.Atom_List:
            f.write("%s\t%f\t%f\t%f\n" % (atom.Element,atom.Position[0],atom.Position[1],atom.Position[2]))
        for atom in self.Bonded_Atoms:
            f.write("%s\t%f\t%f\t%f\n" % ("H",atom.H_Atom.Position[0],atom.H_Atom.Position[1],atom.H_Atom.Position[2]))
        f.close()
        Monomer = Molecule.Molecule("%s.xyz" % self.Name)
        os.system("mkdir ./Optimized_Monomers")
        File_Name = "sub_%s_Optimize_Monomer" % self.Name
        In_File = "%s_Optimize_Monomer.inp" % self.Name
        End_File = "%s_Optimize_Monomer.out" % self.Name
        Job_Type = "Orca"
        Folder_Name = "Optimized_Monomers"
        Job_Name = "%s_Optimize_Monomer" % self.Name
        Cluster_Login = "andrewk@tscc-login.sdsc.edu"
        Base_Cluster_Location = '/oasis/tscc/scratch/andrewk/'
        Scheduler_Type = "TORQUE"
        End_Condition = "Opt_Orca"
        Executable_Location = "/home/andrewk/orca_4_2_0_linux_x86-64_openmpi314"
        OpenMP_Location = "/home/andrewk/openmpi-3.1.4"

        Write_Inputs.Write_Orca_Optimize_Geometry(In_File,Monomer,H_Only = True)
        Write_Submit_Script.Write_TORQUE(File_Name,In_File,Job_Name,1,Cluster_Location,Job_Type,Executable_Path = Executable_Location,OMP_Path = OpenMP_Location)
        Copy_File_List = [File_Name,In_File]
        Cluster_IO.Submit_Job(Copy_File_List,Folder_Name,File_Name,End_File,Job_Name,Cluster_Login,Cluster_Location,Base_Cluster_Location,Scheduler_Type,End_Condition = End_Condition,Analyze_File = End_File,Shared_File_Location = Shared_File_Location)
        Cluster_IO.Return_Info(End_File,End_File,Folder_Name,Job_Type,Cluster_Login,Cluster_Location,End_Condition = End_Condition,Shared_File_Location = Shared_File_Location)
        self.Update_Positions_Orca_Output("./Optimized_Monomers/%s" % End_File)
        self.Set_Up_Bonds("./Optimized_Monomers/%s" % End_File)


        
        #self.Update_H_Positions_Orca("./Optimized_Monomers/%s" % Log_File)
        return "%s.xyz" % self.Name

    #Link Bond_Atoms to other rings and calculate angles with other rings
    def Link_Rings(self,Link_Ring):
        Linked = False
        for b_atom in Link_Ring.Bonded_Atoms:
            for b_atom2 in self.Bonded_Atoms:
                if not b_atom.Is_Linked and not b_atom2.Is_Linked and b_atom.Bonded_Vector[0] == -1*b_atom2.Bonded_Vector[0] and b_atom.Bonded_Vector[1] == -1*b_atom2.Bonded_Vector[1] and b_atom.Bonded_Vector[2] == -1*b_atom2.Bonded_Vector[2]:
                    b_atom.Add_Ring(self)
                    b_atom2.Add_Ring(Link_Ring)
                    center_position = b_atom2.Central_Atom.Position
                    bond_vector = b_atom.Bonded_Vector
                    Linked = True
                    break
                else:
                    continue
            break

        if not Linked:
            alignment = 0
            for b_atom in Link_Ring.Bonded_Atoms:
                for b_atom2 in self.Bonded_Atoms:
                    if not b_atom.Is_Linked and not b_atom2.Is_Linked and np.dot(b_atom.Bonded_Vector/np.linalg.norm(b_atom.Bonded_Vector),-1*b_atom2.Bonded_Vector/np.linalg.norm(b_atom2.Bonded_Vector)) > alignment:
                        alignment = np.dot(b_atom.Bonded_Vector,-1*b_atom2.Bonded_Vector)
                        temp_b_atom = b_atom
                        temp_b_atom2 = b_atom2
                        center_position = b_atom2.Central_Atom.Position
                        bond_vector = b_atom.Bonded_Vector
            temp_b_atom.Add_Ring(self)
            temp_b_atom2.Add_Ring(Link_Ring)
            b_atom = b_atom
       
            if alignment == 0:
                print("Ring 1:")
                for b_atom in Link_Ring.Bonded_Atoms:
                    print(b_atom.Bonded_Vector)
                    print(b_atom.Is_Linked)
                print("Ring 2:")
                for b_atom in self.Bonded_Atoms:
                    print(b_atom.Bonded_Vector)
                    print(b_atom.Is_Linked)
                raise Exception("Rings not properly linked")

        return center_position,bond_vector,b_atom

    def Link_Bonded_Atoms(self):
        for b_atom in self.Bonded_Atoms:
            if b_atom.Is_Linked:
                for b_atom_2 in b_atom.Bonded_Ring.Bonded_Atoms:
                    if b_atom_2.Is_Linked and b_atom_2.Bonded_Ring == self:
                        b_atom.Interring_Bond_Atom = b_atom_2
        for b_atom in self.Bonded_Atoms:
            if b_atom.Is_Linked:
                print(b_atom.Interring_Bond_Atom.Central_Atom.Atom_ID)

    #Update Normal Vector
    def Update_Normal_Vector(self):
        self.Normal_Vector = np.array([0.0,0.0,0.0])
        for cut1,atom1 in enumerate(self.Core_Atom_List):
            for cut2,atom2 in enumerate(self.Core_Atom_List[cut1+1:]):
                for atom3 in self.Core_Atom_List[cut1+cut2+2:]:
                    vec1 = atom1.Position - atom2.Position
                    vec2 = atom2.Position - atom3.Position
                    self.Normal_Vector = self.Normal_Vector + (np.cross(vec1,vec2)/np.linalg.norm(np.cross(vec1,vec2)))
        self.Normal_Vector = self.Normal_Vector/np.linalg.norm(self.Normal_Vector)
    
    def Rotate_Ring(self,Rotation_Matrix):
        for atom in self.Atom_List:
            atom.Position = np.matmul(Rotation_Matrix,atom.Position)
        for b_atom in self.Bonded_Atoms:
            b_atom.H_Atom.Position = np.matmul(Rotation_Matrix,b_atom.H_Atom.Position)
            b_atom.H_Bond_Vector = np.matmul(Rotation_Matrix,b_atom.H_Bond_Vector)
            b_atom.Bonded_Vector = np.matmul(Rotation_Matrix,b_atom.Bonded_Vector)

    def Rotate_Ring_Not_Interring(self,Rotation_Matrix):
        for atom in self.Atom_List:
            atom.Position = np.matmul(Rotation_Matrix,atom.Position)
        for b_atom in self.Bonded_Atoms:
            b_atom.H_Atom.Position = np.matmul(Rotation_Matrix,b_atom.H_Atom.Position)
            b_atom.H_Bond_Vector = np.matmul(Rotation_Matrix,b_atom.H_Bond_Vector)

    def Translate_Ring(self,Center_Position):
        for atom in self.Atom_List:
            atom.Position = atom.Position - Center_Position
        for b_atom in self.Bonded_Atoms:
            b_atom.H_Atom.Position = b_atom.H_Atom.Position - Center_Position     

    def Add_Bond_List(self):
        for atom in self.Atom_List:
            for bond in self.Bond_List:
                if bond.Bond_Master == atom:
                    atom.Bond_List.append(bond.Bond_Slave)
                elif bond.Bond_Slave == atom:
                    atom.Bond_List.append(bond.Bond_Master)
                if len(atom.Bond_List) == 4 or (atom.Element == "H" and len(atom.Bond_List) == 1):
                    break
        for b_atom in self.Bonded_Atoms:
            if b_atom.Is_Linked:
                b_atom.Central_Atom.Bond_List.append(b_atom.Interring_Bond_Atom.Central_Atom)
            else:
                b_atom.Central_Atom.Bond_List.append(b_atom.H_Atom)
                b_atom.H_Atom.Bond_List.append(b_atom.Central_Atom)

