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

class Ring(Molecule):
    
    def Return_Monomer(self):
        f = open("%s.xyz" % self.Name,'w')
        num_atoms = len(self.Atom_List) + len(self.Bonded_Atoms)
        f.write("%d\n\n" % (num_atoms))
        for atom in self.Atom_List:
            f.write("%s\t%f\t%f\t%f\n" % (atom.Element,atom.Position[0],atom.Position[1],atom.Position[2]))
        f.close()
        return "%s.xyz" % self.Name
   
    def __init__(self,Supermolecule,Atom_Numbers,Bond_Atoms,Bond_Atom_Vectors,Name,ID,Symmetric = True,Scheduler = "Torque",Cluster_Location='/oasis/tscc/scratch/andrewk/Optimized_Monomers'):
        self.Atom_List = []
        self.Core_Atom_List = []
        self.Name = Name
        self.Ring_ID = ID
        self.Symmetric = Symmetric
        for atom_num in Atom_Numbers:
            self.Atom_List.append(copy.deepcopy(Supermolecule.Get_Atom(atom_num)))
            if Supermolecule.Get_Atom(atom_num).Element != "H":
                self.Core_Atom_List.append(copy.deepcopy(Supermolecule.Get_Atom(atom_num)))

        #Calculate Normal Vector using non-hydrogen atoms

        self.Normal_Vector = np.array([0.0,0.0,0.0])
        for cut1,atom1 in enumerate(self.Core_Atom_List):
            for cut2,atom2 in enumerate(self.Core_Atom_List[cut1+1:]):
                for atom3 in self.Core_Atom_List[cut1+cut2+2]:
                    vec1 = atom1.Position - atom2.Position
                    vec2 = atom2.Position - atom3.Position
                    self.Normal_Vector = self.Normal_Vector + (np.cross(vec1,vec2)/np.linalg.norm(np.cross(vec1,vec2)))
        self.Normal_Vector = self.Normal_Vector/np.linalg.norm(self.Normal_Vector)

        #Link hydrogens to Bond_Atoms
        self.Bonded_Atoms = np.empty(len(self.Atom_List),dtype = object)
        max_atom_id = len(Supermolecule.Atom_List)
        for i,atom,bond_atom_vec in zip(range(len(Bond_Atoms)),Bond_Atoms,Bond_Atom_Vectors):
            self.Bonded_Atoms[i] = Bonded_Atom.Bonded_Atom(self.Get_Atom(atom),bond_atom_vec,max_atom_id + i + 1,self)

        #Optimize geometry
        In_File = "in.%s_Optimize_Geometry" % self.Name
        Write_Inputs.Write_LAMMPS_Minimize(In_File,Molecule.Molecule(self.Return_Monomer()),self.Name,Coul_Cutoff)

        """Write_Inputs.Write_Orca_Optimize_Geometry(In_File,Molecule.Molecule(self.Return_Monomer()))
        File_Name = "sub_%s_Geometry_Optimize" % self.Name
        nproc = 4
        Cluster_Location = ""
        Job_Type = Orca
        if Scheduler == "Torque":
            Write_Submit_Script.Write_Torque(File_Name,In_File,self.Name,nproc,Cluster_Location,Job_Type,proc_per_node=28,walltime = 2,Executable_Path = "/home/andrewk/orca_4_2_0_linux_x86-64_openmpi314",OMP_Path = "/home/andrewk/openmpi-3.1.4",queue = "condo")
        else:
            Write_Submit_Script.Write_SLURM(File_Name,In_File,self.Name,nproc,Cluster_Location,Job_Type,queue="shared",proc_per_node=32,account = "m3047",walltime = 2,Executable_Path = "/home/andrewk/orca_4_2_0_linux_x86-64_openmpi314",OMP_Path = "/home/andrewk/openmpi-3.1.4",constraint='cori')
        Copy_File_List = [File_Name,In_File]
        Folder_Name = 'Optimized_Monomers'
        Shared_File_Location = '~/Shared_Files_Dihedral_Parameterization/Optimized_Monomers'
        Submit_Script = File_Name
        End_File = "%s.out" % self.Name
        Cluster_Login = "andrewk@tscc-login.sdsc.edu"
        Base_Cluster_Location = '/oasis/tscc/scratch/andrewk/'
        Scheduler_Type = "TORQUE"
        Analyze_File = End_File
        Cluster_IO.Submit_Job(Copy_File_List,Folder_Name,Submit_Script,End_File,Job_Name,Cluster_Login,Cluster_Location,Base_Cluster_Location,Scheduler_Type,Shared_File_Location = Shared_File_Location,End_Condition = "Opt_Orca")
        Cluster_IO.Return_Info(Analyze_File,End_File,Folder_Name,Cluster_Login,Cluster_Location,Shared_File_Location = Shared_File_Location,End_Condition = Opt_Orca)"""
    #Link Bond_Atoms to other rings and calculate angles with other rings

    #Update Normal Vector
    