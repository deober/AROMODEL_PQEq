#! usr/bin/python

# Import relevant modules
import numpy as np
import os
import subprocess
import shutil
import time
import random
import math

# import Class structure
import Atom
import Bond
import Angle
import Dihedral
import Improper
import Configure

from operator import itemgetter

Element_Dict = { 12.011:"C", 1.008:"H", 18.998:"F", 15.999:"O", 32.06:"S", 14.007:"N", 28.086:"Si", 35.453:"Cl"}

class DA_Polymer(object):
    """
        Class defining a DA Polymer from 2015 JACS "Conformational
        Order in Aggregates of ..."
        
        instance variables (self) :
        N = Number of Atoms
        Name = Molecule Name
        Atom_list[N] = List of Atom objects
        Bond_List[Num_Bonds] = list of bond objects
        Angle_List[Num_Angles] = list of angle objects
        Dihedral_List[Num_Dihedrals] = list of Dihedral objects
        Improper_List[Num_Impropers] = list of Improper objects
        MW = Molecular Weight
        COM = Center of Mass
        MOL_ID = ID Number for molecule
        unConverged = Flag for bypassing Orca convergence (Default = False)
        """
    def __init__(self, Filename, Anneal=False, Full = False):
        
        
        self.Name = Filename.split('.')[0].split('_')[0]
        File = open(Filename,'r')
        File_Lines = File.readlines()
        self.Bond_List = []
        self.Angle_List = []
        self.Dihedral_List = []
        self.Improper_List = []
        self.MW = 0.0
        self.COM = np.zeros(3,float)
        self.Mol_ID =0
        self.basis = np.eye(3)
        self.Box_Size = 0

        j = -1
        for Line in File_Lines:
            j += 1
            Line = Line.strip('\n').split(' ')
            #print(Line)
            if len(Line) > 1:
                if len(Line) > 2:
                    if Line[2] == 'xlo':
                        #print(Line[1])
                        #self.Box_Size = float(Line[1].split('e')[0])*(10**(int(float(Line[1].split('+')[1]))))
                        self.Box_Size = float(Line[1])
                if Line[1] == 'atoms':
                    self.N = int(Line[0])
                    #print self.N, "Atoms"
                if Line[1] == 'atom':
                    Atom_Types = int(Line[0])
                    self.Atom_List = np.empty(self.N, dtype=object)
                    self.I_Flags_List = np.empty(self.N, dtype=object)
                    Mass_List = np.empty(Atom_Types, dtype=float)
                    Pair_Coeffs = np.empty((Atom_Types,2), dtype=float)
                if Line[1] == 'bonds':
                    Num_Bonds = Line[0]
                if Line[1] == 'bond':
                    Bond_Types = int(Line[0])
                    Bond_Coeffs = np.empty((Bond_Types,2), dtype=float)
                if Line[1] == 'angles':
                    Num_Angles = Line[0]
                if Line[1] == 'angle':
                    Angle_Types = int(Line[0])
                    Angle_Coeffs = np.empty((Angle_Types,2),dtype=float)
                if Line[1] == 'dihedrals':
                    Num_Dihedrals = Line[0]
                if Line[1] == 'dihedral':
                    Dihedral_Types = int(Line[0])
                    #print('dihedral types')
                    #print (Dihedral_Types)
                    #Dihedral_Coeffs = [np.empty((Dihedral_Types,9),dtype=float)]
                    Dihedral_Coeffs = []
                    Dihedral_Styles = []
                if Line[1] == 'impropers':
                    Num_Impropers = Line[0]
                if Line[1] == 'improper':
                    Improper_Types = int(Line[0])
                    Improper_Coeffs = np.empty((Improper_Types,3),dtype=float)
                if Line[0] == 'Atoms':
                    #print 'Getting Atoms'
                    for i in range(self.N):
                        Temp_Position = [float(File_Lines[i+j+2].split(' ')[4]), float(File_Lines[i+j+2].split(' ')[5]), float(File_Lines[i+j+2].split(' ')[6]) ]
                        try:
                            Temp_I_Flags = [ int(File_Lines[i+j+2].split(' ')[7]), int(File_Lines[i+j+2].split(' ')[8]), int(File_Lines[i+j+2].split(' ')[9])]
                            Temp_I_Flags = np.asarray(Temp_I_Flags, dtype = int)
                            Temp_Position[0] = Temp_Position[0] + Temp_I_Flags[0] * self.Box_Size
                            Temp_Position[1] = Temp_Position[1] + Temp_I_Flags[1] * self.Box_Size
                            Temp_Position[2] = Temp_Position[2] + Temp_I_Flags[2] * self.Box_Size
                            self.I_Flags_List[i] = Temp_I_Flags
                        except:
                            print("No Image Flags")
                        Temp_Position = np.asarray(Temp_Position, dtype=float)

                        #I_Flags = np.asarray(Temp_I_Flags, dtype=int)
                        #print "IMAGE FLAGS", I_Flags
                        #Temp_Position += I_Flags*Box_Length
                        Temp_Charge = float(File_Lines[i+j+2].split(' ')[3])
                        Type = int(File_Lines[i+j+2].split(' ')[2])
                        Mass = Mass_List[Type-1]
                        Element = Element_Dict[Mass]
                        Atom_Id = int(File_Lines[i+j+2].split(' ')[0])
                        #print Element, Mass, Atom_Id
                        self.Atom_List[i] = Atom.Atom(Temp_Position, Element, Atom_Id)
                        self.Atom_List[i].Charge = Temp_Charge
                        self.Atom_List[i].Epsilon = Pair_Coeffs[Type-1,0]
                        self.Atom_List[i].Sigma = Pair_Coeffs[Type-1,1]
                        if Full:
                            self.Atom_List[i].OPLS_Type = 2000 + Type
                            self.Atom_List[i].LAMMPS_Type = 2000 + Type
                        else:
                            self.Atom_List[i].OPLS_Type = Type
                            self.Atom_List[i].LAMMPS_Type = Type
                    self.Atom_List = sorted(self.Atom_List, key=lambda AtomO: AtomO.Atom_ID)
                    #for Atom_Obj in self.Atom_List:
                        #print Atom_Obj.Atom_ID, Atom_Obj.Position
                try:
                    if Line[2] == 'xlo':
                        Box_Length = float(Line[1]) - float(Line[0])
                        #print Box_Length, "Box_Length"
                except:
                    continue

                if Line[0] == 'Pair':
                    #print "Getting Pair Coefficients"
                    for i in range(Atom_Types):
                        Pair_Coeffs[i,0] = float(File_Lines[i+j+2].split(' ')[1])
                        Pair_Coeffs[i,1] = float(File_Lines[i+j+2].split(' ')[2])
                if Line[0] == 'Bond':
                    #print "Getting Bond Coefficients"
                    for i in range(Bond_Types):
                        Bond_Coeffs[i,0] = float(File_Lines[i+j+2].split(' ')[1])
                        Bond_Coeffs[i,1] = float(File_Lines[i+j+2].split(' ')[2])
                if Line[0] == 'Angle':
                    #print "Getting Angle Coefficients"
                    for i in range(Angle_Types):
                        Angle_Coeffs[i,0] = float(File_Lines[i+j+2].split(' ')[1])
                        Angle_Coeffs[i,1] = float(File_Lines[i+j+2].split(' ')[2])
                if Line[0] == 'Dihedral':
                    #print "Getting Dihedral Coefficients"
                    #print(len(Dihedral_Coeffs))
                    for i in range(Dihedral_Types):
                        dummy_coeffs = []
                        for k in range(1,len(File_Lines[i+j+2].split(' '))):
                            dummy_coeffs.append(float(File_Lines[i+j+2].split(' ')[k]))
                        print(dummy_coeffs)
                        Dihedral_Coeffs.append(np.array(dummy_coeffs))
                        Dihedral_Styles.append(File_Lines[i+j+2].split(' ')[0])
                            #Dihedral_Coeffs[i,k-2] = float(File_Lines[i+j+2].split(' ')[k])
                        
                        """Dihedral_Coeffs[i,1] = float(File_Lines[i+j+2].split(' ')[2])
                        Dihedral_Coeffs[i,2] = float(File_Lines[i+j+2].split(' ')[3])
                        Dihedral_Coeffs[i,3] = float(File_Lines[i+j+2].split(' ')[4])
                        Dihedral_Coeffs[i,4] = float(File_Lines[i+j+2].split(' ')[5])
                        Dihedral_Coeffs[i,5] = float(File_Lines[i+j+2].split(' ')[6])
                        Dihedral_Coeffs[i,6] = float(File_Lines[i+j+2].split(' ')[7])
                        Dihedral_Coeffs[i,7] = float(File_Lines[i+j+2].split(' ')[8])
                        Dihedral_Coeffs[i,8] = float(File_Lines[i+j+2].split(' ')[9])"""
                        
                if Line[0] == 'Improper':
                    #print "Getting Improper Coefficients"
                    for i in range(Improper_Types):
                        Improper_Coeffs[i,0] = float(File_Lines[i+j+2].split(' ')[1])
                        Improper_Coeffs[i,1] = int(File_Lines[i+j+2].split(' ')[2])
                        Improper_Coeffs[i,2] = int(File_Lines[i+j+2].split(' ')[3])



            if len(Line) == 1:
                if Line == ['']:
                    continue
                if Line == ['Masses']:
                    #print "Extracting Masses"
                    for i in range(Atom_Types):
                        Mass_List[i] = float(File_Lines[i+j+2].split(' ')[1])
                if Line == ['Bonds']:
                    #print "Extracting Bonds"
                    for i in range(int(Num_Bonds)):
                        Bond_Info = File_Lines[i+j+2].strip('\n').split(' ')
                        Master = self.Atom_List[int(Bond_Info[2])-1]
                        Slave = self.Atom_List[int(Bond_Info[3])-1]
                        Kb = Bond_Coeffs[int(Bond_Info[1])-1, 0]
                        Req = Bond_Coeffs[int(Bond_Info[1])-1,1]
                        Bond_ID = int(Bond_Info[0])
                        self.Bond_List.append(Bond.Bond(Master, Slave,Req))
                        self.Bond_List[i].kb = Kb
                        self.Bond_List[i].Bond_ID = Bond_ID
                        self.Bond_List[i].LAMMPS_Type = Bond_ID = int(Bond_Info[1])
                if Line == ['Angles']:
                    #print "Extracting Angles"
                    for i in range(int(Num_Angles)):
                        Angle_Info = File_Lines[i+j+2].strip('\n').split(' ')
                        Slave1 = self.Atom_List[int(Angle_Info[2])-1]
                        Master = self.Atom_List[int(Angle_Info[3])-1]
                        Slave2 = self.Atom_List[int(Angle_Info[4])-1]
                        Ka = Angle_Coeffs[int(Angle_Info[1])-1,0]
                        Th0 = Angle_Coeffs[int(Angle_Info[1])-1,1]
                        Angle_ID = int(Angle_Info[0])
                        self.Angle_List.append(Angle.Angle(Master, Slave1, Slave2, Th0))
                        self.Angle_List[i].ka = Ka
                        self.Angle_List[i].Angle_ID = Angle_ID
                        self.Angle_List[i].LAMMPS_Type = int(Angle_Info[1])
                print('Check')
                if Line == ['Dihedrals']:
                    #print "Extracting Dihedrals"
                    for i in range(int(Num_Dihedrals)):
                        Dihedral_Info = File_Lines[i+j+2].strip('\n').split(' ')
                        Slave1 = self.Atom_List[int(Dihedral_Info[2])-1]
                        Master1 = self.Atom_List[int(Dihedral_Info[3])-1]
                        Master2 = self.Atom_List[int(Dihedral_Info[4])-1]
                        Slave2 = self.Atom_List[int(Dihedral_Info[5])-1]
                        Coeffs = Dihedral_Coeffs[int(Dihedral_Info[1])-1]
                        #Coeffs = []
                        #Style = Dihedral_Styles[int(Dihedral_Info[1])-1]
                        Style = int(Dihedral_Info[1])
                        Dihedral_ID = int(Dihedral_Info[0])
                        self.Dihedral_List.append(Dihedral.Dihedral(Master1, Master2, Slave1, Slave2, 0.0))
                        self.Dihedral_List[i].Coeffs = Coeffs
                        self.Dihedral_List[i].Dihedral_ID = Dihedral_ID
                        self.Dihedral_List[i].Style = Style
                        self.Dihedral_List[i].LAMMPS_Type = Dihedral_Info[1]

                if Line == ['Impropers']:
                    #print "Extracting Impropers"
                    for i in range(int(Num_Impropers)):
                        Improper_Info = File_Lines[i+j+2].strip('\n').split(' ')
                        Master = self.Atom_List[int(Improper_Info[2])-1]
                        Slave1 = self.Atom_List[int(Improper_Info[3])-1]
                        Slave2 = self.Atom_List[int(Improper_Info[4])-1]
                        Slave3 = self.Atom_List[int(Improper_Info[5])-1]
                        Coeff = Improper_Coeffs[int(Improper_Info[1])-1,0]
                        Improper_ID = int(Improper_Info[0])
                        self.Improper_List.append(Improper.Improper(Master, Slave1, Slave2, Slave3, Coeff, 180.0, Improper_ID))
                        self.Improper_List[i].LAMMPS_Type = int(Improper_Info[1])
                        
                            
        # Compute center of Mass
        Mass_Weighted_Sum = np.zeros(3,dtype=float)
        #print len(self.Atom_List)
        for Atom_Obj in self.Atom_List:
            Mass_Weighted_Sum += Atom_Obj.Position*Atom_Obj.Mass
            self.MW += Atom_Obj.Mass
        
        #print "Molecular Weight is ", self.MW, "Grams/Mole"
        self.COM = Mass_Weighted_Sum/self.MW
        #print "COM is", self.COM
        
        # Zero COM
        for Atom_Obj in self.Atom_List:
            Atom_Obj.Position -= self.COM
        
        #self.COM -= self.COM
        #print self.COM
        
        
        return
    

    def Write_Data_File(self,Filename):
        f = open(Filename,'w')
        f.write("LAMMPS data file\n\n")
        f.write("%d atoms\n" % len(self.Atom_List))
        type_count = 0 
        type_list = []
        Atom_Info = []
        for atom in self.Atom_List:
            if atom.LAMMPS_Type not in type_list:
                type_count += 1
                type_list.append(atom.LAMMPS_Type)
                Atom_Info.append((atom.LAMMPS_Type,atom.Mass,atom.Epsilon,atom.Sigma))
        Atom_Info = sorted(Atom_Info, key=lambda atom_type: atom_type[0])

        f.write("%d atom types\n" % type_count)
        f.write("%d bonds\n" % len(self.Bond_List))
        type_count = 0 
        type_list = []
        Bond_Info = []
        for bond in self.Bond_List:
            if bond.LAMMPS_Type not in type_list:
                type_count += 1
                type_list.append(bond.LAMMPS_Type)
                Bond_Info.append((bond.LAMMPS_Type,bond.kb,bond.req))
        Bond_Info = sorted(Bond_Info, key=lambda bond_type: bond_type[0])
        f.write("%d bond types\n" % type_count)
        f.write("%d angles\n" % len(self.Angle_List))
        type_count = 0 
        type_list = []
        Angle_Info = []
        for ang in self.Angle_List:
            if ang.LAMMPS_Type not in type_list:
                type_count += 1
                type_list.append(ang.LAMMPS_Type)
                Angle_Info.append((ang.LAMMPS_Type,ang.ka,ang.Angle_Eq))
        Angle_Info = sorted(Angle_Info, key=lambda angle_type: angle_type[0])
        f.write("%d angle types\n" % type_count)
        f.write("%d dihedrals\n" % len(self.Dihedral_List))
        type_count = 0 
        type_list = []
        Dihedral_Info = []
        for dih in self.Dihedral_List:
            if dih.LAMMPS_Type not in type_list:
                type_count += 1
                type_list.append(dih.LAMMPS_Type)
                Dihedral_Info.append((int(dih.LAMMPS_Type),dih.Coeffs))
        Dihedral_Info = sorted(Dihedral_Info, key=lambda dih_type: dih_type[0])
        f.write("%d dihedral types\n" % type_count)
        f.write("%d impropers\n" % len(self.Improper_List))
        type_count = 0 
        type_list = []
        Improper_Info = []
        for imp in self.Improper_List:
            if imp.LAMMPS_Type not in type_list:
                type_count += 1
                type_list.append(imp.LAMMPS_Type)
                Improper_Info.append((int(imp.LAMMPS_Type),imp.Ki))
        Improper_Info = sorted(Improper_Info, key=lambda imp_type: imp_type[0])
        f.write("%d improper types\n\n" % type_count)
        f.write("\n\n\n0.0000 300.0000 xlo xhi\n0.0000 300.0000 ylo yhi\n0.0000 300.0000 zlo zhi\n\n")
        f.write("Masses\n\n")
        for atom in Atom_Info:
            f.write("%d %.3f\n" % (atom[0],atom[1]))
        f.write("\nPair Coeffs\n\n")
        for atom in Atom_Info:
            f.write("%d %f %f\n" % (atom[0],atom[2],atom[3]))
        f.write("\nBond Coeffs\n\n")
        for bond in Bond_Info:
            f.write("%d %f %f\n" % (bond[0],bond[1],bond[2]))
        f.write("\nAngle Coeffs\n\n")
        for ang in Angle_Info:
            f.write("%d %f %f\n" % (ang[0],ang[1],ang[2]))
        f.write("\nDihedral Coeffs\n\n")
        for dih in Dihedral_Info:
            f.write("%s" % dih[0])
            for c in dih[1]:
                f.write(" %s" % c)
            f.write("\n")
        f.write("\nImproper Coeffs\n\n")
        for imp in Improper_Info:
            f.write("%d %f -1 2\n" % (imp[0],imp[1]))
        f.write("\nAtoms\n\n")
        for atom in self.Atom_List:
            f.write("%d 1 %d %.16f %.16f %.16f %.16f %d %d %d\n" % (atom.Atom_ID,atom.LAMMPS_Type,atom.Charge,atom.Position[0],atom.Position[1],atom.Position[2],atom.Image_Flags[0],atom.Image_Flags[1],atom.Image_Flags[2]))
        f.write("\nBonds\n\n")
        for bond in self.Bond_List:
            f.write("%d %d %d %d\n" % (bond.Bond_ID,bond.LAMMPS_Type,bond.Bond_Master.Atom_ID,bond.Bond_Slave.Atom_ID))
        f.write("\nAngles\n\n")
        for ang in self.Angle_List:
            f.write("%d %d %d %d %d\n" % (ang.Angle_ID,ang.LAMMPS_Type,ang.Angle_Slave1.Atom_ID,ang.Angle_Master.Atom_ID,ang.Angle_Slave2.Atom_ID))
        f.write("\nDihedrals\n\n")
        for dih in self.Dihedral_List:
            f.write("%d %d %d %d %d %d\n" % (dih.Dihedral_ID,int(dih.LAMMPS_Type),dih.Dihedral_Slave1.Atom_ID,dih.Dihedral_Master1.Atom_ID,dih.Dihedral_Master2.Atom_ID,dih.Dihedral_Slave2.Atom_ID))
        f.write("\nImpropers\n\n")
        for imp in self.Improper_List:
            f.write("%d %d %d %d %d %d\n" % (imp.Improper_ID,int(imp.LAMMPS_Type),imp.Improper_Master.Atom_ID,imp.Improper_Slave1.Atom_ID,imp.Improper_Slave2.Atom_ID,imp.Improper_Slave3.Atom_ID))
        f.close()


    def Add_Bond_List(self):
        for atom in self.Atom_List:
            for bond in self.Bond_List:
                if bond.Bond_Master == atom:
                    atom.Bond_List.append(bond.Bond_Slave)
                elif bond.Bond_Slave == atom:
                    atom.Bond_List.append(bond.Bond_Master)
                if len(atom.Bond_List) == 4 or (atom.Element == "H" and len(atom.Bond_List) == 1):
                    break

    def Get_Atom(self,Atom_ID):
        ID = -1
        for i,atom in enumerate(self.Atom_List):
            if atom.Atom_ID == Atom_ID:
                ID = i
                break
        return self.Atom_List[i]
    
    def Get_Bond(self,Bond_ID):
        ID = -1
        for i,bond in enumerate(self.Bond_List):
            if bond.Bond_ID == Bond_ID:
                ID = i
                break
        return self.Bond_List[i]

    def Unwrap_Coords(self):
        for i,Atom in enumerate(self.Atom_List):
            Atom.Position[0] = Atom.Position[0] - self.I_Flags_List[i][0] * self.Box_Size
            Atom.Position[1] = Atom.Position[1] - self.I_Flags_List[i][1] * self.Box_Size
            Atom.Position[2] = Atom.Position[2] - self.I_Flags_List[i][2] * self.Box_Size

    def Adjust_COM(self):
        # This adjusts the center of mass and gives the molecule a random orientation
        x = random.random()*2*3.1415
        y = random.random()*2*3.1415
        
        C1 = np.cos(x)
        S1 = np.sin(x)
        C2 = np.cos(y)
        S2 = np.cos(y)
        
        for Atom_Obj in self.Atom_List:
            # First rotation
            xt = Atom_Obj.Position[0]
            yt = Atom_Obj.Position[1]
            Atom_Obj.Position[0] = xt*C1 - yt*S1
            Atom_Obj.Position[1] = xt*S1 + yt*C1
            # Second rotation
            #xt = Atom_Obj.Position[0]
            #zt = Atom_Obj.Position[2]
            #Atom_Obj.Position[0] = xt*C2 - zt*S2
            #Atom_Obj.Position[2] = xt*S2 + zt*C2
            
            Atom_Obj.Position += self.COM
        return

    def Separate_Sidechains(self,Sidechain_Atom_List):
        self.Conj_Atom_List = []
        self.Sidechain_Atom_List = []
        for atom in self.Atom_List:
            if atom.OPLS_Type not in Sidechain_Atom_List:
                self.Conj_Atom_List.append(atom)
            else:
                self.Sidechain_Atom_List.append(atom)

    def Sidechain_Grouping(self):
        self.Sidechain_Groups = []
        while len(self.Sidechain_Atom_List) > 0:
            addflag = True
            Group = [self.Sidechain_Atom_List[0]]
            del self.Sidechain_Atom_List[0]
            while addflag:
                addflag = False
                for atom in Group:
                    for bond in self.Bond_List:
                        if bond.Bond_Master.Atom_ID == atom.Atom_ID:
                            if bond.Bond_Slave in self.Sidechain_Atom_List:
                                Group.append(bond.Bond_Slave)
                                addflag = True
                                self.Sidechain_Atom_List.remove(bond.Bond_Slave)
                        if bond.Bond_Slave.Atom_ID == atom.Atom_ID:
                            if bond.Bond_Master in self.Sidechain_Atom_List:
                                Group.append(bond.Bond_Master)
                                addflag = True
                                self.Sidechain_Atom_List.remove(bond.Bond_Master)
            self.Sidechain_Groups.append(Group)



    def Dihedral_Angles(self,Conj_Dih_Types):
        self.Conj_Dihs = []
        Conj_Dih_Angles = []
        Conj_Dih_Energies = []
        Conj_Dih_Atoms = []
        for dih_type in Conj_Dih_Types:
            for dih in self.Dihedral_List:
                if dih.Style == dih_type:
                    self.Conj_Dihs.append(dih)

        for dih in self.Conj_Dihs:
            Vector1 = dih.Dihedral_Master1.Position - dih.Dihedral_Master2.Position
            Vector1 = Vector1/np.linalg.norm(Vector1)
            Vector2 = dih.Dihedral_Slave1.Position - dih.Dihedral_Master1.Position
            Vector2 = Vector2/np.linalg.norm(Vector2)
            Vector3 = dih.Dihedral_Master2.Position - dih.Dihedral_Slave2.Position
            Vector3 = Vector3/np.linalg.norm(Vector3)
            Plane1 = np.cross(Vector1,Vector2)
            Plane1 = Plane1/np.linalg.norm(Plane1)
            Plane2 = np.cross(Vector1,Vector3)
            Plane2 = Plane2/np.linalg.norm(Plane2)
            Theta = math.acos(np.dot(Plane1,Plane2))
            Conj_Dih_Atoms.append([dih.Dihedral_Slave1.Atom_ID,dih.Dihedral_Master1.Atom_ID,dih.Dihedral_Master2.Atom_ID,dih.Dihedral_Slave2.Atom_ID,Theta*180/math.pi])
            Conj_Dih_Angles.append(Theta*180/math.pi)
            E = .5 * dih.Coeffs[0] * (1 + math.cos(Theta)) + .5 * dih.Coeffs[1] * (1 - math.cos(2*Theta)) + .5 * dih.Coeffs[2] * (1 + math.cos(3*Theta)) + .5 * dih.Coeffs[3] * (1 -- math.cos(4*Theta))
            Conj_Dih_Energies.append(E)
        return np.sum(Conj_Dih_Energies)


    def Sidechain_Gyration_Tensor(self):
        SC_Kappas = []
        SC_Rgs = []
        for sidechain in self.Sidechain_Groups:
            SC_COM = np.zeros(3,dtype = float)
            SC_MW = 0.0
            for atom in sidechain:
                SC_COM += atom.Position * atom.Mass
                SC_MW += atom.Mass
            SC_COM = SC_COM/SC_MW
            Tensor = np.array([[0,0,0],[0,0,0],[0,0,0]],dtype = float)
            for atom in sidechain:
                for i in range(3):
                    for j in range(3):
                        Temp_Position = atom.Position - SC_COM
                        Tensor[i][j] += Temp_Position[i] * Temp_Position[j]
            Tensor = Tensor/len(self.Atom_List)
            eigval,eigvec = np.linalg.eig(Tensor)
            SC_Rg = (eigval[0]+eigval[1]+eigval[2])**0.5
            SC_kappa = 1.5*(eigval[0]**2+eigval[1]**2+eigval[2]**2)/((eigval[0]+eigval[1]+eigval[2])**2)-0.5
            SC_Rgs.append(SC_Rg)
            SC_Kappas.append(SC_kappa)
        return (np.sum(SC_Rgs)/len(SC_Rgs)),(np.sum(SC_Kappas)/len(SC_Kappas))


    def Center_Chain(self):
        Weighted_Position_Sum = np.array([0.0,0.0,0.0],dtype = float)
        Chain_MW = 0.0
        for atom in self.Conj_Atom_List:
            Weighted_Position_Sum += atom.Position*atom.Mass
            Chain_MW += atom.Mass
        Chain_Center = Weighted_Position_Sum/Chain_MW
        for atom in self.Atom_List:
            atom.Position -= Chain_Center

    def Total_Gyration_Tensor(self):
        Tensor = np.array([[0,0,0],[0,0,0],[0,0,0]],dtype = float)

        for atom in self.Atom_List:
            for i in range(3):
                for j in range(3):
                    Tensor[i][j] += atom.Position[i] * atom.Position[j]
        Tensor = Tensor/len(self.Atom_List)
        eigval,eigvec = np.linalg.eig(Tensor)
        self.Rg = (eigval[0]+eigval[1]+eigval[2])**0.5
        kappa = 1.5*(eigval[0]**2+eigval[1]**2+eigval[2]**2)/((eigval[0]+eigval[1]+eigval[2])**2)-0.5
        return self.Rg,kappa

    def Main_Chain_Gyration_Tensor(self):
        Tensor = np.array([[0,0,0],[0,0,0],[0,0,0]],dtype = float)
        for atom in self.Conj_Atom_List:
            for i in range(3):
                for j in range(3):
                    Tensor[i][j] += atom.Position[i] * atom.Position[j]
        Tensor = Tensor/len(self.Conj_Atom_List)
        eigval,eigvec = np.linalg.eig(Tensor)
        Chain_Rg = (eigval[0]+eigval[1]+eigval[2])**0.5
        kappa = 1.5*(eigval[0]**2+eigval[1]**2+eigval[2]**2)/((eigval[0]+eigval[1]+eigval[2])**2)-0.5
        return Chain_Rg,kappa

    def Create_Backbone_Vectors(self,List1,List2):
        self.Backbone_Vectors = []
        for i in range(len(List1)):
            temp_vec = self.Get_Atom(List1[i]).Position - self.Get_Atom(List2[i]).Position
            temp_vec = temp_vec/np.linalg.norm(temp_vec)
            self.Backbone_Vectors.append(temp_vec)
            if i != (len(List1)-1):
                temp_vec = self.Get_Atom(List2[i]).Position - self.Get_Atom(List1[i+1]).Position
                temp_vec = temp_vec/np.linalg.norm(temp_vec)
                self.Backbone_Vectors.append(temp_vec)

    def Calculate_QTensor(self):
        Tensor = np.array([[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]],dtype = float)
        for vec in self.Backbone_Vectors:
            for i in range(3):
                for j in range(3):
                    Tensor[i][j] = Tensor[i][j] + vec[i] * vec[j]
        Tensor = Tensor/len(self.Backbone_Vectors)
        for i in range(3):
            Tensor[i][i] = Tensor[i][i] - (.333333333333)
        eigval,eigvec = np.linalg.eig(Tensor)
        eigval = eigval *1.5
        eigval = sorted(eigval,reverse = True)
        return eigval

    def Calculate_Coord_Number_ChainOnly(self,r0,Monomer_Length):
        coord = 0
        for atom1 in self.Conj_Atom_List:
            for atom2 in self.Conj_Atom_List:
                if abs(atom1.Atom_ID - atom2.Atom_ID) > Monomer_Length:
                    rij = np.linalg.norm(atom1.Position - atom2.Position)
                    coord += (1-(rij/r0)**6)/(1-(rij/r0)**12)
        return coord

    def Calculate_Coord_Number_SidechainsOnly(self,r0,Monomer_Length):
        coord = 0
        for atomlist1 in self.Sidechain_Groups:
            for atom1 in atomlist1:
                for atomlist2 in self.Sidechain_Groups:
                    for atom2 in atomlist2:
                        if abs(atom1.Atom_ID - atom2.Atom_ID) > Monomer_Length:
                            rij = np.linalg.norm(atom1.Position - atom2.Position)
                            coord += (1-(rij/r0)**6)/(1-(rij/r0)**12)
        return coord

    def Calculate_Coord_Number_BetweenChainSidechain(self,r0,Monomer_Length):
        coord = 0
        for atom1 in self.Conj_Atom_List:
            for atomlist in self.Sidechain_Groups:
                for atom2 in atomlist:
                    if abs(atom1.Atom_ID - atom2.Atom_ID) > Monomer_Length:
                        rij = np.linalg.norm(atom1.Position - atom2.Position)
                        coord += (1-(rij/r0)**6)/(1-(rij/r0)**12)
        return coord

    def Write_MetaD_File(self):
        Atom_String = ""
        In_Temp = Configure.Template_Path + "plumed.dat"
        In_File = "MetaD_plumed.dat"
        for i,atom in enumerate(self.Conj_Atom_List):
            if i==0:
                Atom_String = Atom_String + str(atom.Atom_ID) 
            else:
                Atom_String = Atom_String + "," + str(atom.Atom_ID) 
        with open(In_Temp) as f:
            template = f.read()
        s = template.format(conj=Atom_String)
        with open(In_File,'w') as f:
            f.write(s)

    def Write_MetaD_Test_File(self,Deg_Polym,Ring_Diff,Period,e2e_atom,Planes):
        Conj_Atom_String = ""
        Conj_Dih_List = []
        Conj_Dih_String = ""
        Sidechain_String = ""
        In_File = "MetaD_plumed_Test.dat"

        for i,atom in enumerate(self.Conj_Atom_List):
            if i==0:
                Conj_Atom_String = Conj_Atom_String + str(atom.Atom_ID) 
            else:
                Conj_Atom_String = Conj_Atom_String + "," + str(atom.Atom_ID) 

        for i,atom in enumerate(self.Sidechain_Atom_List):
            if i==0:
                Sidechain_String = Sidechain_String + str(atom.Atom_ID) 
            else:
                Sidechain_String = Sidechain_String + "," + str(atom.Atom_ID) 

        for dih in self.Conj_Dihs:
            Conj_Dih_List.append([dih.Dihedral_Slave1.Atom_ID,dih.Dihedral_Master1.Atom_ID,dih.Dihedral_Master2.Atom_ID,dih.Dihedral_Slave2.Atom_ID])
        Conj_Dih_List = sorted(Conj_Dih_List,key=itemgetter(0))
        for i,dih1,dih2 in zip(range(len(Conj_Dih_List)),Conj_Dih_List[0:-2],Conj_Dih_List[1:-1]):
            Conj_Dih_String = Conj_Dih_String + "\tATOMS%d=%d,%d,%d,%d,%d,%d,%d,%d\n" % (i+1,dih1[0],dih1[1],dih1[2],dih1[3],dih2[0],dih2[1],dih2[2],dih2[3])

        with open(In_File,'w') as f:
            f.write("GYRATION TYPE=RADIUS ATOMS=%s LABEL=conj_rg\n\n" % Conj_Atom_String)
            f.write("GYRATION TYPE=KAPPA2 ATOMS=%s LABEL=conj_anis\n\n" % Conj_Atom_String)
            f.write("GYRATION TYPE=ASPHERICITY ATOMS=%s LABEL=conj_asphere\n\n" % Conj_Atom_String)
            f.write("GYRATION TYPE=ACYLINDRICITY ATOMS=%s LABEL=conj_acylind\n\n" % Conj_Atom_String)
            f.write("GYRATION TYPE=RGYR_1 ATOMS=%s LABEL=conj_rg1\n\n" % Conj_Atom_String)
            f.write("GYRATION TYPE=RGYR_2 ATOMS=%s LABEL=conj_rg2\n\n" % Conj_Atom_String)
            f.write("GYRATION TYPE=RGYR_3 ATOMS=%s LABEL=conj_rg3\n\n" % Conj_Atom_String)
            f.write("Q6 SPECIES=%s R_0=5.0 MEAN VMEAN LABEL=q6\n\n" % Conj_Atom_String)
            f.write("Q3 SPECIES=%s R_0=5.0 MEAN VMEAN LABEL=q3\n\n" % Conj_Atom_String)
            f.write("Q4 SPECIES=%s R_0=5.0 MEAN VMEAN LABEL=q4\n\n" % Conj_Atom_String)
            f.write("DIHCOR ...\n%s\tLABEL=dih\n... DIHCOR\n\n" % Conj_Dih_String)
            f.write("COORDINATION GROUPA=%s R_0=0.5 NLIST NL_CUTOFF=10.0 NL_STRIDE=100 LABEL=main_coord\n\n" % Conj_Atom_String)
            f.write("COORDINATION GROUPA=%s R_0=0.5 NLIST NL_CUTOFF=10.0 NL_STRIDE=100 LABEL=side_coord\n\n" % Sidechain_String)
            f.write("PLANES ")
            p = 1
            for deg_poly in range(Deg_Polym):
                for Ring_List in Planes:
                    f.write("MOL%d=%d,%d,%d " % (p,Ring_List[0]+deg_poly*Period,Ring_List[1]+deg_poly*Period,Ring_List[2]+deg_poly*Period))
                    p+=1
            f.write(" LABEL=plane_circ\n\n")
            f.write("LOCAL_AVERAGE SPECIES=plane_circ R_0=1.0 VSUM LABEL=ave_plane_circ\n\n")
            f.write("SMAC ...\n\tSPECIES=plane_circ\n\tR_0=2.50\n\tKERNEL1={GAUSSIAN CENTER=0 SIGMA=0.175}\n\tKERNEL2={GAUSSIAN CENTER=0 SIGMA=0.175}\n\tSWITCH_COORD={EXP R_0=2.50}\n\tMEAN\n\tLABEL=smac_test_1\n... SMAC\n\n")
            f.write("SMAC ...\n\tSPECIES=plane_circ\n\tR_0=2.50\n\tKERNEL1={GAUSSIAN CENTER=0 SIGMA=0.175}\n\tKERNEL2={GAUSSIAN CENTER=0 SIGMA=0.35}\n\tSWITCH_COORD={EXP R_0=2.50}\n\tMEAN\n\tLABEL=smac_test_2\n... SMAC\n\n")
            f.write("SMAC ...\n\tSPECIES=plane_circ\n\tR_0=2.50\n\tKERNEL1={GAUSSIAN CENTER=0 SIGMA=0.175}\n\tKERNEL2={GAUSSIAN CENTER=0 SIGMA=0.35}\n\tKERNEL3={GAUSSIAN CENTER=0 SIGMA=0.70}\n\tSWITCH_COORD={EXP R_0=2.50}\n\tMEAN\n\tLABEL=smac_test_3\n... SMAC\n\n")
            f.write("SMAC ...\n\tSPECIES=plane_circ\n\tR_0=2.50\n\tKERNEL1={GAUSSIAN CENTER=0 SIGMA=0.175}\n\tKERNEL2={GAUSSIAN CENTER=0 SIGMA=0.35}\n\tKERNEL3={GAUSSIAN CENTER=0 SIGMA=0.70}\n\tKERNEL4={GAUSSIAN CENTER=0 SIGMA=1.4}\n\tSWITCH_COORD={EXP R_0=2.50}\n\tMEAN\n\tLABEL=smac_test_4\n... SMAC\n\n")
            f.write("SMAC ...\n\tSPECIES=plane_circ\n\tR_0=2.50\n\tKERNEL1={GAUSSIAN CENTER=pi SIGMA=0.175}\n\tSWITCH_COORD={EXP R_0=2.50}\n\tMEAN\n\tLABEL=smac_test_5\n... SMAC\n\n")
            f.write("SMAC ...\n\tSPECIES=plane_circ\n\tR_0=2.50\n\tKERNEL1={GAUSSIAN CENTER=pi SIGMA=0.175}\n\tKERNEL2={GAUSSIAN CENTER=pi SIGMA=0.35}\n\tSWITCH_COORD={EXP R_0=2.50}\n\tMEAN\n\tLABEL=smac_test_6\n... SMAC\n\n")
            f.write("SMAC ...\n\tSPECIES=plane_circ\n\tR_0=2.50\n\tKERNEL1={GAUSSIAN CENTER=pi SIGMA=0.175}\n\tKERNEL2={GAUSSIAN CENTER=pi SIGMA=0.35}\n\tKERNEL3={GAUSSIAN CENTER=pi SIGMA=0.70}\n\tSWITCH_COORD={EXP R_0=2.50}\n\tMEAN\n\tLABEL=smac_test_7\n... SMAC\n\n")
            f.write("SMAC ...\n\tSPECIES=plane_circ\n\tR_0=2.50\n\tKERNEL1={GAUSSIAN CENTER=pi SIGMA=0.175}\n\tKERNEL2={GAUSSIAN CENTER=pi SIGMA=0.35}\n\tKERNEL3={GAUSSIAN CENTER=pi SIGMA=0.70}\n\tKERNEL4={GAUSSIAN CENTER=pi SIGMA=1.4}\n\tSWITCH_COORD={EXP R_0=2.50}\n\tMEAN\n\tLABEL=smac_test_8\n... SMAC\n\n")
            f.write("SMAC ...\n\tSPECIES=plane_circ\n\tR_0=2.50\n\tKERNEL1={GAUSSIAN CENTER=pi SIGMA=0.175}\n\tKERNEL2={GAUSSIAN CENTER=pi SIGMA=0.35}\n\tKERNEL3={GAUSSIAN CENTER=pi SIGMA=0.70}\n\tKERNEL4={GAUSSIAN CENTER=pi SIGMA=1.4}\n\tKERNEL5={GAUSSIAN CENTER=0 SIGMA=0.175}\n\tKERNEL6={GAUSSIAN CENTER=0 SIGMA=0.35}\n\tKERNEL7={GAUSSIAN CENTER=0 SIGMA=0.70}\n\tKERNEL8={GAUSSIAN CENTER=0 SIGMA=1.4}\n\tSWITCH_COORD={EXP R_0=2.50}\n\tMEAN\n\tLABEL=smac_test_9\n... SMAC\n\n")
            f.write("SMAC ...\n\tSPECIES=plane_circ\n\tR_0=1.0\n\tKERNEL1={GAUSSIAN CENTER=0 SIGMA=0.175}\n\tKERNEL2={GAUSSIAN CENTER=0 SIGMA=0.175}\n\tSWITCH_COORD={EXP R_0=1.0}\n\tMEAN\n\tLABEL=smac_test_10\n... SMAC\n\n")
            f.write("SMAC ...\n\tSPECIES=plane_circ\n\tR_0=1.0\n\tKERNEL1={GAUSSIAN CENTER=0 SIGMA=0.175}\n\tKERNEL2={GAUSSIAN CENTER=0 SIGMA=0.35}\n\tSWITCH_COORD={EXP R_0=1.0}\n\tMEAN\n\tLABEL=smac_test_11\n... SMAC\n\n")
            f.write("SMAC ...\n\tSPECIES=plane_circ\n\tR_0=1.0\n\tKERNEL1={GAUSSIAN CENTER=0 SIGMA=0.175}\n\tKERNEL2={GAUSSIAN CENTER=0 SIGMA=0.35}\n\tKERNEL3={GAUSSIAN CENTER=0 SIGMA=0.70}\n\tSWITCH_COORD={EXP R_0=1.0}\n\tMEAN\n\tLABEL=smac_test_12\n... SMAC\n\n")
            f.write("SMAC ...\n\tSPECIES=plane_circ\n\tR_0=1.0\n\tKERNEL1={GAUSSIAN CENTER=0 SIGMA=0.175}\n\tKERNEL2={GAUSSIAN CENTER=0 SIGMA=0.35}\n\tKERNEL3={GAUSSIAN CENTER=0 SIGMA=0.70}\n\tKERNEL4={GAUSSIAN CENTER=0 SIGMA=1.4}\n\tSWITCH_COORD={EXP R_0=1.0}\n\tMEAN\n\tLABEL=smac_test_13\n... SMAC\n\n")
            f.write("SMAC ...\n\tSPECIES=plane_circ\n\tR_0=1.0\n\tKERNEL1={GAUSSIAN CENTER=pi SIGMA=0.175}\n\tSWITCH_COORD={EXP R_0=1.0}\n\tMEAN\n\tLABEL=smac_test_14\n... SMAC\n\n")
            f.write("SMAC ...\n\tSPECIES=plane_circ\n\tR_0=1.0\n\tKERNEL1={GAUSSIAN CENTER=pi SIGMA=0.175}\n\tKERNEL2={GAUSSIAN CENTER=pi SIGMA=0.35}\n\tSWITCH_COORD={EXP R_0=1.0}\n\tMEAN\n\tLABEL=smac_test_15\n... SMAC\n\n")
            f.write("SMAC ...\n\tSPECIES=plane_circ\n\tR_0=1.0\n\tKERNEL1={GAUSSIAN CENTER=pi SIGMA=0.175}\n\tKERNEL2={GAUSSIAN CENTER=pi SIGMA=0.35}\n\tKERNEL3={GAUSSIAN CENTER=pi SIGMA=0.70}\n\tSWITCH_COORD={EXP R_0=1.0}\n\tMEAN\n\tLABEL=smac_test_16\n... SMAC\n\n")
            f.write("SMAC ...\n\tSPECIES=plane_circ\n\tR_0=1.0\n\tKERNEL1={GAUSSIAN CENTER=pi SIGMA=0.175}\n\tKERNEL2={GAUSSIAN CENTER=pi SIGMA=0.35}\n\tKERNEL3={GAUSSIAN CENTER=pi SIGMA=0.70}\n\tKERNEL4={GAUSSIAN CENTER=pi SIGMA=1.4}\n\tSWITCH_COORD={EXP R_0=1.0}\n\tMEAN\n\tLABEL=smac_test_17\n... SMAC\n\n")
            f.write("SMAC ...\n\tSPECIES=plane_circ\n\tR_0=1.0\n\tKERNEL1={GAUSSIAN CENTER=pi SIGMA=0.175}\n\tKERNEL2={GAUSSIAN CENTER=pi SIGMA=0.35}\n\tKERNEL3={GAUSSIAN CENTER=pi SIGMA=0.70}\n\tKERNEL4={GAUSSIAN CENTER=pi SIGMA=1.4}\n\tKERNEL5={GAUSSIAN CENTER=0 SIGMA=0.175}\n\tKERNEL6={GAUSSIAN CENTER=0 SIGMA=0.35}\n\tKERNEL7={GAUSSIAN CENTER=0 SIGMA=0.70}\n\tKERNEL8={GAUSSIAN CENTER=0 SIGMA=1.4}\n\tSWITCH_COORD={EXP R_0=1.0}\n\tMEAN\n\tLABEL=smac_test_18\n... SMAC\n\n")
            f.write("CUSTOM LABEL=smac_div_1 ARG=smac_test_3.mean,smac_test_7.mean FUNC=x/y PERIODIC=NO\n\n")
            f.write("CUSTOM LABEL=smac_div_2 ARG=smac_test_12.mean,smac_test_16.mean FUNC=x/y PERIODIC=NO\n\n")
            f.write("CUSTOM LABEL=smac_div_3 ARG=smac_test_3.mean,smac_test_9.mean FUNC=x/y PERIODIC=NO\n\n")
            f.write("CUSTOM LABEL=smac_div_4 ARG=smac_test_7.mean,smac_test_9.mean FUNC=x/y PERIODIC=NO\n\n")
            f.write("CUSTOM LABEL=smac_div_5 ARG=smac_test_12.mean,smac_test_18.mean FUNC=x/y PERIODIC=NO\n\n")
            f.write("CUSTOM LABEL=smac_div_6 ARG=smac_test_16.mean,smac_test_18.mean FUNC=x/y PERIODIC=NO\n\n")
            j = 1
            Tors_String = ""
            for dih in Conj_Dih_List:
                f.write("TORSION ATOMS=%d,%d,%d,%d LABEL=t%d\n\n" % (dih[0],dih[1],dih[2],dih[3],j))
                Tors_String = Tors_String + "t%d," % j
                j+=1
            j = 1
            Tors_String.strip(',')
            Torsion_Helix_Groups = []
            for i in range(len(Conj_Dih_List)-5):
                Torsion_Helix_Groups.append(Conj_Dih_List[i:i+6])
                f.write("circ%d: MATHEVAL ARG=t%d,t%d,t%d,t%d,t%d,t%d VAR=a%d,b%d,c%d,d%d,e%d,f%d FUNC=1-(1-((cos(a%d)+cos(b%d)+cos(c%d)+cos(d%d)+cos(e%d)+cos(f%d))/3)^6)/(1-((cos(a%d)+cos(b%d)+cos(c%d)+cos(d%d)+cos(e%d)+cos(f%d))/3)^12) PERIODIC=NO\n\n" % (i+1,i+1,i+2,i+3,i+4,i+5,i+6,i+1,i+1,i+1,i+1,i+1,i+1,i+1,i+1,i+1,i+1,i+1,i+1,i+1,i+1,i+1,i+1,i+1,i+1))
            f.write("circ: COMBINE ARG=circ1")
            for i in range(2,len(Conj_Dih_List)-4):
                f.write(",circ%d" % i)
            f.write(" PERIODIC=NO\n\n")
            f.write("ene: ENERGY\n\n")

            q = 1
            r = 1
            for i in range(Deg_Polym):
                f.write("d%d: DISTANCE ATOMS=%d,%d COMPONENTS\n\n" % (q,self.Conj_Atom_List[0].Atom_ID+i*Period,self.Conj_Atom_List[1].Atom_ID+i*Period))
                f.write("d%d_norm: DISTANCE ATOMS=%d,%d\n\n" % (q,self.Conj_Atom_List[0].Atom_ID+i*Period,self.Conj_Atom_List[1].Atom_ID+i*Period))
                f.write("d%d_x_norm: MATHEVAL ARG=d%d.x,d%d_norm FUNC=x/y PERIODIC=NO\n\n" % (q,q,q))
                f.write("d%d_y_norm: MATHEVAL ARG=d%d.y,d%d_norm FUNC=x/y PERIODIC=NO\n\n" % (q,q,q))
                f.write("d%d_z_norm: MATHEVAL ARG=d%d.z,d%d_norm FUNC=x/y PERIODIC=NO\n\n" % (q,q,q))
                q+=1
                f.write("d%d: DISTANCE ATOMS=%d,%d COMPONENTS\n\n" % (q,self.Conj_Atom_List[1].Atom_ID+i*Period,self.Conj_Atom_List[2].Atom_ID+i*Period))
                f.write("d%d_norm: DISTANCE ATOMS=%d,%d\n\n" % (q,self.Conj_Atom_List[1].Atom_ID+i*Period,self.Conj_Atom_List[2].Atom_ID+i*Period))
                f.write("d%d_x_norm: MATHEVAL ARG=d%d.x,d%d_norm FUNC=x/y PERIODIC=NO\n\n" % (q,q,q))
                f.write("d%d_y_norm: MATHEVAL ARG=d%d.y,d%d_norm FUNC=x/y PERIODIC=NO\n\n" % (q,q,q))
                f.write("d%d_z_norm: MATHEVAL ARG=d%d.z,d%d_norm FUNC=x/y PERIODIC=NO\n\n" % (q,q,q))
                q+=1
                f.write("cross_x_%d: MATHEVAL ARG=d%d_y_norm,d%d_z_norm,d%d_y_norm,d%d_z_norm VAR=a,b,c,d FUNC=a*d-b*c PERIODIC=NO\n\n" % (r,q-2,q-2,q-1,q-1))
                f.write("cross_y_%d: MATHEVAL ARG=d%d_x_norm,d%d_z_norm,d%d_x_norm,d%d_z_norm VAR=a,b,c,d FUNC=b*c-a*d PERIODIC=NO\n\n" % (r,q-2,q-2,q-1,q-1))
                f.write("cross_z_%d: MATHEVAL ARG=d%d_x_norm,d%d_y_norm,d%d_x_norm,d%d_y_norm VAR=a,b,c,d FUNC=a*d-b*c PERIODIC=NO\n\n" % (r,q-2,q-2,q-1,q-1))
                f.write("cross_%d_norm: MATHEVAL ARG=cross_x_%d,cross_y_%d,cross_z_%d FUNC=sqrt(x*x+y*y+z*z) PERIODIC=NO\n\n" % (r,r,r,r))
                f.write("cross_x_%d_norm: MATHEVAL ARG=cross_x_%d,cross_%d_norm FUNC=x/y PERIODIC=NO\n\n" % (r,r,r))
                f.write("cross_y_%d_norm: MATHEVAL ARG=cross_y_%d,cross_%d_norm FUNC=x/y PERIODIC=NO\n\n" % (r,r,r))
                f.write("cross_z_%d_norm: MATHEVAL ARG=cross_z_%d,cross_%d_norm FUNC=x/y PERIODIC=NO\n\n" % (r,r,r))
                r+=1
                f.write("d%d: DISTANCE ATOMS=%d,%d COMPONENTS\n\n" % (q,self.Conj_Atom_List[0].Atom_ID+i*Period+Ring_Diff,self.Conj_Atom_List[1].Atom_ID+i*Period+Ring_Diff))
                f.write("d%d_norm: DISTANCE ATOMS=%d,%d\n\n" % (q,self.Conj_Atom_List[0].Atom_ID+i*Period+Ring_Diff,self.Conj_Atom_List[1].Atom_ID+i*Period+Ring_Diff))
                f.write("d%d_x_norm: MATHEVAL ARG=d%d.x,d%d_norm FUNC=x/y PERIODIC=NO\n\n" % (q,q,q))
                f.write("d%d_y_norm: MATHEVAL ARG=d%d.y,d%d_norm FUNC=x/y PERIODIC=NO\n\n" % (q,q,q))
                f.write("d%d_z_norm: MATHEVAL ARG=d%d.z,d%d_norm FUNC=x/y PERIODIC=NO\n\n" % (q,q,q))
                q+=1
                f.write("d%d: DISTANCE ATOMS=%d,%d COMPONENTS\n\n" % (q,self.Conj_Atom_List[1].Atom_ID+i*Period+Ring_Diff,self.Conj_Atom_List[2].Atom_ID+i*Period+Ring_Diff))
                f.write("d%d_norm: DISTANCE ATOMS=%d,%d\n\n" % (q,self.Conj_Atom_List[1].Atom_ID+i*Period+Ring_Diff,self.Conj_Atom_List[2].Atom_ID+i*Period+Ring_Diff))
                f.write("d%d_x_norm: MATHEVAL ARG=d%d.x,d%d_norm FUNC=x/y PERIODIC=NO\n\n" % (q,q,q))
                f.write("d%d_y_norm: MATHEVAL ARG=d%d.y,d%d_norm FUNC=x/y PERIODIC=NO\n\n" % (q,q,q))
                f.write("d%d_z_norm: MATHEVAL ARG=d%d.z,d%d_norm FUNC=x/y PERIODIC=NO\n\n" % (q,q,q))
                q+=1
                f.write("cross_x_%d: MATHEVAL ARG=d%d_y_norm,d%d_z_norm,d%d_y_norm,d%d_z_norm VAR=a,b,c,d FUNC=a*d-b*c PERIODIC=NO\n\n" % (r,q-2,q-2,q-1,q-1))
                f.write("cross_y_%d: MATHEVAL ARG=d%d_x_norm,d%d_z_norm,d%d_x_norm,d%d_z_norm VAR=a,b,c,d FUNC=b*c-a*d PERIODIC=NO\n\n" % (r,q-2,q-2,q-1,q-1))
                f.write("cross_z_%d: MATHEVAL ARG=d%d_x_norm,d%d_y_norm,d%d_x_norm,d%d_y_norm VAR=a,b,c,d FUNC=a*d-b*c PERIODIC=NO\n\n" % (r,q-2,q-2,q-1,q-1))
                f.write("cross_%d_norm: MATHEVAL ARG=cross_x_%d,cross_y_%d,cross_z_%d FUNC=sqrt(x*x+y*y+z*z) PERIODIC=NO\n\n" % (r,r,r,r))
                f.write("cross_x_%d_norm: MATHEVAL ARG=cross_x_%d,cross_%d_norm FUNC=x/y PERIODIC=NO\n\n" % (r,r,r))
                f.write("cross_y_%d_norm: MATHEVAL ARG=cross_y_%d,cross_%d_norm FUNC=x/y PERIODIC=NO\n\n" % (r,r,r))
                f.write("cross_z_%d_norm: MATHEVAL ARG=cross_z_%d,cross_%d_norm FUNC=x/y PERIODIC=NO\n\n" % (r,r,r))
                r+=1
            for cross_1,cross_2,cross_3,cross_4,cross_5,cross_6 in zip(range(Deg_Polym*2-5),range(1,Deg_Polym*2-4),range(2,Deg_Polym*2-3),range(3,Deg_Polym*2-2),range(4,Deg_Polym*2-1),range(5,Deg_Polym*2)):
                f.write("curve_vec_norm_%d: MATHEVAL ARG=cross_x_%d_norm,cross_y_%d_norm,cross_z_%d_norm,cross_x_%d_norm,cross_y_%d_norm,cross_z_%d_norm,cross_x_%d_norm,cross_y_%d_norm,cross_z_%d_norm,cross_x_%d_norm,cross_y_%d_norm,cross_z_%d_norm,cross_x_%d_norm,cross_y_%d_norm,cross_z_%d_norm,cross_x_%d_norm,cross_y_%d_norm,cross_z_%d_norm VAR=c_x_1,c_y_1,c_z_1,c_x_2,c_y_2,c_z_2,c_x_3,c_y_3,c_z_3,c_x_4,c_y_4,c_z_4,c_x_5,c_y_5,c_z_5,c_x_6,c_y_6,c_z_6 FUNC=sqrt((c_x_1+c_x_2+c_x_3+c_x_4+c_x_5+c_x_6)^2+(c_y_1+c_y_2+c_y_3+c_y_4+c_y_5+c_y_6)^2+(c_z_1+c_z_2+c_z_3+c_z_4+c_z_5+c_z_6)^2) PERIODIC=NO\n\n" % (cross_1+1,cross_1+1,cross_1+1,cross_1+1,cross_2+1,cross_2+1,cross_2+1,cross_3+1,cross_3+1,cross_3+1,cross_4+1,cross_4+1,cross_4+1,cross_5+1,cross_5+1,cross_5+1,cross_6+1,cross_6+1,cross_6+1))
            f.write("vec_circuity: COMBINE ARG=")
            f.write("curve_vec_norm_1")
            for curve_norm in range(2,Deg_Polym*2-4):
                f.write(",curve_vec_norm_%d" % curve_norm)
            f.write(" PERIODIC=NO\n\n")
            r = 1
            for e2e1,e2e2 in zip(range(0,Deg_Polym-3),range(3,Deg_Polym)):
                f.write("e2e%d: DISTANCE ATOMS=%d,%d\n\n" % (r,e2e_atom + e2e1*Period,e2e_atom + e2e2*Period))
                r+=1
                f.write("e2e%d: DISTANCE ATOMS=%d,%d\n\n" % (r,e2e_atom + Ring_Diff + e2e1*Period,e2e_atom + Ring_Diff + e2e2*Period))
                r+=1
            f.write("e2e_plumed: DISTANCES ")
            r=1
            for e2e1,e2e2 in zip(range(0,Deg_Polym-3),range(3,Deg_Polym)):
                f.write("ATOMS%d=%d,%d " % (r,e2e_atom + e2e1*Period,e2e_atom + e2e2*Period))
                r+=1
                f.write("ATOMS%d=%d,%d " % (r,e2e_atom + Ring_Diff + e2e1*Period,e2e_atom + Ring_Diff + e2e2*Period))
                r+=1
            f.write("MEAN\n\n")
            f.write("e2e: COMBINE ARG=e2e1")
            for e2e3 in range(1,(Deg_Polym-2)*2-1):
                f.write(",e2e%d" % e2e3)
            f.write(" PERIODIC=NO\n\n")
            f.write("PRINT ARG=conj_rg,conj_anis,conj_asphere,conj_acylind,conj_rg1,conj_rg2,conj_rg3,q6.vmean,q6.mean,q3.vmean,q3.mean,q4.vmean,q4.mean,dih,main_coord,circ,ene,vec_circuity,e2e_plumed.mean,side_coord,smac_test_1.mean,smac_test_2.mean,smac_test_3.mean,smac_test_4.mean,smac_test_5.mean,smac_test_6.mean,smac_test_7.mean,smac_test_8.mean,smac_test_9.mean,smac_test_10.mean,smac_test_11.mean,smac_test_12.mean,smac_test_13.mean,smac_test_14.mean,smac_test_15.mean,smac_test_16.mean,smac_test_17.mean,smac_test_18.mean,smac_div_1,smac_div_2,smac_div_3,smac_div_4,smac_div_5,smac_div_6")
            j = 1
            for dih in Conj_Dih_List:
                f.write(",t%d" % j)
                j+=1
            f.write(" FILE=COLVAR_TEST STRIDE=1\n\nDUMPMULTICOLVAR DATA=plane_circ FILE=Planecirc.xyz")

    def Write_MetaD_Circuity(self):
        Conj_Atom_String = ""
        Conj_Dih_List = []
        Conj_Dih_String = ""
        Tors_String = ""
        In_File = "MetaD_plumed_Test.dat"

        for i,atom in enumerate(self.Conj_Atom_List):
            if i==0:
                Conj_Atom_String = Conj_Atom_String + str(atom.Atom_ID) 
            else:
                Conj_Atom_String = Conj_Atom_String + "," + str(atom.Atom_ID) 

        for dih in self.Conj_Dihs:
            Conj_Dih_List.append([dih.Dihedral_Slave1.Atom_ID,dih.Dihedral_Master1.Atom_ID,dih.Dihedral_Master2.Atom_ID,dih.Dihedral_Slave2.Atom_ID])
        Conj_Dih_List = sorted(Conj_Dih_List,key=itemgetter(0))
        for i,dih1,dih2 in zip(range(len(Conj_Dih_List)),Conj_Dih_List[0:-2],Conj_Dih_List[1:-1]):
            Conj_Dih_String = Conj_Dih_String + "\tATOMS%d=%d,%d,%d,%d,%d,%d,%d,%d\n" % (i+1,dih1[0],dih1[1],dih1[2],dih1[3],dih2[0],dih2[1],dih2[2],dih2[3])

        with open(In_File,'w') as f:
            f.write("COORDINATION GROUPA=%s R_0=0.5 NLIST NL_CUTOFF=10.0 NL_STRIDE=100 LABEL=main_coord\n\n" % Conj_Atom_String)
            j = 1
            for dih in Conj_Dih_List:
                f.write("TORSION ATOMS=%d,%d,%d,%d LABEL=t%d\n\n" % (dih[0],dih[1],dih[2],dih[3],j))
                Tors_String = Tors_String + "t%d," % j
                j+=1
            Tors_String.strip(',')
            Torsion_Helix_Groups = []
            for i in range(len(Conj_Dih_List)-6):
                Torsion_Helix_Groups.append(Conj_Dih_List[i:i+6])
                f.write("circ%d: MATHEVAL ARG=t%d,t%d,t%d,t%d,t%d,t%d FUNC=1-(1-(cos(t%d)+(cos(t%d)+(cos(t%d)+(cos(t%d)+(cos(t%d)+(cos(t%d)/3)^6)/(1-(cos(t%d)+(cos(t%d)+(cos(t%d)+(cos(t%d)+(cos(t%d)+(cos(t%d)/3)^12) PERIODIC=NO\n\n" % (i+1,i+1,i+2,i+3,i+4,i+5,i+6,i+1,i+2,i+3,i+4,i+5,i+6,i+1,i+2,i+3,i+4,i+5,i+6))

            f.write("circ: COMBINE ARG=circ1")
            for i in range(2,len(Conj_Dih_List)-5):
                f.write(",circ%d" % i)
            f.write(" PERIODIC=NO")


            f.write("PRINT ARG=main_coord,circ")
            j = 1
            for dih in Conj_Dih_List:
                f.write(",t%d" % j)
                j+=1
            f.write(" FILE=COLVAR_CIRC STRIDE=1")

    def Separate_Rings(self,Dividing_Atoms,Deg_Polym,Sep_Num):
        Rings = []
        Min = 0
        for mult in range(Deg_Polym):
            for div in Dividing_Atoms:
                Max = div + mult*Sep_Num
                Ring = []
                for atom in self.Conj_Atom_List:
                    if atom.Atom_ID > Min and atom.Atom_ID < Max:
                        Ring.append(atom)
                Rings.append(Ring)
                Min = Max
        return Rings











