#! usr/bin/python


# Import relevant modules
import numpy as np
# Atom Class

Mass_Dict = {"C":12.011, "H":1.008, "Cl":35.453, "Si":28.086, "O":15.9994, "S":32.06, "F":18.998, "N":14.007, "I":126.905}
Reverse_Dict = {12.011:"C", 1.008:"H", 35.453:"Cl", 28.086:"Si", 15.9994:"O", 32.06:"S", 18.998:"F", 14.007:"N", 126.905:"I"}

class Atom(object):
    """
    Class defining an atom
    instance variables:
        Position = Position Vector numpy[3,float]
        Charge = Float
        Element = String
        Bond_List = list of atom objects
        Atom_ID = 0
    """
    def __init__(self, Position, Element, Atom_ID):
        self.Position = Position # Numpy[3,float]
        self.Element = Element # String
        self.Charge = np.longdouble(0.0) # Float
        self.Bond_List = [] # List of atom objects
        self.OPLS_Type = 0
        self.OPLS_Class = 0
        self.Mass = Mass_Dict[self.Element]
        self.Sigma = 0.0
        self.Epsilon = 0.0
        self.Atom_ID = Atom_ID
        self.LAMMPS_Type = 0
        self.System_ID = 0
        #self.COM_Position = np.zeros(3,dtype=float)
        self.Image_Flags = np.zeros(3,dtype=int)
        self.Unwrapped_Position = np.zeros(3,dtype=float)
        self.ring= False # Assigned in Ring.py
        return

    def Assign_OPLS_ID(self, Fullerene = False):
    # Function for assigning OPLS Types and Classes based on connectivity of the atom
        if self.Element == "C":
            Temp_Bond_List = sorted([ Atomobj.Element for Atomobj in self.Bond_List ])
            if len(Temp_Bond_List) == 2:
                if Temp_Bond_List == ['C','N']:
                    #Acetonitrile C->C-N
                    self.OPLS_Type = 695
                    self.OPLS_Class = 19
            if len(Temp_Bond_List) == 3:
                if Temp_Bond_List == ['C', 'C', 'H']:
                    # Aromatic Carbon
                    self.OPLS_Type = 90
                    self.OPLS_Class = 48
                elif Temp_Bond_List == ['C', 'C', 'Cl']:
                    # Chlorobenzene >C-Cl
                    self.OPLS_Type = 205
                    self.OPLS_Class = 48
                elif Temp_Bond_List == ['C', 'C', 'C']:
                    """if not Fullerene:
                        Fullerene = raw_input("Is this a fullerene? True or False: ")
                    if Fullerene:
                        self.OPLS_Type = 907
                        self.OPLS_Class = 48"""
                    # Aromatic Carbon
                    #else:
                    self.OPLS_Type = 90
                    self.OPLS_Class = 48
                elif Temp_Bond_List == ['C', 'C', 'Si']:
                    self.OPLS_Type = 907
                    self.OPLS_Class = 48
                elif Temp_Bond_List == ['C', 'C', 'I']:
                    self.OPLS_Type = 672
                    self.OPLS_Class = 48
                elif Temp_Bond_List == ['C', 'O', 'O']:
                    # Ester -COOR"
                    self.OPLS_Type = 406
                    self.OPLS_Class = 3
                elif Temp_Bond_List == ['C','C','S']:
                    self.OPLS_Type =909
                    self.OPLS_Class = 48
                elif Temp_Bond_List == ['C','H','S']:
                    self.OPLS_Type = 911
                    self.OPLS_Class =48
                elif Temp_Bond_List == ['C','N','O']:
                    self.OPLS_Type = 853
                    self.OPLS_Class = 3
                elif Temp_Bond_List == ['C','C', 'N']:
                    # Quinoline C2
                    self.OPLS_Type = 545
                    self.OPLS_Class = 48
                elif Temp_Bond_List == ['C','C','O']:
                    # Furan C2
                    self.OPLS_Type = 508
                    self.OPLS_Class = 84
                elif Temp_Bond_List == ['C', 'H', 'O']:
                    # Furan C3
                    self.OPLS_Type = 509
                    self.OPLS_Class = 87
        
            elif len(Temp_Bond_List) == 4:
                if Temp_Bond_List == ['C', 'H', 'H', 'I']:
                    self.OPLS_Type = 835
                    self.OPLS_Class = 13
                elif Temp_Bond_List == ['C', 'C', 'H', 'H']:
                    # Alkane C-CH2-C
                    self.OPLS_Type = 81
                    self.OPLS_Class = 13
                elif Temp_Bond_List == ['C', 'H', 'H', 'H']:
                    Acetonitrile = False
                    for bonded_atom in self.Bond_List:
                        if bonded_atom.Element == 'C':
                            if len(bonded_atom.Bond_List) == 2:
                                Acetonitrile = True
                    if Acetonitrile:
                        #Acetonitrile >C-C-N
                        self.OPLS_Type = 696
                        self.OPLS_Class = 13
                    # Alkane C-CH3
                    else:
                        self.OPLS_Type = 80
                        self.OPLS_Class = 13
                elif Temp_Bond_List == ['C', 'C', 'C', 'H']:
                    # Alkane C-CH3
                    self.OPLS_Type = 82
                    self.OPLS_Class = 13
                elif Temp_Bond_List == ['C', 'C', 'H', 'Si']:
                    self.OPLS_Type = 873
                    self.OPLS_Class = 13
                elif Temp_Bond_List == ['Cl', 'Cl', 'Cl', 'H']:
                    # CH-Cl3
                    self.OPLS_Type = 46
                    self.OPLS_Class = 13
                elif Temp_Bond_List == ['H', 'H', 'H', 'O']:
                    Methanol = False
                    for bonded_atom in Atom.Bond_List:
                        if bonded_atom.Element == 'O':
                            for bonded_to_O in bonded_atom.Bond_List:
                                if bonded_to_O.Element == 'H':
                                    Methanol = True
                    if Methanol:
                    #Methanol CH3OH
                        self.OPLS_Type = 99
                        self.OPLS_Class = 13
                    else:
                        # Methyl Ester CH3-O-R
                        self.OPLS_Type = 409
                        self.OPLS_Class = 13
                elif Temp_Bond_List == ['C', 'C', 'C', 'C']:
                    # Alkcane >C<
                    self.OPLS_Type = 84
                    self.OPLS_Class = 13
                elif Temp_Bond_List == ['C', 'C', 'H', 'Si']:
                    # Alkyl silane
                    self.OPLS_Type = 873
                    self.OPLS_Class = 13
                elif Temp_Bond_List == ['H','H','H', 'N']:
                    # 1-methylimidazole
                    self.OPLS_Type = 603
                    self.OPLS_Class = 13
                elif Temp_Bond_List == ['H','H','H', 'H']:
                    # 1-methylimidazole
                    self.OPLS_Type = 83
                    self.OPLS_Class = 13
                elif Temp_Bond_List == ['H','H','H','Si']:
                    self.OPLS_Type = 871
                    self.OPLS_Class = 13
                        
                        

        elif self.Element == "H":
            Bonded_Atom = self.Bond_List[0]
            if Bonded_Atom.OPLS_Type == 90:
                # Aromatic Hydrogen C-H
                print("Then why %d %d" % (self.Atom_ID,Bonded_Atom.OPLS_Class))
                self.OPLS_Type = 91
                self.OPLS_Class = 49
            elif Bonded_Atom.OPLS_Type == 99:
                self.OPLS_Type = 98
                self.OPLS_Class = 46
            elif Bonded_Atom.OPLS_Type == 696:
                self.OPLS_Type = 700
                self.OPLS_Class = 46
            elif Bonded_Atom.OPLS_Type == 96:
                self.OPLS_Type = 97
                self.OPLS_Class = 7
            elif Bonded_Atom.OPLS_Type == 46:
                self.OPLS_Type = 802
                self.OPLS_Class = 46
            elif Bonded_Atom.OPLS_Type == 835:
                self.OPLS_Type = 839
                self.OPLS_Class = 46
            elif Bonded_Atom.OPLS_Class == 13:
                print("The Fuck? %d %d" % (self.Atom_ID,Bonded_Atom.OPLS_Class))
                self.OPLS_Type = 85
                self.OPLS_Class = 46
            elif Bonded_Atom.OPLS_Class == 48:
                print("Okay Cool %d %d" % (self.Atom_ID,Bonded_Atom.OPLS_Class))
                self.OPLS_Type = 91
                self.OPLS_Class = 49
            elif Bonded_Atom.OPLS_Class == 835:
                self.OPLS_Type = 839
                self.OPLS_Class = 46
            else:
                print(Bonded_Atom.OPLS_Type)
                print(Bonded_Atom.OPLS_Class)
                raise Exception("Could not find OPLS Class")
            print(self.OPLS_Type)
            print(self.OPLS_Class)
        elif self.Element == "Cl":
            Bonded_Atom = self.Bond_List[0]
            if Bonded_Atom.OPLS_Class == 48:
                # ChloroBenzene
                self.OPLS_Type = 206
                self.OPLS_Class = 21
            if Bonded_Atom.OPLS_Class == 13:
                # CH-Cl3
                self.OPLS_Type = 47
                self.OPLS_Class = 21

        elif self.Element == "I":
            Bonded_Atom = self.Bond_List[0]
            if Bonded_Atom.OPLS_Class == 48:
                # IodoBenzene
                self.OPLS_Type = 673
                self.OPLS_Class = 66
            if Bonded_Atom.OPLS_Class == 13:
                # Alkyl Iodide
                self.OPLS_Type = 838
                self.OPLS_Class = 66

        elif self.Element == "Si":
            Temp_Bond_List = sorted([ Atomobj.Element for Atomobj in self.Bond_List])
            if Temp_Bond_List == ['C', 'C', 'C', 'C']:
                # Alkyl Silane
                self.OPLS_Type = 866
                self.OPLS_Class = 108

        elif self.Element == "O":
            Temp_Bond_List = sorted([ Atomobj.Element for Atomobj in self.Bond_List])
            if len(Temp_Bond_List) == 2:
                if Temp_Bond_List == ['C', 'C']:
                    # Dialkyl Ether
                    self.OPLS_Type = 408
                    self.OPLS_Class = 20
                if Temp_Bond_List == ['C', 'H']:
                    #Alcohol
                    self.OPLS_Type = 96
                    self.OPLS_Class = 5
            elif len(Temp_Bond_List) ==  1:
                if Temp_Bond_List == ['C']:
                    #   Ester C=O
                    self.OPLS_Type = 407
                    self.OPLS_Class = 4
        elif self.Element == "S":
            Temp_Bond_List = sorted([ Atomobj.Element for Atomobj in self.Bond_List])
            if len(Temp_Bond_List) == 2:
                if Temp_Bond_List == ['C', 'C']:
                    #Sulfide Thiophene
                    self.OPLS_Type = 908
                    self.OPLS_Class = 16
        elif self.Element == "N":
            Temp_Bond_List = sorted([ Atomobj.Element for Atomobj in self.Bond_List])
            if len(Temp_Bond_List) == 3:
                if Temp_Bond_List == ['C','C','C']:
                    # DPP
                    self.OPLS_Type = 847
                    self.OPLS_Class = 107
            if len(Temp_Bond_List) == 1:
                if Temp_Bond_List == ['C']:
                    if len(self.Bond_List[0].Bond_List) == 2:
                        self.OPLS_Type = 694
                        self.OPLS_Class = 18
                    elif len(self.Bond_List) == 3:
                        self.OPLS_Type = 204
                        self.OPLS_Class = 18


# Functions operating on sets of Atom objects


def Find_OPLS_ID(Molecule,Atom, Fullerene):
    # Funcition for assigning OPLS Types and Classes based on connectivity of the atom
    if Atom.Element == "C":
        Temp_Bond_List = sorted([ Molecule.Get_Atom(Atomobj).Element for Atomobj in Atom.Bond_List ])
        if len(Temp_Bond_List) == 2:
            if Temp_Bond_List == ['C','N']:
                #Acetonitrile C->C-N
                Atom.OPLS_Type = 695
                Atom.OPLS_Class = 19
        if len(Temp_Bond_List) == 3:
            if Temp_Bond_List == ['C', 'C', 'H']:
                # Aromatic Carbon
                Atom.OPLS_Type = 90
                Atom.OPLS_Class = 48
            elif Temp_Bond_List == ['C', 'C', 'Cl']:
                # Chlorobenzene >C-Cl
                Atom.OPLS_Type = 205
                Atom.OPLS_Class = 48
            elif Temp_Bond_List == ['C', 'C', 'C']:
                """if not Fullerene:
                    Fullerene = raw_input("Is this a fullerene? True or False: ")
                if Fullerene:
                    Atom.OPLS_Type = 907
                    Atom.OPLS_Class = 48"""
                # Aromatic Carbon
                #else:
                Atom.OPLS_Type = 90
                Atom.OPLS_Class = 48
            elif Temp_Bond_List == ['C', 'C', 'Si']:
                Atom.OPLS_Type = 907
                Atom.OPLS_Class = 48
            elif Temp_Bond_List == ['C', 'C', 'I']:
                Atom.OPLS_Type = 672
                Atom.OPLS_Class = 48
            elif Temp_Bond_List == ['C', 'O', 'O']:
                # Ester -COOR"
                Atom.OPLS_Type = 406
                Atom.OPLS_Class = 3
            elif Temp_Bond_List == ['C','C','S']:
                Atom.OPLS_Type =909
                Atom.OPLS_Class = 48
            elif Temp_Bond_List == ['C','H','S']:
                Atom.OPLS_Type = 911
                Atom.OPLS_Class =48
            elif Temp_Bond_List == ['C','N','O']:
                Atom.OPLS_Type = 853
                Atom.OPLS_Class = 3
            elif Temp_Bond_List == ['C','C', 'N']:
                # Quinoline C2
                Atom.OPLS_Type = 545
                Atom.OPLS_Class = 48
            elif Temp_Bond_List == ['C','C','O']:
                # Furan C2
                Atom.OPLS_Type = 508
                Atom.OPLS_Class = 84
            elif Temp_Bond_List == ['C', 'H', 'O']:
                # Furan C3
                Atom.OPLS_Type = 509
                Atom.OPLS_Class = 87
    
        elif len(Temp_Bond_List) == 4:
            if Temp_Bond_List == ['C', 'H', 'H', 'I']:
                Atom.OPLS_Type = 835
                Atom.OPLS_Class = 13
            elif Temp_Bond_List == ['C', 'C', 'H', 'H']:
                # Alkane C-CH2-C
                Atom.OPLS_Type = 81
                Atom.OPLS_Class = 13
            elif Temp_Bond_List == ['C', 'H', 'H', 'H']:
                Acetonitrile = False
                for bonded_atom in Atom.Bond_List:
                    if Molecule.Get_Atom(bonded_atom).Element == 'C':
                        if len(Molecule.Get_Atom(bonded_atom).Bond_List) == 2:
                            Acetonitrile = True
                if Acetonitrile:
                    #Acetonitrile >C-C-N
                    Atom.OPLS_Type = 696
                    Atom.OPLS_Class = 13
                # Alkane C-CH3
                else:
                    Atom.OPLS_Type = 80
                    Atom.OPLS_Class = 13
            elif Temp_Bond_List == ['C', 'C', 'C', 'H']:
                # Alkane C-CH3
                Atom.OPLS_Type = 82
                Atom.OPLS_Class = 13
            elif Temp_Bond_List == ['C', 'C', 'H', 'Si']:
                Atom.OPLS_Type = 873
                Atom.OPLS_Class = 13
            elif Temp_Bond_List == ['Cl', 'Cl', 'Cl', 'H']:
                # CH-Cl3
                Atom.OPLS_Type = 46
                Atom.OPLS_Class = 13
            elif Temp_Bond_List == ['H', 'H', 'H', 'O']:
                Methanol = False
                for bonded_atom in Atom.Bond_List:
                    if Molecule.Get_Atom(bonded_atom).Element == 'O':
                        for bonded_to_O in Molecule.Get_Atom(bonded_atom).Bond_List:
                            if Molecule.Get_Atom(bonded_to_O).Element == 'H':
                                Methanol = True
                if Methanol:
                #Methanol CH3OH
                    Atom.OPLS_Type = 99
                    Atom.OPLS_Class = 13
                else:
                    # Methyl Ester CH3-O-R
                    Atom.OPLS_Type = 409
                    Atom.OPLS_Class = 13
            elif Temp_Bond_List == ['C', 'C', 'C', 'C']:
                # Alkcane >C<
                Atom.OPLS_Type = 84
                Atom.OPLS_Class = 13
            elif Temp_Bond_List == ['C', 'C', 'H', 'Si']:
                # Alkyl silane
                Atom.OPLS_Type = 873
                Atom.OPLS_Class = 13
            elif Temp_Bond_List == ['H','H','H', 'N']:
                # 1-methylimidazole
                Atom.OPLS_Type = 603
                Atom.OPLS_Class = 13
            elif Temp_Bond_List == ['H','H','H', 'H']:
                # 1-methylimidazole
                Atom.OPLS_Type = 83
                Atom.OPLS_Class = 13
                    
                    

    elif Atom.Element == "H":
        Bonded_Atom = Molecule.Get_Atom(Atom.Bond_List[0])
        if Bonded_Atom.OPLS_Type == 90:
            # Aromatic Hydrogen C-H
            Atom.OPLS_Type = 91
            Atom.OPLS_Class = 49
        elif Bonded_Atom.OPLS_Type == 99:
            Atom.OPLS_Type = 98
            Atom.OPLS_Class = 46
        elif Bonded_Atom.OPLS_Type == 696:
            Atom.OPLS_Type = 700
            Atom.OPLS_Class = 46
        elif Bonded_Atom.OPLS_Type == 96:
            Atom.OPLS_Type = 97
            Atom.OPLS_Class = 7
        elif Bonded_Atom.OPLS_Type == 46:
            Atom.OPLS_Type = 802
            Atom.OPLS_Class = 46
        elif Bonded_Atom.OPLS_Type == 835:
            Atom.OPLS_Type = 839
            Atom.OPLS_Class = 46
        elif Bonded_Atom.OPLS_Class == 13:
            Atom.OPLS_Type = 85
            Atom.OPLS_Class = 46
        elif Bonded_Atom.OPLS_Class == 48:
            Atom.OPLS_Type = 91
            Atom.OPLS_Class = 49
        elif Bonded_Atom.OPLS_Class == 835:
            Atom.OPLS_Type = 839
            Atom.OPLS_Class = 46
        else:
            Atom.OPLS_Type = 85
            Atom.OPLS_Class = 46
    elif Atom.Element == "Cl":
        Bonded_Atom = Molecule.Get_Atom(Atom.Bond_List[0])
        if Bonded_Atom.OPLS_Class == 48:
            # ChloroBenzene
            Atom.OPLS_Type = 206
            Atom.OPLS_Class = 21
        if Bonded_Atom.OPLS_Class == 13:
            # CH-Cl3
            Atom.OPLS_Type = 47
            Atom.OPLS_Class = 21

    elif Atom.Element == "I":
        Bonded_Atom = Molecule.Get_Atom(Atom.Bond_List[0])
        if Bonded_Atom.OPLS_Class == 48:
            # IodoBenzene
            Atom.OPLS_Type = 673
            Atom.OPLS_Class = 66
        if Bonded_Atom.OPLS_Class == 13:
            # Alkyl Iodide
            Atom.OPLS_Type = 838
            Atom.OPLS_Class = 66

    elif Atom.Element == "Si":
        Temp_Bond_List = sorted([ Molecule.Get_Atom(Atomobj).Element for Atomobj in Atom.Bond_List])
        if Temp_Bond_List == ['C', 'C', 'C', 'C']:
            # Alkyl Silane
            Atom.OPLS_Type = 866
            Atom.OPLS_Class = 108

    elif Atom.Element == "O":
        Temp_Bond_List = sorted([ Molecule.Get_Atom(Atomobj).Element for Atomobj in Atom.Bond_List])
        if len(Temp_Bond_List) == 2:
            if Temp_Bond_List == ['C', 'C']:
                # Dialkyl Ether
                Atom.OPLS_Type = 408
                Atom.OPLS_Class = 20
            if Temp_Bond_List == ['C', 'H']:
                #Alcohol
                Atom.OPLS_Type = 96
                Atom.OPLS_Class = 5
        elif len(Temp_Bond_List) ==  1:
            if Temp_Bond_List == ['C']:
                #   Ester C=O
                Atom.OPLS_Type = 407
                Atom.OPLS_Class = 4
    elif Atom.Element == "S":
        Temp_Bond_List = sorted([ Molecule.Get_Atom(Atomobj).Element for Atomobj in Atom.Bond_List])
        if len(Temp_Bond_List) == 2:
            if Temp_Bond_List == ['C', 'C']:
                #Sulfide Thiophene
                Atom.OPLS_Type = 908
                Atom.OPLS_Class = 16
    elif Atom.Element == "N":
        Temp_Bond_List = sorted([ Molecule.Get_Atom(Atomobj).Element for Atomobj in Atom.Bond_List])
        if len(Temp_Bond_List) == 3:
            if Temp_Bond_List == ['C','C','C']:
                # DPP
                Atom.OPLS_Type = 847
                Atom.OPLS_Class = 107
        if len(Temp_Bond_List) == 1:
            if Temp_Bond_List == ['C']:
                if len(Molecule.Get_Atom(Atom.Bond_List[0]).Bond_List) == 2:
                    Atom.OPLS_Type = 694
                    Atom.OPLS_Class = 18
                elif len(Molecule.Get_Atom(Atom.Bond_List[0]).Bond_List) == 3:
                    Atom.OPLS_Type = 204
                    Atom.OPLS_Class = 18


    return Fullerene

def Find_PQEq_Params(self,Mass):
    Element = Reverse_Dict[Mass]
    f = open('pqeq.par.txt','r')
    lines = f.readlines()
    for line in lines:
        if line.strip().split(' ')[0] == Element:
            PQEq_P = float(line.strip().split(' ')[1])
            PQEq_X0 = float(line.strip().split(' ')[2])
            PQEq_J0 = float(line.strip().split(' ')[3])
            PQEq_Z = float(line.strip().split(' ')[4])
            PQEq_Rc = float(line.strip().split(' ')[5])
            PQEq_Rs = float(line.strip().split(' ')[6])
            PQEq_Ks = float(line.strip().split(' ')[7])
    if PQEQ_X0 == null:
        raise ValueError('OPLS Parameter Not Found')
    return [PQEq_X0,PQEq_J0,PQEq_Rc,PQEq_P,PQEq_Z,PQEq_Rs,PQEq_Ks]






