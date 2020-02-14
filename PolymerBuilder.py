import sys
import os
import random
import Molecule

#This script can build a linear polymer, either a random copolymer or a block copolymer. Polymers will be regioregular and run along the x-axis.

class Atom(object):
    #This class defines an atom. Atom_Num signifies the atom's identifying number (as used in LAMMPS), Atom_Type is the element, x, y, and z are the respective coordinates.
    def __init__(self,Atom_Num,Atom_Type,x,y,z,Start,Link):
        Mass_List = {'H':1.01,'C':12.011,'N':14.007, 'O':15.999,'F':18.998, 'S':32.06, 'Se':78.96}
        self.Atom_Num = Atom_Num
        self.Atom_Type = Atom_Type
        self.Mass = Mass_List[self.Atom_Type]
        self.Coord = [x,y,z]
        self.Start = Start
        self.Link = Link

    def print_atom(self):
        print 'A %s atom, mass %f, at position %f, %f, %f' % (self.Atom_Type,self.Mass,self.Coord[0],self.Coord[1],self.Coord[2])

class Monomer(object):
    #Describes a full, unique monomer. Atom_List contains all the atoms present in the monomer, Bond_List contains the Atom_Num of bonded atoms, ID is the number used to represent the monomer in Mono_Order (begins at 0), Start_ and End_Atom give the coordinates of the connective atoms to adjacent monomers. End_Atom will not actually be present in the Monomer (the connecting atom of an adjacent monomer).
    
    def __init__(self,Name,Atom_List,Bond_List,ID,Start_Atom,End_Atom,Start_Group):
        self.Name = Name
        self.Atom_List = Atom_List
        self.Bond_List = Bond_List
        self.ID = ID
        self.Start_Atom = Start_Atom
        self.End_Atom = End_Atom
        self.Start_Group = Start_Group

Script, Deg_Poly, Block, Mono, Ratio_dum, Name = sys.argv #Deg_Poly is the total number of monomers in the polymer, Block is whether or not the polymer is a block copolymer (0 false, 1 true), Mono is the xyz files of the different monomers being used (comma separated, no spaces or brackets), Ratio_dum is the ratio of each monomer in the copolymer (comma separated, no spaces or brackets, decimal form, will get transformed into Ratio), Name is the name of the folder/file that will be created to write data into.

Deg_Poly = int(Deg_Poly)
Block = int(Block)
Mono = Mono.split(',')
Ratio_dum = Ratio_dum.split(',')
Ratio = []

for r in Ratio_dum:
    Ratio.append(float(r))

Monomer_List = [] #List of Monomer objects, extracted from the input files

for i,m in enumerate(Mono):
    Atom_List = []
    Bond_List = []
    Num_Bonds = 0
    j = 1
    f = open(m,'r')
    lines = f.readlines()
    #f2 = open('%s.txt' %m.split('.')[0], 'r')
    #bond_lines = f2.readlines()
    for l in lines:
        l = l.strip().split()
        if len(l) ==4:
            Atom_List.append(Atom(j,l[0].strip(),float(l[1].strip()),float(l[2].strip()),float(l[3].strip()),0,0))
            j += 1
        if len(l) ==5:
            if l[4] == '*':
                Start_Atom = [float(l[1].strip()),float(l[2].strip()),float(l[3].strip())]
                Atom_List.append(Atom(j,l[0].strip(),float(l[1].strip()),float(l[2].strip()),float(l[3].strip()),1,0))
                j += 1
            if l[4] == '**':
                End_Atom = [l[0].strip(),float(l[1].strip()),float(l[2].strip()),float(l[3].strip())]

            if l[4] == '***':
                Atom_List.append(Atom(j,l[0].strip(),float(l[1].strip()),float(l[2].strip()),float(l[3].strip()),0,1))
            
            if l[4] == '****':
                Start_Group = Atom(j,l[0].strip(),float(l[1].strip()),float(l[2].strip()),float(l[3].strip()),0,0)
    """for l in bond_lines:
        l = l.strip().split()
        Bond_List.append([int(l[0].strip()),int(l[1].strip()),int(l[2].strip())])
        if int(l[0].strip()) > Num_Bonds:
            Num_Bonds = int(l[0].strip())"""
    f.close()
    #f2.close()
    Monomer_List.append(Monomer(m.split('.')[0],Atom_List,Bond_List,i,Start_Atom,End_Atom,Start_Group))

Mono_Order = [] #Order of monomers, randomly generated, with numbers corresponding to the order put in Mono (starting with 0)


if Block:
    for i in range(len(Monomer_List)):
        for j in range(int(Ratio[i]*Deg_Poly)):
            Mono_Order.append(i)

else:
    for i in range(Deg_Poly):
        random.seed()
        r = random.random()
        p = 0
        for j,q in enumerate(Ratio):
            if r < q+p:
                Mono_Order.append(j)
                break
            p += q

Start_Position = [0.0,0.0,0.0]
Mol_Atom_List = [] #List of atoms for the final polymer, will be used for final file
Mol_Bond_List = [] #List of bonds for the final polymer, will be used for final file
k = 1 # counter for atoms
l = 1 # counter for bonds

Last_Atom = -1
Last_Atom_dum = -1
Start_Atom = -1


for q in Monomer_List:
    if Mono_Order[0] == q.ID:
        Mol_Atom_List.append(Atom(k,q.Start_Group.Atom_Type,q.Start_Group.Coord[0],q.Start_Group.Coord[1],q.Start_Group.Coord[2],0,0))
        k+=1
        break
for z,m in enumerate(Mono_Order):
    for q in Monomer_List:
        if m == q.ID:
            for b in q.Bond_List:
                Mol_Bond_List.append([l,b[1]+k-1,b[2]+k-1])
                l+=1
            for a in q.Atom_List:
                Mol_Atom_List.append(Atom(k,a.Atom_Type,a.Coord[0]+Start_Position[0],a.Coord[1]+Start_Position[1],a.Coord[2]+Start_Position[2],0,0))
                if a.Start:
                    Start_Atom = k
                if a.Link:
                    Last_Atom_dum = k
                k += 1
            if Last_Atom != -1:
                Mol_Bond_List.append([l,Last_Atom,Start_Atom])
                l+=1
            if z == len(Mono_Order) - 1:
                Mol_Atom_List.append(Atom(k,q.End_Atom[0],q.End_Atom[1] + Start_Position[0],q.End_Atom[2] + Start_Position[1],q.End_Atom[3] + Start_Position[2],0,0))
            Last_Atom = Last_Atom_dum
            Start_Position[0] += q.End_Atom[1]
            Start_Position[1] += q.End_Atom[2]
            Start_Position[2] += q.End_Atom[3]
            break

#os.system('mkdir ~/%s' % Name)

f = open('/Users/andrewkleinschmidt/AROMODEL/AROMODEL/%s.xyz' % Name,'w')
f.write('%d\n\n' % Mol_Atom_List[-1].Atom_Num)

for m in Mol_Atom_List:
    f.write('%s\t%f\t%f\t%f\n' % (m.Atom_Type,m.Coord[0],m.Coord[1],m.Coord[2]))
                          
f.write('\n')

f.close()

"""Polymer = Molecule.Molecule('%s.xyz' % Name)
Polymer.Set_Up_FF(local = False)
for m in Mol_Bond_List:
    f.write('%d\t%d\t%d\n' % (m[0],m[1],m[2]))
"""
