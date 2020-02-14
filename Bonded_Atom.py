import Atom
import numpy as np
import Ring

class Bonded_Atom(object):

	def __init__(self, Central_Atom, Bonded_Vector, H_Atom_Number,Self_Ring):
		self.Central_Atom = Central_Atom
		self.Bonded_Vector = Bonded_Vector
		self.H_Bond_Vector = self.Bonded_Vector/np.linalg.norm(self.Bonded_Vector)*1.08
		self.H_Atom = Atom.Atom(self.Central_Atom + self.H_Bond_Vector,'H',H_Atom_Number)
		self.Self_Ring = Self_Ring

	def Add_Ring(Bonded_Ring):
		self.Bonded_Ring = Bonded_Ring
		self.Bonded_Ring_ID = Bonded_Ring.Ring_ID
		self.Bonded_Ring_Name = Bonded_Ring.Name

	def Add_Ring_Bonded_Atom():
		for bond_atom in self.Bonded_Ring.Bonded_Atoms:
			if bond_atom.Bonded_Ring == self.Self_Ring:
				self.Interring_Bond_Atom = bond_atom

