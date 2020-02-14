import Atom
import numpy as np
import Ring
import math

class Bonded_Atom(object):

	def __init__(self, Central_Atom, Bonded_Vector, Bonded_Atom_ID,H_Atom_Number,Self_Ring,Same_Ring_Bonded_Atom_List):
		self.Central_Atom = Central_Atom
		self.Bonded_Vector = Bonded_Vector
		self.Bonded_Atom_ID = Bonded_Atom_ID
		self.H_Bond_Vector = self.Bonded_Vector/np.linalg.norm(self.Bonded_Vector)*1.08
		self.H_Atom = Atom.Atom(self.Central_Atom.Position + self.H_Bond_Vector,'H',H_Atom_Number)
		self.Self_Ring = Self_Ring
		self.Same_Ring_Angle_List = []
		self.Is_Linked = False
		for b_atom_id in Same_Ring_Bonded_Atom_List:
			temp_atom = self.Self_Ring.Get_Atom(b_atom_id)
			temp_angle = math.acos(np.dot((self.Central_Atom.Position - temp_atom.Position)/np.linalg.norm(self.Central_Atom.Position - temp_atom.Position),self.Bonded_Vector/np.linalg.norm(self.Bonded_Vector)))
			self.Same_Ring_Angle_List.append((temp_atom,temp_angle))

	def Add_Ring(self,Bonded_Ring):
		self.Bonded_Ring = Bonded_Ring
		self.Bonded_Ring_ID = Bonded_Ring.Ring_ID
		self.Bonded_Ring_Name = Bonded_Ring.Name
		self.Is_Linked = True

	def Add_Bonded_Atom(self,Bonded_Atom):
		for bond_atom in self.Bonded_Ring.Bonded_Atoms:
			if bond_atom.Bonded_Ring == self.Self_Ring:
				self.Interring_Bond_Atom = bond_atom

	def Check_Alignment(self):
		average_angle = 0
		Check_Average_Angle_List = []
		for atom_and_angle in self.Same_Ring_Angle_List:
			bonded_vector = self.Bonded_Vector
			bonded_vector = bonded_vector/np.linalg.norm(bonded_vector)
			self_vector = self.Central_Atom.Position - atom_and_angle[0].Position
			self_vector = self_vector/np.linalg.norm(self_vector)
			Check_Average_Angle_List.append(atom_and_angle[1] - math.acos(np.dot(self_vector,bonded_vector)))
			average_angle += atom_and_angle[1] - math.acos(np.dot(self_vector,bonded_vector))
		average_angle = average_angle / len(atom_and_angle)
		for angle1 in Check_Average_Angle_List:
			for angle2 in Check_Average_Angle_List:
				if abs(angle1 - angle2) > 5:
					print(angle1)
					print(angle2)
					raise Exception("Rotation Angles too divergent")

		return average_angle

