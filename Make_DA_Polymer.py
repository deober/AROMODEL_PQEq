import Molecule
import sys

def main():
    DA_polymer = Molecule.Molecule('/Users/andrewkleinschmidt/DPP-HD.xyz')
    DA_polymer.Set_Up_FF()

if __name__=='__main__': main()
