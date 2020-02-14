import Analyze_Traj
"""A function that takes in a trajectory, analyzes it, and fits the RDF, bond lengths, angles, dihedrals, and impropers to a tabulated potential function via Iterative Boltzmann Inversion, then fits the tabulated potential to a 9-6 LJ potential"""
def main():
    Script, Traj = sys.argv
    New_Traj = Trajectory(Traj)


if __name__=='__main__': main()
