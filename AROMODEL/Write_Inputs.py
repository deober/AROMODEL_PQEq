#! usr/bin/python

import os
import subprocess
import time

def Write_Orca_SPE(File_Name,Test_Molecule,Method = "RI-MP2",Basis = "cc-pVTZ",Aux_Basis = "cc-pVTZ/C",convergence = "NORMALSCF",nproc = 4,Implicit_Solvent = "CPCM(Chloroform)",Memory = 3000):
	f = open(File_Name,'w')
	if nproc == 1:
		Parallel = ""
	else:
		Parallel = "PAL%d" % nproc
	f.write("! %s %s %s %s %s %s\n%%maxcore %d\n\n*xyz 0 1\n" % (Method,Basis,Aux_Basis,convergence,Parallel,Implicit_Solvent,Memory))
	for atom in Test_Molecule.Atom_List:
		f.write('%s\t%f\t%f\t%f\n' % (atom.Element,atom.Position[0],atom.Position[1],atom.Position[2]))
	f.write("*")
	f.close()

def Write_Orca_Optimize_Geometry(File_Name,Test_Molecule,Method = "RI BP86",Basis = "def2-SVP",Aux_Basis = "def2/J",Dispersion_Correction = "D3BJ",convergence = "TIGHTSCF",nproc = 1,Memory = 3000,H_Only = False,Planarize = False,Planarize_Atoms = [],Bond_Eq = 0.0,Bond_Atoms = []):
	f = open(File_Name,'w')
	if nproc == 1:
		Parallel = ""
	else:
		Parallel = "PAL%d" % nproc
	if H_Only:
		H_Opt = "%geom\noptimizehydrogens true\nend\n\n"
	else:
		H_Opt = ""
	if Planarize:
		Plane = "%geom\nConstraints\n"
		for atoms in Planarize_Atoms:
			Add_Line = "{D %d %d %d %d 0.0 C}\n" % (atoms[0]-1,atoms[1]-1,atoms[2]-1,atoms[3]-1)
			Plane = Plane + Add_Line
		Plane = Plane + "end\nend\n"
	else:
		Plane = ""
	if Bond_Eq == 0.0:
		Bond_Stretch = ""
	else:
		Bond_Stretch = "%%geom Scan\nB %d %d = %.2f,%.2f,12" % (Bond_Atoms[0]-1,Bond_Atoms[1]-1,Bond_Eq-.5,Bond_Eq+.5)
	f.write("! %s %s %s %s %s %s Opt Grid3 FinalGrid5\n%%maxcore %d\n%s%s%s\n*xyz 0 1\n" % (Method,Basis,Aux_Basis,Dispersion_Correction,convergence,Parallel,Memory,H_Opt,Plane,Bond_Stretch))
	for atom in Test_Molecule.Atom_List:
		f.write('%s\t%f\t%f\t%f\n' % (atom.Element,atom.Position[0],atom.Position[1],atom.Position[2]))
	f.write("*")
	f.close()

def Write_Orca_ChelpG(File_Name,Test_Molecule,Method = "RI BP86",Basis = "def2-SVP",Aux_Basis = "def2/J",Dispersion_Correction = "D3BJ",convergence = "TIGHTSCF",nproc = 4,Memory = 3000,H_Only = False,Planarize = False,Planarize_Atoms = []):
	f = open(File_Name,'w')
	if nproc == 1:
		Parallel = ""
	else:
		Parallel = "PAL%d" % nproc
	f.write("! RKS B3LYP 6-31+G** NormalSCF NOSOSCF CHELPG\n%%maxcore %d\n\n*xyz 0 1\n" % (Memory))
	for atom in Test_Molecule.Atom_List:
		f.write('%s\t%f\t%f\t%f\n' % (atom.Element,atom.Position[0],atom.Position[1],atom.Position[2]))
	f.write("*")
	f.close()

def Write_QChem_SPE(File_Name,Test_Molecule,Exchange_Method = "HF",Correlation_Method = "pRIMP2",Basis = "cc-pvtz",Method = "rimp2",Aux_Basis = "rimp2-cc-pvtz",Memory = 12000,convergence = 6,Implicit_Solvent_Method = "PCM",Implicit_Solvent_Dielectric = 4.9):
	f = open(File_Name,'w')
	f.write("$molecule\n\t0 1\n")
	for atom in Test_Molecule.Atom_List:
		f.write('%s\t%f\t%f\t%f\n' % (atom.Element,atom.Position[0],atom.Position[1],atom.Position[2]))
	f.write("$end\n\n$rem\n\tJOBTYPE\t\t\tSP\n\tEXCHANGE\t%s\n\tCORRELATION%s\t\n\tBASIS\t\t\t%s\n\tMETHOD\t\t\t%s\n\tAUX_BASIS\t%s\nSOLVENT_METHOD\t%s\n\tPURECART\t11111\n\tSYMMETRY\tfalse\n\tMEM_TOTAL\t%d\n\tSCF_CONVERGENCE = %d\n\tTHRESH=%d\nNBO\t1\n$end" % (Exchange_Method,Correlation_Method,Basis,Method,Aux_Basis,Implicit_Solvent_Method,Memory,convergence,convergence + 4))
	if Implicit_Solvent != "":
		f.write("\n\n$solvent\n\tDielectric\t\t\t%.1f\n$end" % Implicit_Solvent_Dielectric)
	f.close()

def Write_QChem_Optimize_Geometry(File_Name,Test_Molecule,Basis = "def2-SVP",Method = "BP86",Approximations = "\n\tRI-J\t\t\tTRUE",Aux_Basis = "def2-SVP-J",Memory = 12000,convergence = 6,Implicit_Solvent_Method = "PCM",Implicit_Solvent_Dielectric = 4.9):
	f = open(File_Name,'w')
	f.write("$molecule\n\t0 1\n")
	for atom in Test_Molecule.Atom_List:
		f.write('%s\t%f\t%f\t%f\n' % (atom.Element,atom.Position[0],atom.Position[1],atom.Position[2]))
	f.write("$end\n\n$rem\n\tJOBTYPE\t\t\tOPT\n\tBASIS\t\t\t%s\n\tMETHOD\t\t\t%s\n\tDFT_D\t\t\tD3_BJ\n\tAUX_BASIS\t%s\n\tPURECART\t11111\n\tSYMMETRY\tfalse\n\tMEM_TOTAL\t%d\n\tSCF_CONVERGENCE = %d\n\tTHRESH=%d\n$end" % (Basis,Method,Aux_Basis,Memory,convergence,convergence + 4))
	if Implicit_Solvent != "":
		f.write("\n\n$solvent\n\tDielectric\t\t\t%.1f\n$end" % Implicit_Solvent_Dielectric)
	f.close()


def Write_NWChem_SPE(File_Name,Test_Molecule,Name,Method = "rimp2",charge = 0,Basis = "cc-pVTZ",Aux_Basis = "cc-pVTZ",Implicit_Solvent_Dielectric = 4.9,ChelpG = False,DFT_Method = 'xc becke88 perdew86\n vdw 3'):
	f = open(File_Name,'w')
	f.write("start %s\ncharge %d\ngeometry\n" % (Name,charge))
	for atom in Test_Molecule.Atom_List:
		f.write(' %s\t%f\t%f\t%f\n' % (atom.Element,atom.Position[0],atom.Position[1],atom.Position[2]))
	f.write("end\nbasis\n * library %s\nend\n" % (Basis))
	if Method == "rimp2":
		f.write("basis \"ri-mp2 basis\"\n * library %s\nend\n" % (Aux_Basis))
	if Method == "dft":
		f.write("dft\n %s\nend" % DFT_Method)
	if Implicit_Solvent_Dielectric != 0.0:
		f.write("cosmo\n dielec %d\nend\n" % (Implicit_Solvent_Dielectric))
	f.write("task %s energy" % Method)
	if ChelpG:
		f.write("\ntask esp")
	f.close()
	return

def Write_NWChem_Optimize_Geometry(File_Name,Test_Molecule,Name,Method = "dft",charge = 0,Basis = "def2-SVP",Aux_Basis = "cc-pVTZ",Implicit_Solvent_Dielectric = 0.0):
	return

def Write_LAMMPS_Minimize(File_Name,Data_File,Name,Coul_Cutoff):
	f = open(File_Name,'w')
	f.write("units\t\t\treal\natom_style\tfull\npair_style\tlj/cut/coul/cut 10.0 %.1f\nbond_style\tharmonic\nangle_style\tharmonic\ndihedral_style\topls\nspecial_bonds\tlj/coul 0 0 0.5\nimproper_style\tcvff\npair_modify mix geometric\nread_data %s\n\ntimestep 1.0\nminimize 1.0e-8 1.0e-8 1000 100000\nwrite_data %s_Optimized_Geometry.data" % (Coul_Cutoff,Data_File,Name))
	f.close()

def Write_LAMMPS_Nonbonded(File_Name,Data_File,Name,Coul_Cutoff):
	f = open(File_Name,'w')
	f.write("units\t\t\treal\natom_style\tfull\npair_style\tlj/cut/coul/cut 10.0 %.1f\nbond_style\tharmonic\nangle_style\tharmonic\ndihedral_style\topls\nspecial_bonds\tlj/coul 0 0 0.5\nimproper_style\tcvff\npair_modify mix geometric\nread_data %s\n\ntimestep 1.0\ngroup frozen molecule == 1\ngroup solvent molecule != 1\ncompute nonbond frozen pe/atom pair\ncompute totalnonbond frozen reduce sum c_nonbond\nthermo 10\nthermo_style custom step temp density press etotal pe ke c_nonbond[*] c_totalnonbond\nneighbor 2.0 bin\nrun_style verlet\nfix 1 solvent npt temp 300.0 300.0 100.0 iso 1.0 1.0 1000.0\ndump 1 all custom 1000 %s.lammpstrj id type mol x y z ix iy iz\nrun 2000000\nwrite_data Final_%s.data" % (Coul_Cutoff,Data_File,Name,Name))
	f.close()

def Write_LAMMPS_Single_Point_Energy_No_Nonbonded(File_Name,Data_File,Name,Coul_Cutoff):
	f = open(File_Name,'w')
	f.write("units\t\t\treal\natom_style\tfull\npair_style\tlj/cut/coul/cut 10.0 %.1f\nbond_style\tharmonic\nangle_style\tharmonic\ndihedral_style\topls\nspecial_bonds\tlj/coul 0 0 0.5\nimproper_style\tcvff\npair_modify mix geometric\nread_data %s.data\n\ntimestep 1.0\ncompute bonded all pe bond angle dihedral improper\nthermo 10\nthermo_style custom step temp density press etotal pe ke c_bonded\nneighbor 2.0 bin\nrun_style verlet\nfix 1 solvent npt temp 300.0 300.0 100.0 iso 1.0 1.0 1000.0\ndump 1 all custom 1000 %s.lammpstrj id type mol x y z ix iy iz\nrun 2000000\nwrite_data Final_%s.data" % (Coul_Cutoff,Name,Name,Name))
	f.close()

def Write_LAMMPS_PQEq(File_Name,Data_File,Name,Final_Temperature,Full_System,Coul_Cutoff=10,Implicit_Solvent=False):
	'''
	Function to write LAMMPS input file with polarizable charge equilibration (PQEq) charge definitions. PQEq pair coefficients HAVE 
	to appear in the LAMMPS control file (in.<systemName>) and NOT in the data file (<systemName>.data)
	
	File_Name = input file name 
	Data_File = .data file name generated by "System"
	Name = suffix name to append to trajectory names, etc.
	Full_System = Object from System class in System.py: need this to collect PQEq parameters and write pair coefficients to Lammps control file
	'''

	f = open(File_Name,'w')
	f.write('#LAMMPS input file')
	f.write('units\t\t\treal\natom_style\t\t\tpqeq\nboundary\t\t\tp p p\nbond_style\t\t\tharmonic\nangle_style\t\t\tharmonic\ndihedral_style\t\t\thybrid opls multi/harmonic8\n')
	f.write('improper_style\t\t\tcvff\nkspace_style\t\t\tnone')

	#Add option for implicit solvent
	f.write('dielectric 1\n')

	#data source
	f.write('read_data %s\n' % Data_File)

	f.write('special_bonds\t\t\tlj 1e-08 1e-08 0.5 coul 1e-08 1e-08 0.5\npair_style\t\t\thybrid/overlay coul/pqeqgauss 0.00 12.50 lj/cut\t\t\t%f\n'% (Coul_Cutoff))
	f.write('\nread_data\t\t\t%s\nthermo_style\t\t\tmulti\n' % (Data_File))

	#coul/pqeqgauss pair coefficient informaiton
	f.write('#coul/pqeqgauss\n')
	f.write('pair_coeff           * * coul/pqeqgauss 0.0000 0.0000 0.000 0 0.00 0.000 0.0 0 #dummy\n')
	for AtomParam, PQEq_Param in Full_System.Atom_Params, Full_System.PQEq_Params:
		f.write('pair_coeff\t\t\t%d lj/cut/coul/long %.3f %.3f\n' % (i,Atom_Param[2],Atom_Param[3]))
		f.write('pair_coeff\t\t\t%d %d coul/pqeqgauss %.4f %.4f %.3f %d %.2f %.3f %.1f %d\n' % tuple([i]+[i]+PQEq_Param) )

	#pair mixing rules, modifying output format, etc
	f.write('pair_modify\t\t\tmix geometric\nneighbor\t\t\t2.0 bin\nneigh_modify\t\t\tevery 2 delay 4 check yes\nthermo_style\t\t\tmulti\nthermo_modify\t\t\tline multi format float %14.6f flush yes\n')
	f.write('fix\t\t\tpqeq all pqeq method 2 nevery 1 charge 0.0 tolerance 1.0e-6 damp 1.0\n')
	
	#minimization steps
	f.write('\n\n\nprint\t\t\t.\nprint\t\t\t____________________________\nprint\t\t\t"500 steps Conjugate Gradient minimization"\nprint\t\t\t____________________________\nprint\t\t\t.')
	f.write('dump\t\t\t1 all custom 25 ${sname}.min.lammps id type xu yu zu vx vy vz q\nthermo\t\t\t10\nmin_style\t\t\tsd\nminimize\t\t\t1.0e-4 1.0e-4 500 5000\nmin_style\t\t\tcg\nminimize\t\t\t1.0e-4 1.0e-4 500 5000\nminimize\t\t\t1.0e-4 1.0e-4 500 5000\nundump\t\t\t1')
	f.write('variable\t\t\tinput index %s\nvariable\t\t\tsname index %s\nvariable rtemp equal %f\n' % File_Name,Name,Final_Temperature)
	restOfTemplate = '''print                .
	print                =====================================
	print                "NVT dynamics to heat system"
	print                =====================================
	print                .
						
	reset_timestep       0
	timestep             1.0
	fix                  shakeH all shake 0.0001 20 500 m 1.0079
	velocity             all create 0.0 12345678 dist uniform
	thermo               100
	thermo_style         multi
	dump                 1 all custom 1000 ${sname}.heat.lammpstrj id type xu yu zu vx vy vz q
	fix                  4 all nvt temp 1.0 ${rtemp} 100.0
	run                  10000
	unfix                4
	undump               1
						
	print                .
	print                ================================================
	print                "NPT dynamics with an isotropic pressure of 1atm."
	print                ================================================
	print                .
						
	timestep             2.0
	fix                  2 all npt temp ${rtemp} ${rtemp} 100.0 iso 1.0 1.0 2000.0
	thermo               100
	thermo_style         multi
	dump                 1 all custom 5000 ${sname}.${rtemp}K.npt.lammps id type xu yu zu vx vy vz q
	run                  500000 # run for 15 ns
						
	variable             latx equal lx
	variable             laty equal ly
	variable             latz equal lz
	fix                  lxavg all ave/time 1 250000 250000 v_latx
	fix                  lyavg all ave/time 1 250000 250000 v_laty
	fix                  lzavg all ave/time 1 250000 250000 v_latz
	run                  510000 # run for 15 ns
	variable             xavg equal f_lxavg
	variable             yavg equal f_lyavg
	variable             zavg equal f_lzavg
	undump               1
	unfix                2
	print                "current cell: ${latx} ${laty} ${latz} cell avgs: ${xavg} ${yavg} ${zavg}"
	print                "deforming cell"
	fix                  2 all nvt temp ${rtemp} ${rtemp} 100.0 nreset 10000
	dump                 1 all custom 5000 ${sname}.${rtemp}K.deform.lammps id type xu yu zu vx vy vz q
	fix                  1 all deform 100 x final 0 ${xavg} y final 0 ${yavg} z final 0 ${zavg} units box
	undump               1
	unfix                lxavg
	unfix                lyavg
	unfix                lzavg
	run                  100000
	unfix                1
	unfix                2
	fix                  2 all nvt temp ${rtemp} ${rtemp} 100.0 nreset 10000
	dump                 1 all custom 5000 ${sname}.${rtemp}K.nvt.lammps id type xu yu zu vx vy vz q
	run                  500000
	undump               1
	unfix                2
	reset_timestep       1000000
						
	print                .
	print                ================================================
	print                "NVT production dynamics "
	print                ================================================
	print                .
						
	fix                  2 all nvt temp ${rtemp} ${rtemp} 100.0 tloop 10 ploop 10
	thermo               100
	restart              50000 ${sname}.${rtemp}K.*.restart
	dump                 1 all custom 1000 ${sname}.${rtemp}K.prod.lammps id type xu yu zu vx vy vz q
	run                  5000000 # run for 15 ns
	unfix                2
	undump               1
	'''


	


	f.close()



	


def Write_LAMMPS_Test_Condensed():
	return
