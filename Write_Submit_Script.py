#! usr/bin/python
import math

def Write_Qchem(f,nproc,In_File,Name):
	f.write('export QCMACHINEFILE=`generate_pbs_nodefile`\ncat $QCMACHINEFILE|uniq >.tmp\nmv .tmp $QCMACHINEFILE\nmodule load qchem\nexport QCSCRATCH=/scratch/$USER/$SLURM_JOBID\nexport QCMACHINEFILE=`generate_pbs_nodefile`\ncat $QCMACHINEFILE|uniq >.tmp\nmv .tmp $QCMACHINEFILE\nmodule load qchem\nexport QCSCRATCH=/scratch/$USER/$SLURM_JOBID\nqchem -nt %d %s %s.out' % (nproc,In_File,Name))

def Write_Qchem_Run_Only(f,nproc,In_File,Name):
	f.write('&\nqchem -nt %d %s %s.out' % (nproc,In_File,Name))

def Write_Orca(f,nproc,In_File,Name,Executable_Path,OMP_Path):
	f.write('export PATH=%s:%s/bin:$PATH\nexport RSH_COMMAND="/usr/bin/ssh -x"\nexport LD_LIBRARY_PATH=%s/lib\n%s/orca %s >> %s.out\n' % (OMP_Path,Executable_Path,OMP_Path,Executable_Path,In_File,Name))

def Write_Orca_Run_Only(f,nproc,In_File,Name,Executable_Path,OMP_Path):
	f.write('orca %s >> %s.out\n' % (In_File,Name))

def Write_LAMMPS(f,nproc,In_File,Name):
	f.write('module load lammps\n\nsrun -np %d lammps -in %s -log log.%d' % (nproc,In_File,Name))

def Write_LAMMPS_Run_Only(f,nproc,In_File,Name):
	f.write('srun -np %d lammps -in %s -log log.%d' % (nproc,In_File,Name))

def Write_NWChem():
	return

def Write_SLURM(File_Name,In_File,Name,nproc,Cluster_Location,Job_Type,queue="shared",proc_per_node=32,account = "m3047",walltime = 2,Executable_Path = "",OMP_Path = "",constraint='cori'):
	f = open(File_Name,'w')
	"""nodes = math.floor(nproc/proc_per_node)
	if nodes == 0:
		nodes = 1
		queue = "shared"
		ppn = nproc
	else:
		queue = "compute"
		ppn = proc_per_node"""
	f.write('#!/bin/bash\n#SBATCH --job-name="%s"\n#SBATCH --output=%s\n#SBATCH --partition=%s\n#SBATCH --nodes=%d\n#SBATCH --ntasks-per-node=%d\n#SBATCH -A %s\n#SBATCH --export=ALL\n#SBATCH -t %d:00:00\n\ncd %s\n' % (Name,Name,queue,nodes,ppn,account,walltime,Cluster_Location))
	if Job_Type == "QChem":	
		if Executable_Path == "":
			Write_Qchem(f,nproc,In_File,Name)

	if Job_Type == "Orca":
		if Executable_Path == "":
			raise Exception("Executable_Path required for Orca")
		if OMP_Path == "":	
			raise Exception("OMP_Path required for Orca")
		Write_Orca(f,nproc,In_File,Name,Executable_Path,OMP_Path)


	if Job_Type == "LAMMPS":
		if Executable_Path == "":
			Write_LAMMPS(f,nproc,In_File,Name)
			
	f.close()

def Write_SLURM_Batch(File_Name,In_File_List,Name,nproc,Cluster_Location,Job_Type,queue="",proc_per_node=32,account = "m3047",walltime = 2,Executable_Path = "",OMP_Path = "",constraint = 'knl'):
	f = open(File_Name,'w')
	"""nodes = math.floor(nproc/proc_per_node)
	if nodes == 0:
		nodes = 1
		queue = "shared"
		ppn = nproc
	else:
		queue = "compute"
		ppn = proc_per_node"""
	f.write('#!/bin/bash\n#SBATCH --job-name="%s"\n#SBATCH --output=%s\n#SBATCH --partition=%s\n#SBATCH --nodes=%d\n#SBATCH --ntasks-per-node=%d\n#SBATCH -A %s\n#SBATCH --export=ALL\n#SBATCH -t %d:00:00\n\ncd %s\n' % (Name,Name,queue,nodes,ppn,account,walltime,Cluster_Location))
	if Job_Type == "QChem":	
		if Executable_Path == "":
			Write_Qchem(f,nproc,In_File_List[0],Name)
		for in_file in In_File_List[1:]:
			Write_Qchem_Run_Only(f,nproc,in_file,Name)
		f.write("wait")

	if Job_Type == "Orca":
		if Executable_Path == "":
			raise Exception("Executable_Path required for Orca")
		if OMP_Path == "":	
			raise Exception("OMP_Path required for Orca")
		Write_Orca(f,nproc,In_File_List[0],Name,Executable_Path,OMP_Path)
		for in_file in In_File_List[1:]:
			Write_Orca_Run_Only(f,nproc,in_file,Name,Executable_Path,OMP_Path)
		f.write("wait")


	if Job_Type == "LAMMPS":
		if Executable_Path == "":
			Write_LAMMPS(f,nproc,In_File_List[0],Name)
		for in_file in In_File_List[1:]:
			Write_LAMMPS_Run_Only(f,nproc,in_file,Name)
		f.write("wait")
			
	f.close()

def Write_TORQUE(File_Name,In_File,Name,nproc,Cluster_Location,Job_Type,proc_per_node=28,walltime = 2,Executable_Path = "",OMP_Path = "",queue = "condo"):
	f = open(File_Name,'w')
	nodes = math.floor(nproc/proc_per_node)
	if nodes == 0:
		nodes = 1
		ppn = nproc
	else:
		ppn = proc_per_node
	f.write('#!/bin/bash\n#PBS -N %s\n#PBS -o %s\n#PBS -q %s\n#PBS -l nodes=%d:ppn=%d\n#PBS -l walltime=%d:00:00\n\ncd %s\n\n' % (Name,Name,queue,nodes,ppn,walltime,Cluster_Location))
	if Job_Type == "QChem":	
		if Executable_Path == "":
			Write_Qchem(f,nproc,In_File,Name)


	if Job_Type == "Orca":
		if Executable_Path == "":
			raise Exception("Executable_Path required for Orca")
		if OMP_Path == "":	
			raise Exception("OMP_Path required for Orca")
		Write_Orca(f,nproc,In_File,Name,Executable_Path,OMP_Path)


	if Job_Type == "LAMMPS":
		if Executable_Path == "":
			Write_LAMMPS(f,nproc,In_File,Name)
			
	f.close()