#! usr/bin/python
import os
import subprocess
import time

def Quantum_Orca(file):
#Returns True if an Orca job has finished running, and False if not
	Finished = False
	Energy = 0.0
	Ended_Flag = False
	f = open("./%s" % (file),'r')
	lines = f.readlines()
	for line in lines:
		if len(line.strip().split()) > 2:
			if line.strip().split()[0] == "FINAL" and line.strip().split()[1] == "SINGLE":
				Energy = float(line.strip().split()[4]) * 627.5
			if line.strip().split()[0] == "TOTAL" and line.strip().split()[1] == "RUN" and line.strip().split()[2] == "TIME:":
				Ended_Flag = True
	f.close()
	if Energy != 0.0 and Ended_Flag:
		Finished = True

	return Finished

def Opt_Orca(file):
#Returns True if an Orca job has finished running, and False if not
	Finished = False
	Ended_Flag = False
	f = open("./%s" % (file),'r')
	lines = f.readlines()
	for line in lines:
		if len(line.strip().split()) > 2:
			if line.strip().split()[0] == "TOTAL" and line.strip().split()[1] == "RUN" and line.strip().split()[2] == "TIME:":
				Ended_Flag = True
	f.close()
	if Ended_Flag:
		Finished = True

	return Finished

def Quantum_QChem(file):
#Returns True if a Qchem job has finished running, and False if not
	Finished = False
	Energy = 0.0
	Ended_Flag = False
	f = open("./%s" % (file),'r')
	lines = f.readlines()
	for line in lines:
		if len(line.strip().split()) > 1:
			if line.strip().split()[0] == "RIMP2" and line.strip().split()[1] == "total":
				Energy = float(line.strip().split()[4]) * 627.5
			if line.strip().split()[0] == "Total" and line.strip().split()[1] == "job" and line.strip().split()[2] == "time":
				Ended_Flag = True
	f.close()
	if Energy != 0.0 and Ended_Flag:
		Finished = True

	return Finished

def Quantum_NWChem():
#Returns True if an NWChem job has finished running, and False if not
	return False

def Energy_Orca(Analyze_File,Folder_Name):
#Returns the single point energy in kcal/mol for an Orca calculation
	Energy = 0.0
	f = open("./%s/%s" % (Folder_Name,Analyze_File),'r')
	lines = f.readlines()
	for line in lines:
		if len(line.strip().split()) > 2:
			if line.strip().split()[0] == "FINAL" and line.strip().split()[1] == "SINGLE":
				Energy = float(line.strip().split()[4]) * 627.5
	f.close()

	return Energy

def Loewdin_Charges_Orca(Analyze_File,Folder_Name):
#Returns the Loewdin point charges for each atom for an Orca calculation
	charges = False
	Charge_List = []
	f = open("./%s/%s" % (Folder_Name,Analyze_File),'r')
	lines = f.readlines()
	count = 0
	for line in lines:
		if len(line.strip().split()) > 1:
			if line.strip().split()[0] == "LOEWDIN" and line.strip().split()[1] == "ATOMIC":
				charges = True
				count = 0
			if count > 1:
				charges = False
			if charges and line.strip().split()[0].isdigit():
				Charge_List.append(float(line.strip().split()[-1]))
		else:
			count += 1
	f.close()

	return Charge_List

def Energy_Qchem(Analyze_File,Folder_Name):
#Returns the single point energy in kcal/mol for a QChem calculation 
	Energy = 0.0
	f = open("./%s/%s" % (Folder_Name,Analyze_File),'r')
	lines = f.readlines()
	for line in lines:
		if len(line.strip().split()) > 1:
			if line.strip().split()[0] == "RIMP2" and line.strip().split()[1] == "total":
				Energy = float(line.strip().split()[4]) * 627.5
	f.close()

	return Energy

def NBO_Charges_QChem(Analyze_File,Folder_Name):
#Returns the NBO point charges for each atom for a QChem calculation
	charges = False
	Charge_List = []
	f = open("./%s/%s" % (Folder_Name,Analyze_File),'r')
	lines = f.readlines()
	for line in lines:
		if len(line.strip().split()) > 1:
			if line.strip().split()[0] == "Ground-State" and line.strip().split()[0] == "Mulliken":
				charges = True
			if line.strip().split()[0] == "Sum":
				charges = False
			if charges and line.strip().split()[0].isdigit():
				Charge_List.append(float(line.strip().split()[2]))
	f.close()

	return Charge_List

def Energy_NWChem(Analyze_File,Folder_Name):
#Returns the single point energy in kcal/mol for an NWChem calculation 
	return

def Energy_Averaged_LAMMPS(Analyze_File,Folder_Name,Start_Point):
#Returns the average of the energy (in kcal/mol) in the last column of the thermo output in a LAMMPS log file from the timestep Start_Point to the end
	f = open("./%s/%s" % (Folder_Name,Analyze_File),'r')
	lines = f.readlines()
	f.close()
	Energy_Flag = False
	Nonbonded_Energy_List = []
	for line in lines:
		if len(line.strip().split()) >= 4:
			if line.strip().split()[0] == "Loop":
				Energy_Flag = False
			if Energy_Flag and int(line.strip().split()[0].strip()) >= Start_Point:
				Nonbonded_Energy_List.append(float(line.strip().split()[-1].strip()))
			if line.strip().split()[0] == "Step":
				Energy_Flag = True
		else:
			Energy_Flag = False
	Nonbonded_Energy_List = np.asarray(Nonbonded_Energy_List)
	Nonbonded_Energy = np.sum(Nonbonded_Energy_List)/len(Nonbonded_Energy_List)

	return Nonbonded_Energy

def Check_Finished(End_File,Folder_Name,Job_Name,Cluster_Login,Cluster_Location,End_Condition = "",Analyze_File = "",Shared_File_Location = ""):
	#Determines whether the output file already exists in the desired location, in the shared folder, or on the supercomputer
	Finished = False

	#Determine whether the file already exists locally
	if os.path.isfile("./%s/%s" % (Folder_Name,End_File)):
		#Check to see if local file is complete
		if End_Condition == "":
				print("Job Exists")
				Finished = True
		if End_Condition == "SPE_Orca":
			if Quantum_Orca("./%s/%s" % (Folder_Name,Analyze_File)):
				print("Job Exists")
				Finished = True
		if End_Condition == "Opt_Orca":
			if Opt_Orca("./%s/%s" % (Folder_Name,Analyze_File)):
				print("Job Exists")
				Finished = True
		if End_Condition == "SPE_QChem":
			if Quantum_QChem("./%s/%s" % (Folder_Name,Analyze_File)):
				print("Job Exists")
				Finished = True
		if End_Condition == "SPE_NWChem":
			if Quantum_NWChem("./%s/%s" % (Folder_Name,Analyze_File)):
				print("Job Exists")
				Finished = True

	#Determine whether the output file already exists in the shared folder, and if so copy it locally
	if not Finished and os.path.isfile("%s/%s" % (Shared_File_Location,End_File)):
		#Check whether file exists in shared folder
		if End_Condition == "":
				print("Job Exists")
				Finished = True
		if End_Condition == "SPE_Orca":
			if Quantum_Orca("%s/%s" % (Shared_File_Location,Analyze_File)):
				print("Job Exists")
				os.system("scp %s/%s ./%s" % (Shared_File_Location,End_File,Folder_Name))
				os.system("scp %s/%s ./%s" % (Shared_File_Location,Analyze_File,Folder_Name))
				Finished = True
		if End_Condition == "Opt_Orca":
			if Opt_Orca("./%s/%s" % (Folder_Name,Analyze_File)):
				print("Job Exists")
				Finished = True
		if End_Condition == "SPE_QChem":
			if Quantum_QChem("%s/%s" % (Shared_File_Location,Analyze_File)):
				print("Job Exists")
				os.system("scp %s/%s ./%s" % (Shared_File_Location,End_File,Folder_Name))
				os.system("scp %s/%s ./%s" % (Shared_File_Location,Analyze_File,Folder_Name))
				Finished = True
		if End_Condition == "SPE_NWChem":
			if Quantum_NWChem("%s/%s" % (Shared_File_Location,Analyze_File)):
				print("Job Exists")
				os.system("scp %s/%s ./%s" % (Shared_File_Location,End_File,Folder_Name))
				os.system("scp %s/%s ./%s" % (Shared_File_Location,Analyze_File,Folder_Name))
				Finished = True


	if not Finished:
		#Check whether file exists on cluster/supercomputer and copy it locally and to the shared folder
		os.system("scp %s:%s/%s ./%s" % (Cluster_Login,Cluster_Location,End_File,Folder_Name))
		if Analyze_File != End_File:
			os.system("scp %s:%s/%s ./%s" % (Cluster_Login,Cluster_Location,Analyze_File,Folder_Name))
		if os.path.isfile("./%s/%s" % (Folder_Name,End_File)):
			#Check to see if remote file is complete
			if End_Condition == "":
				print("Job Exists")
				if Shared_File_Location != "":
					os.system("scp ./%s/%s %s" % (Folder_Name,End_File,Shared_File_Location))
					if Analyze_File != End_File:
						os.system("scp ./%s/%s %s" % (Folder_Name,Analyze_File,Shared_File_Location))
				Finished = True
			if End_Condition == "Quantum_Orca":
				if Quantum_Orca("./%s/%s" % (Folder_Name,Analyze_File)):
					if Shared_File_Location != "":
						os.system("scp ./%s/%s %s" % (Folder_Name,End_File,Shared_File_Location))
						if Analyze_File != End_File:
							os.system("scp ./%s/%s %s" % (Folder_Name,Analyze_File,Shared_File_Location))
					print("Job Exists")
					Finished = True
			if End_Condition == "Opt_Orca":
				if Opt_Orca("./%s/%s" % (Folder_Name,Analyze_File)):
					print("Job Exists")
					Finished = True
			if End_Condition == "Quantum_QChem":
				if Quantum_QChem("./%s/%s" % (Folder_Name,Analyze_File)):
					if Shared_File_Location != "":
						os.system("scp ./%s/%s %s" % (Folder_Name,End_File,Shared_File_Location))
						if Analyze_File != End_File:
							os.system("scp ./%s/%s %s" % (Folder_Name,Analyze_File,Shared_File_Location))
					print("Job Exists")
					Finished = True
			if End_Condition == "Quantum_NWChem":
				if Quantum_NWChem("./%s/%s" % (Folder_Name,Analyze_File)):
					if Shared_File_Location != "":
						os.system("scp ./%s/%s %s" % (Folder_Name,End_File,Shared_File_Location))
						if Analyze_File != End_File:
							os.system("scp ./%s/%s %s" % (Folder_Name,Analyze_File,Shared_File_Location))
					print("Job Exists")
					Finished = True

	return Finished

def Submit_Job(Copy_File_List,Folder_Name,Submit_Script,End_File,Job_Name,Cluster_Login,Cluster_Location,Base_Cluster_Location,Scheduler_Type,End_Condition = "",Analyze_File = "",Shared_File_Location = "",Force_Copy = False,Force_New_Job = False):
	#Checks to see whether a job has already been completed, and if not submits a job to the cluster/supercomputer. Takes in: Copy_File_List: A List of files to be copied to the supercomputer to make the job run (usually a submit script, input file, and possible a data file), include path from ./ if in a subfolder; Folder_Name: The subdirectory to write all files to on the cluster/supercomputer and the shared folder; Submit_Script: script script file name (include path from ./ if in a subfolder); End_File: The script that will be output at the end of the run (for MD run, usually a data file with the final positions; for quantum file, the output file); Job_Name: the local location to copy files to. Will be made as a subdirectory under ./; Cluster_Login: String specifying your login for the cluster/supercomputer (e.g. yourname@comet.sdsc.edu); Cluster_Location: String specifying the folder to copy files to/from on the cluster; Base_Cluster_Location: String specifying cluster location for project as a whole; End_Condition: Optional string indicating if there is a special end condition that must be met or else resubmit the job. Options include "Quantum_Orca","Quantum_QChem", and "Quantum_NWChem" (checks that calculation converged and terminated for each quantum chemistry packaged); Analyze_File: The output file that will be parsed to get relevant information, required if End_Condition is defined; Shared_File_Location: String indicating a location where files can be shared between different runs of the program. Particularly useful for quantum calculations to avoid rerunning calculations between identical rings; Force_Copy: Boolean variable that if True forces the program to recopy from the cluster/supercomputer even if the file already exists locally; Force_New_Job: Boolean variable that if True forces the program to recopy from the cluster/supercomputer even if the job has already been run

	if End_Condition != "" and Analyze_File == "":
		raise Exception("Analyze_File must be defined if End_Condition is defined")
	#Make local directories
	os.system('mkdir ./%s' % Folder_Name)
	if Shared_File_Location != "":
		os.system('mkdir %s/%s' % (Shared_File_Location,Folder_Name))
		os.system("rm -f %s/%s.*" % (Shared_File_Location,Job_Name))

	if Force_New_Job or Force_Copy:
		os.system('rm -f ./%s/%s' % (Folder_Name,End_File))
		os.system('rm -f ./%s/%s' % (Folder_Name,Analyze_File))
		if Shared_File_Location != "":
			os.system("rm -f %s/%s" % (Shared_File_Location,End_File))
			os.system("rm -f %s/%s" % (Shared_File_Location,Analyze_File))
			os.system("rm -f %s/%s.*" % (Shared_File_Location,Job_Name))
		if Force_New_Job:
			subprocess.call(["ssh", "%s" % Cluster_Login, "rm -f %s/%s" % (Cluster_Location,End_File)])
			subprocess.call(["ssh", "%s" % Cluster_Login, "rm -f %s/%s" % (Cluster_Location,Analyze_File)])
			subprocess.call(["ssh", "%s" % Cluster_Login, "rm -f %s/%s.*" % (Cluster_Location,Job_Name)])

	if not Force_New_Job:
		Finished = Check_Finished(End_File,Folder_Name,Job_Name,Cluster_Login,Cluster_Location,End_Condition,Analyze_File,Shared_File_Location)
	else:
		Finished = False
	#Check to see whether the file already exists locally

	#run job
	if not Finished:
		#Make all subdirectories for location specified on cluster/supercomputer
		Continuing_Cluster_Location = Base_Cluster_Location
		Sub_Dir_List = Cluster_Location.split(Base_Cluster_Location)[-1].split('/')
		for sub_dir in Sub_Dir_List:
			Continuing_Cluster_Location = Continuing_Cluster_Location + ("%s/" % sub_dir)
			subprocess.call(["ssh", "%s" % Cluster_Login, "mkdir %s" % Continuing_Cluster_Location])
		#Copy files to supercomputer and run
		for file in Copy_File_List:
			os.system("scp %s %s:%s" % (file,Cluster_Login,Cluster_Location))
		if Scheduler_Type == "SLURM":
			subprocess.call(["ssh", "%s" % Cluster_Login, "sbatch %s/%s" % (Cluster_Location,Submit_Script)])
		elif Scheduler_Type == "TORQUE":
			subprocess.call(["ssh", "%s" % Cluster_Login, "qsub %s/%s" % (Cluster_Location,Submit_Script)])


def Return_Info(Analyze_File,End_File,Folder_Name,Job_Name,Cluster_Login,Cluster_Location,End_Condition = "",Shared_File_Location = "",Return_Energy_Orca = False,Return_Loewdin_Charges_Orca = False,Return_Energy_QChem = False,Return_NBO_Charges_QChem=False,Return_Energy_NWChem = False,Return_Averaged_Energy_Lammps = False,Start_Point = 0):
	#Checks for completed runs from the cluster/supercomputer and returns desired values. Takes in: Analyze_File: The file which will be used to return values from (may be the same as End_File); End_File: The script that will be output at the end of the run (for MD run, usually a data file with the final positions; for quantum file, the output file); Folder_Name: the local location to copy files to. Will be made as a subdirectory under ./; Cluster_Login: String specifying your login for the cluster/supercomputer (e.g. yourname@comet.sdsc.edu); Cluster_Location: String specifying the folder to copy files to/from on the cluster; End_Condition: Optional string indicating if there is a special end condition that must be met or else resubmit the job. Options include "Quantum_Orca","Quantum_QChem", and "Quantum_NWChem" (checks that calculation converged and terminated for each quantum chemistry packaged); Shared_File_Location: String indicating a location where files can be shared between different runs of the program. Particularly useful for quantum calculations to avoid rerunning calculations between identical rings; The Return flags each are boolean values which, if True, make the function return the requested value for the requested program (e.g. Return_Energy_Orca will return the single point energy from an Orca output file). Values are always returned in a tuple in the same order as the flags are listed above.

	Finished = Check_Finished(End_File,Folder_Name,Job_Name,Cluster_Login,Cluster_Location,End_Condition,Analyze_File,Shared_File_Location)
	#check if file is in shared folder
	if not Finished:
		#loop until job completes
		while not Finished:
			#Check whether file exists on cluster/supercomputer and copy it locally
			os.system("scp %s:%s/%s ./%s" % (Cluster_Login,Cluster_Location,End_File,Folder_Name))
			if os.path.isfile("./%s/%s" % (Folder_Name,End_File)):
				#Check to see if remote file is complete
				if End_Condition == "Quantum_Orca":
					if Quantum_Orca("./%s/%s" % (Folder_Name,End_File)):
						print("Job Exists")
						Finished = True
				if End_Condition == "Quantum_QChem":
					if Quantum_QChem("./%s/%s" % (Folder_Name,End_File)):
						print("Job Exists")
						Finished = True
				if End_Condition == "Quantum_NWChem":
					if Quantum_NWChem("./%s/%s" % (Folder_Name,End_File)):
						print("Job Exists")
						Finished = True
			if not Finished:
				time.sleep(300)

		#copy file to be analyzed to the local folder and the shared folder
		os.system("scp %s:%s/%s ./%s" % (Cluster_Login,Cluster_Location,Analyze_File,Folder_Name))
		if Shared_File_Location != "":
			os.system("scp %s:%s/%s %s/%s" % (Cluster_Login,Cluster_Location,Analyze_File,Shared_File_Location,Folder_Name))
	
	Return_Values = []
	#parse the analyze_file for values
	if Return_Energy_Orca:
		Return_Values.append(Energy_Orca(Analyze_File,Folder_Name))
	if Return_Loewdin_Charges_Orca:
		Return_Values.append(Loewdin_Charges_Orca(Analyze_File,Folder_Name))
	if Return_Energy_QChem:
		Return_Values.append(Energy_Qchem(Analyze_File,Folder_Name))
	if Return_NBO_Charges_QChem:
		Return_Values.append(NBO_Charges_QChem(Analyze_File,Folder_Name))
	if Return_Energy_NWChem:
		Return_Values.append(Energy_NWChem(Analyze_File,Folder_Name))
	if Return_Averaged_Energy_Lammps:
		Return_Values.append(Energy_Averaged_LAMMPS(Analyze_File,Folder_Name,Start_Point))
	Return_Values = tuple(Return_Values)
	return Return_Values
