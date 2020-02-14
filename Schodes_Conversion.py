import sys

def main():
    Script, FileName = sys.argv
    AtomWeights = {'H':'1.008','C':'12.011','N':'14.007','O':'15.999','F':'18.998','S':'32.065'}
    f = open(FileName,'r')
    lines = f.readlines()
    BondFlag = False
    AngleFlag = False
    DihedralFlag = False
    ImproperFlag = False
    NonbondedFlag = False
    TypeFlag = False
    BondParams = []
    BondDefs = []
    BondNum = 1
    AngleDefs = []
    AngleParams = []
    AngleNum = 1
    DihDefs = []
    DihParams = []
    DihNum = 1
    NonDefs = []
    NonParams = []
    NonNum = 1
    ImDefs = []
    ImParams = []
    ImNum = 1
    Masses = []
    TypeNum = 1
    for line in lines:
        if line.strip().split() != [] and line.strip().split()[0] == "BOND":
            BondFlag = True
        if BondFlag:
            if line.strip().split() != [] and line.strip().split()[0] == "ANGLES,":
                BondFlag = False
                AngleFlag = True
            elif line.strip().split() != [] and line.strip().split()[0].isdigit() :
                BondParams.append(str(BondNum)+" "+line.strip().split()[2]+" "+line.strip().split()[3])
                BondDefs.append(str(BondNum)+" "+line.strip().split()[0]+" "+line.strip().split()[1])
                BondNum+=1

        if AngleFlag:
            if line.strip().split() != [] and line.strip().split()[0] == "DIHEDRAL":
                AngleFlag = False
                DihedralFlag = True
            elif line.strip().split() != [] and line.strip().split()[0].isdigit():
                AngleParams.append(str(AngleNum)+" "+line.strip().split()[3]+" "+line.strip().split()[5])
                AngleDefs.append(str(AngleNum)+" "+line.strip().split()[0]+" "+line.strip().split()[1]+" "+line.strip().split()[2])
                AngleNum+=1
        if DihedralFlag:
            if len(line.strip().split()) >1 and line.strip().split()[1] == "Improper":
                DihedralFlag = False
                ImproperFlag = True
            elif line.strip().split() != [] and line.strip().split()[0].isdigit():
                DihParams.append(str(DihNum)+" "+line.strip().split()[4]+" "+line.strip().split()[5]+line.strip().split()[6]+" "+line.strip().split()[7])
                DihDefs.append(str(DihNum)+" "+line.strip().split()[0]+" "+line.strip().split()[1]+" "+line.strip().split()[2]+" "+line.strip().split()[3])
                DihNum+=1
        if ImproperFlag:
            if line.strip().split() != [] and line.strip().split()[0] == "NONBONDED":
                ImproperFlag = False
                NonbondedFlag = True
            elif line.strip().split() != [] and line.strip().split()[0].isdigit():
                ImParams.append(str(ImNum)+" "+line.strip().split()[2]+" "+line.strip().split()[3])
                ImDefs.append(str(ImNum)+" "+line.strip().split()[0]+" "+line.strip().split()[1])
                ImNum+=1
        if NonbondedFlag:
            if line.strip().split() != [] and line.strip().split()[0] == "TIME:":
                NonbondedFlag = False
            elif len(line.strip().split()) >1 and line.strip().split()[0].isdigit():
                print line.strip().split()[0]
                NonParams.append(str(NonNum)+" "+line.strip().split()[2]+" "+line.strip().split()[3])
                NonDefs.append(str(NonNum)+" "+line.strip().split()[0]+" "+line.strip().split()[1])
                NonNum+=1
        if line.strip().split() != [] and line.strip().split()[0] == "Connection":
            TypeFlag = True
            print line
        if TypeFlag:
            if len(line.strip().split()) > 1 and line.strip().split()[0][1].isdigit():
                Masses.append(str(TypeNum)+" "+AtomWeights[line.strip().split()[0][0]])
                TypeNum+=1
    f.close()
    f = open('Test.txt','w')
    f.write("Masses\n\n")
    f.write("\n".join(Masses))
    f.write("\n\nPair Coeffs #LJ\n\n")
    f.write("\n".join(NonParams))
    f.write("\n\nBond Coeffs #harmonic\n\n")
    f.write("\n".join(BondParams))
    f.write("\n\nAngle Coeffs #harmonic\n\n")
    f.write("\n".join(AngleParams))
    f.write("\n\nDihedral Coeffs #opls\n\n")
    f.write("\n".join(DihParams))
    f.write("\n\nImproper Coeffs #harmonic\n\n")
    f.write("\n".join(ImParams))



if __name__=='__main__': main()
