#!usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import time
import pickle
import matplotlib as mpl
from pylab import rcParams
import seaborn as sns
#sns.set_style('ticks')
from Analyze_Traj import Output
#from Analyze_Traj import RDF_OBJ
#

mpl.rcParams['axes.color_cycle'] = ['b', 'k','r',  'c', 'y', 'm']

class RDF_OBJ:
    def __init__(self, Type1, Type2):
        self.Radius = np.zeros(182,dtype=float)
        self.RDF = np.zeros(182, dtype=float)
        self.Num_Snaps = 0.0
        self.Type1 = Type1
        self.Type2 = Type2
    def Add_Snap(self, RDF, radius):
        self.RDF += RDF
        self.Radius = radius
        self.Num_Snaps += 1.0
    def Normalize(self):
        self.RDF = self.RDF/self.Num_Snaps





def Plot_Dist(Pickle_File, Bins):
    File = open(Pickle_File,'rb')
    Output = pickle.load(File)
    plt.hist(Output.C_Ratio, bins = Bins , normed=False, histtype = 'bar')
    plt.legend( loc = 'upper right', frameon = False, fontsize= 40)
    plt.tick_params( labelsize = 30, width=4, length=10)
    plt.xlabel("Characteristic Ratio", fontsize=25)
    plt.ylabel("Fraction", fontsize=25)
    plt.show()

def Plot_Conj(Pickle_File, Bins):
    File = open(Pickle_File,'rb')
    Output = pickle.load(File)
    print Output.Conjugation_Dist
    Conjugation_Dist, bin_edges = np.histogram(Output.Conjugation_Dist, bins= (np.max(Output.Conjugation_Dist)-1), normed=True)
    plt.yscale('log')
    #plt.hist(Output.Conjugation_Dist, bins = Bins , normed=True, histtype = 'bar')
    plt.plot(bin_edges[0:-1], Conjugation_Dist, marker='o', fillstyle = 'right')
    plt.legend( loc = 'upper right', frameon = False, fontsize= 25)
    plt.tick_params( labelsize = 20, width=2, length=7)
    plt.xlabel("Conjugation Length", fontsize=25)
    plt.ylabel("Probability", fontsize=25)
    plt.show()
    return


def Plot_RDF(Pickle_File):
    File = open(Pickle_File,'rb')
    Output = pickle.load(File)
    for RDF1 in Output.RDF:
        print len(RDF1.Radius[0:-1]), len(RDF1.RDF)
        RDF1.RDF[0:10] = 0.0
        plt.plot(RDF1.Radius[0:-1], RDF1.RDF, label = RDF1.Type1 + "-" + RDF1.Type2, linewidth=5)
    plt.legend( loc = 'lower right', frameon = False, fontsize= 40)
    plt.tick_params( labelsize = 30, width=4, length=10)
    plt.xlabel("Distance ($\AA$)", fontsize=40)
    plt.ylabel("Pair Distribution Function", fontsize=40)
    plt.axvline(x=3.9, linewidth = 20, alpha= 0.4, color = 'g')
    plt.xlim((0,25))
    plt.ylim((0,1.2))
    X_Max = 25
    Y_Max = 1.2
    plt.axhline(y=Y_Max,linewidth=4, color='k'); plt.axhline(linewidth=4, color='k'); plt.axvline(linewidth=4, color='k');plt.axvline( x=X_Max,linewidth=4, color='k');
    plt.show()


def Plot_Tangent_Correlation(Pickle_File):
    File = open(Pickle_File,'rb')
    Output = pickle.load(File)
    plt.plot(Output.Tangent_Correlation, linewidth=3)
    plt.ylim((-.2, 1.0))
    plt.xlim ((0, 24))
    plt.axhline(y=-.2,linewidth=4, color='k');
    plt.axhline(y=1.0,linewidth=4, color='k');
    plt.axvline(linewidth=4, color='k');
    plt.axvline( x=24 ,linewidth=4, color='k');
    plt.tick_params( labelsize = 30, width=2, length=7)
    plt.ylabel("Tangent Correlation Function", fontsize=25)
    plt.xlabel("Distance along chain (N)", fontsize=25)
    plt.show()

def Plot_Binormal_Correlation(Pickle_File):
    File = open(Pickle_File,'rb')
    Output = pickle.load(File)
    plt.plot(Output.Binormal_Correlation, linewidth=3)
    print Output.Binormal_Correlation
    plt.ylim((-0.1, 1.0))
    plt.xlim ((0, 24))
    plt.axhline(y=-.2,linewidth=4, color='k');
    plt.axhline(y=1.0,linewidth=4, color='k');
    plt.axvline(linewidth=4, color='k');
    plt.axvline( x=24 ,linewidth=4, color='k');
    plt.tick_params( labelsize = 30, width=2, length=7)
    plt.ylabel("Binormal Correlation Function", fontsize=25)
    plt.xlabel("Distance along chain (N)", fontsize=25)
    plt.show()



def Plot_Structure(Pickle_File):
    File = open(Pickle_File,'rb')
    Output = pickle.load(File)
    for RDF1 in Output.RDF:
        Q = 1.0/RDF1.Radius
        Struct = np.fft.fft(np.ones(len(RDF1.RDF), dtype=float) - RDF1.RDF)
        plt.plot(Q, Struct, label = RDF1.Type1 + "-" + RDF1.Type2, linewidth=5)
    plt.legend( loc = 'upper left', frameon = False, fontsize= 25)
    X_Max = 0.5
    Y_Max = 2
    plt.ylim((0,Y_Max))
    plt.xlim((.0667,X_Max))
    plt.tick_params( labelsize = 20, width=2, length=7)
    plt.xlabel("Distance ($\AA$)", fontsize=25)
    plt.ylabel("Pair Distribution Function", fontsize=25)
    plt.show()

def Plot_Ramachadran(Pickle_File):
    File = open(Pickle_File, 'rb')
    Output = pickle.load(File)
    Bin = 10
    Rama = np.zeros([360/Bin+1,360/Bin+1  ], dtype=float)
    X = np.asarray(range(-180,180+Bin, Bin))
    Y = np.asarray(range(-180,180+Bin, Bin))
    print Rama
    #print Output.Dihedral[0]

    for i in range(len(Output.Dihedral[0])):
            dih_1 = int(Output.Dihedral[0][i]/Bin) + 180/Bin
            dih_2 = int(Output.Dihedral[1][i]/Bin) + 180/Bin
            print dih_2
            Rama[ dih_1, dih_2] += 1

    Rama = Rama / float(len(Output.Dihedral[0]))
    cp = plt.contourf(X, Y, Rama, cmap='hot')
    cbar = plt.colorbar(cp,format='%.0e')
    fsize = 30
    cbar.set_label("Probability", size=fsize)
    cbar.ax.tick_params(labelsize=fsize)
    plt.xticks([-180, -135, -90, -45, 0, 45, 90, 135, 180])
    plt.yticks([-180, -135, -90, -45, 0, 45, 90, 135, 180])
    plt.tick_params( labelsize = fsize, width=2, length=7)
    plt.xlabel('A-D Dihedral Angle (degrees)', fontsize=fsize)
    plt.ylabel('D-A Dihedral Angle (degrees)', fontsize=fsize)
    #plt.tick_params( width=2, length=7)
    plt.show()

    return

def Plot_Dihedral(Pickle_File):
    File = open(Pickle_File, 'rb')
    Output = pickle.load(File)
    plt.axvspan(-180, -140, color='g', alpha=0.5)
    plt.axvspan(140, 180, color='g', alpha=0.5)
    plt.axvspan(-40, 40, color='g', alpha=0.5)
    plt.hist([Output.Dihedral[1], Output.Dihedral[0]], bins=90, normed=1, histtype='bar', label=["DA", "AD"], color=['k', 'r'])
    plt.legend( loc = 'upper left', frameon = False, fontsize= 40)
    plt.tick_params( labelsize = 45, width=2, length=7)
    plt.xlabel("Angle (Degrees)", fontsize=40)
    plt.ylabel("Probability", fontsize=40)
    plt.ylim(0,.008)
    plt.xlim((-180,180))
    plt.yticks([0, 0.002, 0.004, 0.006, 0.008])
    plt.xticks([-180, -135, -90, -45, 0, 45, 90, 135, 180])
    plt.axvline(x=-180,linewidth=4, color='k')
    plt.axvline( x=180, linewidth=4, color='k')
    plt.axhline(y=.008, linewidth=4, color='k')
    plt.axhline(linewidth=4, color='k')
    print Output.Dihedral
    syn = 0
    anti = 0
    nonplanar = 0
    for i in range(len(Output.Dihedral[0])):
        Angle = Output.Dihedral[1][i]
        if Angle <= -140 or Angle >= 140:
            anti += 1
        elif Angle >= -40 and Angle <= 40:
            syn += 1
        else:
            nonplanar += 1
    Total = syn + anti + nonplanar
    syn = float(syn)/float(Total)
    anti = float(anti)/float(Total)
    nonplanar = float(nonplanar)/float(Total)
    print "syn", syn
    print "Anti", anti
    print "Nonplanar", nonplanar

    
    plt.show()

#Plot_Binormal_Correlation("Output_600_BN.pickle")
Plot_Ramachadran("Output_600_BN.pickle")
#Plot_Dihedral("Output_420.pickle")
#Plot_RDF("Output_420.pickle")
#Plot_Tangent_Correlation("Output_420.pickle")
