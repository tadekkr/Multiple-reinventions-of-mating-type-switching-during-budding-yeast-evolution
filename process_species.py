import sys
import os
import glob
import subprocess
import math

def samespec(s1, s2):
    if (s1.lower() in s2.lower()) or (s2.lower() in s1.lower()):
        return True
    else:
        return False

def main(args):
    clades = {}
    with open("clades.txt") as clads:
        k = -1
        for line in clads:
            k = k + 1
            tabs = line.strip().split("\t")
            if k >= 1:
                clades[tabs[1]] = tabs[3]
    l1 = os.listdir("./y1000_genomes/300_genomes/")
    l2 = []
    for ent in l1:
        if ent[-4:] == ".fas" and ent not in l2 and "ania_huma" not in ent:
            l2.append(ent.split(".")[0])
    l3 = []
    cladesofl2 = {}
    ofclades = {}
    for ent in l2:
        found = 0
        for ent2 in clades.keys():
            if samespec(ent2, ent):
                cladesofl2[ent] = clades[ent2]
                if clades[ent2] not in ofclades.keys():
                    ofclades[clades[ent2]] = [ent]
                else:
                    ofclades[clades[ent2]].append(ent)
                found = 1
        #if found != 1:
        #    print ent
    for ent in l2:
        l3.append([len(ent.split("_")), ent])
    l3.sort()
    l4 = []
    ofclades2 = {}
    newtoold = {}
    oldtonew = {}
    for ent in l3:
        num = ent[0]
        nam = ent[1]
        clad = cladesofl2[nam]
        if num == 4:
            nam2 = nam.split("_")
            nnam = "_".join([nam2[1], nam2[2], nam2[3], nam2[0]])
            newtoold[nnam] = nam
            oldtonew[nam] = nnam
            if clad not in ofclades2.keys():
                ofclades2[clad] = [nnam]
            else:
                ofclades2[clad].append(nnam)
        else:
            nam2 = nam.split("_")
            nnam = nam
            newtoold[nnam] = nam
            oldtonew[nam] = nnam
            if clad not in ofclades2.keys():
                ofclades2[clad] = [nnam]
            else:
                ofclades2[clad].append(nnam)
    for clad in ofclades2.keys():
        ofclades2[clad].sort()
    with open("addition.txt", "w") as ad:
        for clad in ofclades2.keys():
            clad2 = clad.split("/")[0]
            with open(clad2 + ".txtx", "w") as cladfil:
                for spec in ofclades2[clad]:
                    filname = "Y1000/" + newtoold[spec] + ".fas"
                    cladfil.write(filname + "\n")
                    ad.write(filname + "\t" + clad + " " + spec + "\n")
                ad.write("\n\n")



    #print preparegenestowrite
if __name__ == "__main__":
    sys.exit(main(sys.argv))
