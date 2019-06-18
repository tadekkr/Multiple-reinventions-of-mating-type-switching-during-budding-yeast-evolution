import sys
import os
import glob
import subprocess
import math

def findgennam(currentnam):
    nam = currentnam
    #print nam
    found = 0
    found2 = 0
    saccer = "---"
    anc = "---"
    with open("pillars.tab", "r") as fil:
        for line in fil:
            if nam in line:
                lin = line.strip().split("\t")
                #print lin
                saccer = lin[11]
                anc = lin[12]
                #print anc
                if saccer == "---":
                    pass
                else:
                    found = 1
                break
    if found == 1:
        with open("association.tab", "r") as fil:
            for line in fil:
                if saccer in line:
                    lin = line.strip().split("\t")
                    saccer = lin[2]
                    found2 = 1
                    break
    if saccer != "---":
        if anc != "---":
            gennam = saccer + "_(" + anc + ")"
        else:
            gennam = saccer
    else:
        if anc != "---":
            gennam = nam + "_(" + anc + ")"
        else:
            gennam = nam
    #print gennam
    return gennam

def inrange(num, rang):
    if num >= rang[0] and num <= rang[1]:
        return True
    else:
        return False

def findorfs(seq, stpsfor, stpsbck):
    #print "FINDING ORFS"
    l = len(seq)
    forframestps = {}
    forframestps[0] = []
    forframestps[1] = []
    forframestps[2] = []
    for frame in [0, 1, 2]:
        k = -3
        while k + frame + 3 <= l:
            k = k + 3
            cod = seq[(k + frame):(k + frame + 3)].upper()
            if cod in stpsfor:
                forframestps[frame].append(k + frame)
    bckframestps = {}
    bckframestps[0] = []
    bckframestps[1] = []
    bckframestps[2] = []
    for frame in [0, 1, 2]:
        k = -3
        while k + frame + 3 <= l:
            k = k + 3
            cod = seq[(k + frame):(k + frame + 3)].upper()
            if cod in stpsbck:
                bckframestps[frame].append(k + frame)
    minlen = 500
    answs = []
    for frame in [0, 1, 2]:
        if len(forframestps[frame]) > 0:
            if forframestps[frame][0] > minlen:
                answs.append([[frame, forframestps[frame][0] + 3], 0, seq[frame:(forframestps[frame][0] + 3)]])
            if forframestps[frame][-1] < l - minlen:
                answs.append([[forframestps[frame][-1] + 3, l - ((l - frame) % 3)], 0, seq[(forframestps[frame][-1] + 3):(l - ((l - frame) % 3))]])
            l2 = len(forframestps[frame])
            for i in range(l2 - 1):
                if forframestps[frame][i + 1] - forframestps[frame][i] > minlen:
                    answs.append([[forframestps[frame][i] + 3, forframestps[frame][i + 1] + 3], 0, seq[(forframestps[frame][i] + 3):(forframestps[frame][i + 1] + 3)]])
        if len(bckframestps[frame]) > 0:
            if bckframestps[frame][0] > minlen:
                answs.append([[frame, bckframestps[frame][0]], 1, seq[frame:(bckframestps[frame][0])]])
            if bckframestps[frame][-1] < l - minlen:
                answs.append([[bckframestps[frame][-1], l - ((l - frame) % 3)], 1, seq[(bckframestps[frame][-1]):(l - ((l - frame) % 3))]])
            l2 = len(bckframestps[frame])
            for i in range(l2 - 1):
                if bckframestps[frame][i + 1] - bckframestps[frame][i] > minlen:
                    answs.append([[bckframestps[frame][i], bckframestps[frame][i + 1]], 1, seq[(bckframestps[frame][i]):(bckframestps[frame][i + 1])]])
    #print seq
    #print forframestps
    #print bckframestps
    #print answs
    return answs



def findindex(scafind, ordofscaf):
    k = -1
    for i in ordofscaf:
        k = k + 1
        if scafind in i:
            return k
    return -1

def findnam(scafind, ordofscaf):
    k = -1
    for i in ordofscaf:
        k = k + 1
        if scafind in i:
            return i
    return -1

def extractseqaround(gen, lengg, scafs, ordof):
    howfar = int(lengg / 2)
    scaf = gen[1]
    locs = gen[6]
    mid = locs[0] + int((locs[1] - locs[0])/2)
    bigseq = scafs[findnam(scaf, ordof)]
    stt = max(mid - howfar, 0)
    enn = min(mid + howfar, len(bigseq))
    smallseq = bigseq[(stt):(enn)]
    return [stt, enn, scaf, smallseq]

def clearannotations(fil):
    lins = []
    with open(fil, "r") as f:
        for line in f:
            lins.append(line)
    write = 0
    with open(fil, "w") as f:
        for l in lins:
            if ">" in l:
                write = 1
            if write == 1:
                f.write(l)

def writegenes(writinglist):
    linestowrite = []
    for ent in writinglist:
        typp = ent[0]
        scaf = ent[1]
        st = ent[2]
        en = ent[3]
        spec = ent[4]
        mode = ent[5]
        if mode == 1:
            en = en - 1
            endlin1 = "complement(" + str(st) + ".." + str(en) + ")\n"
        else:
            endlin1 = str(st) + ".." + str(en) + "\n"
        endlin2 = typp + '"\n'
        lin1 = "FT   CDS             "
        lin2 = 'FT                   /gene="'
        wri = lin1 + endlin1 + lin2 + endlin2
        linestowrite.append(wri)
    fil = []
    wher = "../editing/" + spec.split(".")[0] + ".fas"
    clearannotations(wher)
    with open(wher, "r") as f:
        for lin in f:
            fil.append(lin)
    with open(wher, "w") as f:
        for i in linestowrite:
            f.write(i)
        for j in fil:
            f.write(j)

def writemisc(misclist):
    linestowrite = []
    for ent in misclist:
        typp = ent[0]
        scaf = ent[1]
        st = ent[2]
        en = ent[3]
        spec = ent[4]
        mode = ent[5]
        if mode == 1:
            en = en - 1
            endlin1 = "complement(" + str(st) + ".." + str(en) + ")\n"
        else:
            endlin1 = str(st) + ".." + str(en) + "\n"
        endlin2 = typp + '"\n'
        lin1 = "FT   misc_feature    "
        lin2 = 'FT                   /gene="'
        wri = lin1 + endlin1 + lin2 + endlin2
        linestowrite.append(wri)
    fil = []
    wher = "../editing/" + spec.split(".")[0] + ".fas"
    with open(wher, "r") as f:
        for lin in f:
            fil.append(lin)
    with open(wher, "w") as f:
        for i in linestowrite:
            f.write(i)
        for j in fil:
            f.write(j)


def findlongest(sequence, position, backward, stps):
    savedbeg = 0
    savedend = 0
    savedlnght = 0
    for start in [position, position - 1, position - 2]:
        jazda = start
        cod = sequence[start:(start + 3)]
        while cod not in stps:
            jazda = jazda + 3
            if jazda >= len(sequence):
                jazda = jazda - 3
                break
            cod = sequence[jazda:(jazda + 3)].upper()
        end = jazda
        jazda = start
        cod = sequence[start:(start + 3)]
        while cod not in stps:
            jazda = jazda - 3
            if jazda < 0:
                jazda = jazda + 3
                break
            cod = sequence[jazda:(jazda + 3)].upper()
        beg = jazda
        if (end - beg) > savedlnght:
            savedbeg = beg
            savedend = end
            savedlnght = end - beg
    return (savedbeg, savedend)
def main(args):
    blastingmode = 0
    typy = ["a1", "a2", "alpha1", "alpha2"]
    fakesknown = ['HPODL_1823', 'OPOL_17053', 'PKUD0AY00100', 'AFR301C', 'AFL049C', 'ZYRO0B08998g', 'OPOL_15913', 'CAGL0C02563g', 'CAGL0H02959g', 'CAGL0C01551g', 'NABA0s09e00792g', 'KAFR0G03650', 'CJAD0D00390', 'CJAD0D04350', 'NDAI0A04660', 'Ecym_2247', 'Ecym_3378', 'KLLA0D10043g', 'Kpol_1004.50', 'SU7_0508', 'KNAG0C02770', 'KMAR0S00220', 'Kwal_33.12977', 'KLTH0D02222g', 'Kwal_27.11092', 'Skud_4.640', 'TPHA0E02340', 'SAKL0C03652g', 'SAKL0A05610g', 'KLTH0F09636g', 'NDAI0K00400', 'KLTH0F00352g', 'KLTH0H05236g', 'TDEL0D00460', 'BN7_4777', 'CABR0s05e01705g', 'CACA0s15e03630g', 'CAGL0E00363g', 'NADE0s27e03487g', 'NADE0s27e04466g', 'CANI0s17e03619g', 'CAGL0G09295g', 'Ptan_3805', 'NCAS0H01130', 'NDAI0F02170', 'KAFR0E02870', 'SU7_3505', 'YPL177C', 'Suva_16.132', 'Skud_16.103', 'Smik_6.374', 'Smik_7.181', 'YGL096W', 'CJAD0A08660', 'BN7_141', 'TBLA0C03740', 'PP7435_Chr1-0672', 'TPHA0B02580', 'CABR0s29e02662g', 'TDEL0G02680', 'TDEL0F01810', 'Kpol_1036.56', 'CAGL0L07436g', 'BN7_5448', 'TBLA0B07440', 'KMAR0H01540', 'ZYRO0G22044g', 'Ptan_1165', 'KAFR0A04370', 'KNAG0M00500', 'KNAG0F02210', 'YLR356W', 'CJAD0E01080', 'CABR0s11e00836g', 'ZYBA0C04566g', 'OPOL_103301', 'Ptan_3465', 'CJAD0B04080', 'CAGL0M01914g', 'TBLA0A01840', 'Ptan_1381', 'NABA0s08e02651g', 'KLTH0F00968g', 'CACA0s12e08217g', 'CABR0s16e04136g', 'TDEL0A04320', 'TDEL0C05830', 'ZYBA0S00100g', 'ZYBA0N00100g', 'TPHA0E03630', 'TDEL0F04400', 'PKUD0BJ00230', 'Skud_10.240', 'OPOL_17154', 'SAKL0G05896g', 'NDAI0J02520', 'HPODL_1214', 'Kpol_1072.26', 'KNAG0I00830', 'CJAD0B06450', 'KLTH0E15642g', 'TBLA0B05050', 'Ptan_1048', 'CJAD0G02800', 'Ptan_3946', 'KLTH0H05984g', 'PP7435_Chr3-0616', 'pusto']
    known = {}
    known["a1"] = ["Smik_3.143", "Suva_3.77", "KAFR0G00180", "KNAG0C00795",
                    "NCAS0B08560", "NDAI0A00670", "TPHA0E03595",
                    "Kpol_2002.48", "TDEL0C05820", "KLLA0C03135g",
                    "AFR643C", "Sklu_YGOB_MATa1", "Kthe_YGOB_MATa1",
                    "BN7_1723", "CJAD0D03400", "PKUD0PD00100",
                    "Ptan_5337", "HPODL_1930", "TDEL0G00190",
                    "ADL394C", "SAKL0C03674g", "SAKL0C03674g",
                    "YCR097W", "CBS432_03G01410", "Suva_3.138",
                    "CANI0s19e05281g", "NADE0s06e02733g", "CABR0s40e02652g",
                    "CAGL0E00341g", "KNAG0C01625", "NCAS0A15350",
                    "TBLA0G01060", "NDAI0H00100", "Kpol_1058.2",
                    "ZYRO0C18348g", "KLLA0B14553g", "KLLA0B14553g",
                    "KLTH0F00242g", "Kwal_YGOB_HMa1", "AER456W",
                    "Skud_80.2", "PP7435_Chr4-0002", "Kwal_33.12991",
                    "CANI0s28e02751g", "CACA0s33e04911g", "HPODL_1923",
                    "OPOL_93623", "Smik_88.2", "Smik_88.1", "TPHA0E04055"]
    known["a2"] = ["TDEL0C05810", "KLLA0C03157g", "AFR643WA", "SAKL0C03696g",
                    "KLTH0F03278g", "BN7_1725", "CJAD0D03390", "BN7_1724",
                    "ZYRO0C18326g", "KLLA0B14575g", "AER455C", "KLTH0F00264g",
                    "Kwal_YGOB_HMa2", "TDEL0G00200", "ADL393W", "KMAR0S00210",
                    "Suva_3.135", "TDEL0E00340"]
    known["alpha1"] = ["YCR040W", "Smik_3.144", "Suva_3.79", "SU7_0315",
                        "CAGL0B01243g", "NABA0s02e00641g", "NABA0s02e00103g",
                        "KAFR0D00710", "TBLA0A07040", "TPHA0E03620",
                        "ZYRO0F15840g", "ZYBA0N00122g", "KMAR0E01330",
                        "Ecym_1114", "Kwal_YGOB_matalpha1", "Ptan_5336",
                        "PAS_chr4_1001", "PP7435_Chr4-0073", "YCL066W",
                        "CBS432_03G00050", "Suva_3.148", "CANI0s02e03669g",
                        "CANI0s02e03229g", "NADE0s06e00103g", "CABR0s26e00135g",
                        "CABR0s26e00612g", "CAGL0B00242g", "CACA0s21e03798g",
                        "CACA0s21e03248g", "KNAG0C00150", "NCAS0B09150",
                        "NDAI0A00100", "TBLA0A07590", "TPHA0E04080",
                        "Kpol_2002.2", "ZYRO0F18590g", "TDEL0C07010",
                        "KLLA0C00352g", "Ecym_1003", "KLTH0F00374g",
                        "Kwal_YGOB_HMalpha1", "CBS432_03G00950",
                        "Skud_3.119", "CAGL0B01243g", "PKUD0OO00110",
                        "Skud_71.1", "Skud_80.1", "Kwal_33.matalpha1",
                        "Suva_69.2", "Smik_88.1"]
    known["alpha2"] = ["YCR039C", "CBS432_03G00940", "Smik_3.142", "Skud_3.118",
                        "Suva_3.75", "CANI0s02e03226g", "CANI0s02e03666g",
                        "CABR0s26e00138g", "CABR0s26e00619g", "CAGL0B01265g",
                        "CACA0s21e03802g", "CACA0s21e03251g", "NABA0s02e00106g",
                        "NABA0s02e00644g", "KAFR0D00720", "TBLA0A07050",
                        "TPHA0E03610", "ZYRO0F15818g", "KMAR0E01340",
                        "Ecym_1115", "Kwal_YGOB_matalpha2", "YCL067C",
                        "CBS432_03G00040", "Suva_3.147", "NADE0s06e00106g",
                        "CAGL0B00264g", "KNAG0C00160", "NCAS0B09140",
                        "NDAI0A00110", "TBLA0A07600", "TPHA0E04070",
                        "Kpol_2002.3", "ZYRO0F18568g", "TDEL0C07000",
                        "KLLA0C00374g", "Ecym_1002", "KLTH0F00396g",
                        "Kwal_YGOB_HMalpha2", "PKUD0OO00100",
                        "OPOL_82229", "Skud_71.2", "Suva_69.1",
                        "PAS_chr4_0878"]
    inp = args[1]
    answs = {}
    answs["a1"] = []
    answs["a2"] = []
    answs["alpha1"] = []
    answs["alpha2"] = []
    stops = {}
    stops[0] = ["TAG", "TAA", "TGA"]
    stops[1] = ["CTA", "TTA", "TCA"]
    spec = inp.strip().split(".")[0]
    print spec
    fasta = spec + ".fas"
    count = 0
    cntlines = 0
    evalfirstblast = "1e-01"
    if blastingmode == 0:
        p = subprocess.Popen("tblastn -query {0} -db {1} -out {2} -outfmt {3} -evalue {4}".format("mats.fasta", fasta, spec + ".blasted2", "6", evalfirstblast).split())
    else:
        p = subprocess.Popen("/opt/bifxapps/ncbi-blast-2.6.0+/bin/tblastn -query {0} -db {1} -out {2} -outfmt {3} -evalue {4}".format("mats.fasta", fasta, spec + ".blasted2", "6", evalfirstblast).split())
    p.communicate()
    inp = spec + ".blasted2"
    bestresul = {}
    additionallist = []
    with open(inp, "r") as inpp:
        for line in inpp:
            with open("intrs.txt", "w") as outt:
                cntlines = cntlines + 1
                #if cntlines % 10 == 0:
                #    print cntlines
                #if cntlines > 2:
                #    break
                spl = line.strip().split("\t")
                typ = spl[0].split("_")[0]
                gen = spl[1]
                querrrr = spl[0]
                #print gen
                q_st = int(spl[6])
                q_en = int(spl[7])
                s_st = int(spl[8])
                s_en = int(spl[9])
                resul = float(spl[10])
                ingore = 0
                if querrrr in bestresul.keys():
                    if resul > math.sqrt(bestresul[querrrr]):
                        ignore = 1
                        continue
                    if resul < bestresul[querrrr]:
                        bestresul[querrrr] = resul
                else:
                    bestresul[querrrr] = resul
                definitelymat = False
                if resul < 1e-25:
                    definitelymat = True
                #print line
                #print cntlines
                with open(fasta, "r") as fast:
                    writee = 0
                    curr = ""
                    found = 0
                    for line2 in fast:
                        if found == 1:
                            break
                        if gen in line2:
                            writee = 1
                        else:
                            if ">" in line2:
                                if writee == 1:
                                    found = 1
                                writee = 0
                        if writee == 1 and ">" not in line2:
                            curr = curr + line2.strip()
                #print len(curr)
                if s_st < s_en:
                    mode = 0
                else:
                    mode = 1
                if mode == 1:
                    c = s_en
                    s_en = s_st
                    s_st = c
                poss = s_st + int((s_en - s_st)/2)
                (savedbeg, savedend) = findlongest(curr, poss, mode, stops[mode])
                #search_st = s_st + 2
                #search = search_st
                #gotit = 0
                #while gotit == 0:
                #    if curr[search:(search + 3)].upper() in stops[mode]:
                #        gotit = 1
                #        earlier = search
                #    search = search - 3
                #    if search < 0:
                #        earlier = search + 3
                #        gotit = 1
                #gotit = 0
                #search = search_st
                #while gotit == 0:
                #    if curr[search:(search + 3)].upper() in stops[mode]:
                #        gotit = 1
                #        later = search
                #    search = search + 3
                #    if search >= len(curr):
                #        later = search - 3
                #        gotit = 1
                earlier = savedbeg
                later = savedend
                theseq = curr[earlier:(later + 3)]
                #print theseq
                #if typ=="a2":
                #    print earlier
                #    print later
                outt.write(">seq\n" + theseq)
            if blastingmode == 0:
                p = subprocess.Popen("blastx -query {0} -db {1} -out {2} -outfmt {3}".format("intrs.txt", "AA.fsa", "blastpython.tab", "6").split())
            else:
                p = subprocess.Popen("/opt/bifxapps/ncbi-blast-2.6.0+/bin/blastx -query {0} -db {1} -out {2} -outfmt {3}".format("intrs.txt", "AA.fsa", "blastpython.tab", "6").split())
            p.communicate()
            #if "136" in gen:
            #    p = subprocess.Popen("blastx -query {0} -db {1} -out {2}".format("intrs.txt", "AA.fsa", "ciekawewyniki.tab").split())
            #    p.communicate()
                #print theseq
            linnumb = 0
            gennam = "pusto"
            with open("blastpython.tab", "r") as blastres:
                for line3 in blastres:
                    linnumb = linnumb + 1
                    if linnumb > 1:
                        break
                    gennam = line3.strip().split("\t")[1]
                    #print gennam
            #print definitelymat
            if ([gen, earlier, (later + 3), gennam, mode, definitelymat] not in answs[typ]) and ([gen, earlier, (later + 3), gennam, mode, not definitelymat] not in answs[typ]):
                answs[typ].append([gen, earlier, (later + 3), gennam, mode, definitelymat])
            else:
                #print "else"
                if ([gen, earlier, (later + 3), gennam, mode, not definitelymat] in answs[typ]) and definitelymat:
                    additionallist.append([gen, earlier, (later + 3), gennam, mode, typ])
    #print additionallist
    for entr in additionallist:
        #print answs
        #print entr
        indee = -1
        for entr2 in answs[entr[5]]:
            indee = indee + 1
            if entr2[0] == entr[0] and entr2[1] == entr[1] and entr2[2] == entr[2] and entr2[3] == entr[3] and entr2[4] == entr[4]:
                del answs[entr[5]][indee]
                answs[entr[5]].append([entr[0], entr[1], entr[2], entr[3], entr[4], True])
                break
    #print answs
    #print "first blats finished"
    #print answs
    #print answs["a1"]
    #print answs
    ct1 = 0
    ct2 = 0
    foundgenes = {}
    foundgenes["a1"] = []
    foundgenes["a2"] = []
    foundgenes["alpha1"] = []
    foundgenes["alpha2"] = []
    for typp in typy:
        # print out numbers of candidate genes
        #if len(answs[typp]) > 0:
            #print typp + " " + str(len(answs[typp]))
        for el in answs[typp]:
            reject = 0
            another = 0
            foundother = ""
            defomat = el[5]
            #print el
            #print el[3] in known[typp]
            #print defomat
            if (el[3] not in known[typp]) and not defomat:
                if el[3] in fakesknown:
                    reject = 1
                else:
                    for typp2 in typy:
                        if el[3] in known[typp2]:
                            another = 1
                            foundother = typp2
                if reject == 0 and another == 0:
                    #foundgenes[typp].append(el)
                    #foundgenes[typp][-1].append("unknown")
                    fakesatm = []
                    with open("fakegenes.txt", "r") as rea:
                        for line9 in rea:
                            li9 = line9.strip()
                            if li9 not in fakesatm:
                                fakesatm.append(li9)
                    if el[3] not in fakesatm:
                        with open("fakegenes.txt", "a") as wri:
                            wri.write(el[3] + "\n")
                if reject == 0 and another == 1:
                    foundgenes[typp].append(el)
                    foundgenes[typp][-1].append(foundother)
                #counting candidate genes
                #ct2 = ct2 + 1
            else:
                foundgenes[typp].append(el)
                foundgenes[typp][-1].append("same")
                #counting candidate genes
                #ct1 = ct1 + 1
    #print candidate gene counts
    #print (ct1, ct2)
    #print foundgenes
    #print len("hjsgbjsgjksgn")
    #print foundgenes
    scaffolds = {}
    ordofscaf = []
    with open(fasta, "r") as fast:
        for line in fast:
            if ">" in line:
                curren = line.strip()
                ordofscaf.append(curren)
                scaffolds[line.strip()] = ""
            else:
                scaffolds[curren] = scaffolds[curren] + line.strip()
    noscaf = len(ordofscaf)
    startofscaf = [1]
    ct = 0
    for k in ordofscaf:
        ct = ct + 1
        if ct < len(ordofscaf):
            sta = startofscaf[-1] + len(scaffolds[k])
            startofscaf.append(sta)
        #print len(scaffolds[k])
    preparegenestowrite = []
    #print ordofscaf
    #print startofscaf
    for typp in foundgenes:
        for ent in foundgenes[typp]:
            toappend = []
            toappend.append(typp)
            toappend.append(ent[0])
            cnt = 0
            idx = findindex(ent[0], ordofscaf)
            startofs = startofscaf[idx]
            toappend = toappend + [ent[1] + startofs, ent[2] + startofs]
            toappend.append(spec)
            toappend.append(ent[4])
            toappend.append([ent[1], ent[2]])
            preparegenestowrite.append(toappend)
    #print preparegenestowrite
    #print "BREAK"
    #print foundgenes
    #print preparegenestowrite
    #with open("y1000_genomes/editing/cyberlindnera_jadinii.fas", "r") as f:
    #    for line in f:
    #        if "FT" in line:
    #            print list(line)
    #            print line
    #lin2 = 'FT                   /gene="'
    #print lin2
    lengg = 16000
    annotating = []
    lengg2 = 10000
    potrepeats = []
    with open("SUMMARY.txt", "a") as summary:
        a1n = 0
        a2n = 0
        alpha1n = 0
        alpha2n = 0
        scaffsdic = {}
        if len(preparegenestowrite) > 0:
            for ent in preparegenestowrite:
                typpik = ent[0]
                if typpik == "a1":
                    a1n = a1n + 1
                if typpik == "a2":
                    a2n = a2n + 1
                if typpik == "alpha1":
                    alpha1n = alpha1n + 1
                if typpik == "alpha2":
                    alpha2n = alpha2n + 1
                scaffik = ent[1]
                if scaffik in scaffsdic.keys():
                    if typpik not in scaffsdic[scaffik]:
                        scaffsdic[scaffik].append(typpik)
                else:
                    scaffsdic[scaffik] = [typpik]
        if (a1n > 0 or a2n > 0) and (alpha1n > 0 or alpha2n > 0):
            verdict = "interesting"
        else:
            verdict = "uninteresting"
        if a1n + a2n + alpha1n + alpha2n == 0:
            verdict = "nothing"
        summary.write("\t".join([spec, str(a1n), str(a2n), str(alpha1n), str(alpha2n), verdict]) + "\n")
        summary.write(str(scaffsdic) + "\n")
    if len(preparegenestowrite) > 0:
        nums2 = 0
        nums = 0
        #writegenes(preparegenestowrite)
        with open("intrs.txt", "w") as f:
            f.write("")
        #print "finding orfs"
        for gen in preparegenestowrite:
            nums = nums + 1
            ex = extractseqaround(gen, lengg, scaffolds, ordofscaf)
            #print gen
            #print ex
            ex2 = extractseqaround(gen, lengg2, scaffolds, ordofscaf)
            potrepeats.append(ex2)
            #print ex
            sseq = ex[3]
            orfs = findorfs(sseq, stops[0], stops[1])
            #print orfs
            #print orfs
            #print ex
            for ent in orfs:
                idx = findindex(ex[2], ordofscaf)
                startofs = startofscaf[idx]
                toap = [ex[2], ent[1], ent[2], startofs + ex[0] + ent[0][0], startofs + ex[0] + ent[0][1]]
                addit = 1
                for ent3 in annotating:
                    if ((toap[0] == ent3[0]) and (toap[1] == ent3[1]) and ((toap[3] == ent3[3]) or (toap[4] == ent3[4]))):
                        addit = 0
                if addit == 1:
                    nums2 = nums2 + 1
                    nammmm = ">seq" + str(nums2)
                    toap.append(nammmm)
                    toap.append([ex[0] + ent[0][0], ex[0] + ent[0][1]])
                    annotating.append(toap)
        #print "second blast"
        #print annotating
        #print len(annotating)
        smallcount = 0
        #print len(annotating)
        for gen in annotating:
            #print gen
            smallcount = smallcount + 1
            #print smallcount
            #print smallcount
            #print gen
            seqnam = gen[5]
            st = gen[3]
            en = gen[4]
            bck = gen[1]
            seq = gen[2]
            with open("intrs.txt", "a") as f:
                f.write(seqnam + "\n" + seq + "\n")
        if blastingmode == 0:
            p = subprocess.Popen("blastx -query {0} -db {1} -out {2} -outfmt {3}".format("intrs.txt", "AA.fsa", "blastpython.tab", "6").split())
        else:
            p = subprocess.Popen("/opt/bifxapps/ncbi-blast-2.6.0+/bin/blastx -query {0} -db {1} -out {2} -outfmt {3}".format("intrs.txt", "AA.fsa", "blastpython.tab", "6").split())
        p.communicate()
        #print "Done with second blast."
        #print annotating
        besteva = {}
        gennam = {}
        borderval = 1e-30
        with open("blastpython.tab", "r") as ff:
            for lin in ff:
                #print lin
                lis = lin.strip().split("\t")
                snam = lis[0]
                gnam = lis[1]
                eva = float(lis[-2])
                if eva < borderval:
                    if snam not in besteva.keys():
                        besteva[snam] = eva
                        gennam[snam] = gnam
                    else:
                        if eva < besteva[snam]:
                            besteva[snam] = eva
                            gennam[snam] = gnam
    #print besteva
    #print gennam
    #print preparegenestowrite
    #print annotating
        for ent in annotating:
            mod = ent[1]
            sca = ent[0]
            dlis = ent[6]
            snam = ent[5][1:]
            st = ent[3]
            en = ent[4]
            isitmat = 0
            if snam in besteva.keys():
                gnam = gennam[snam]
                for typ in ["a1", "a2", "alpha1", "alpha2"]:
                    if gnam in known[typ]:
                        isitmat = 1
                if isitmat == 0:
                    preparegenestowrite.append([findgennam(gnam), sca, st, en, spec, mod, dlis])
                    #print "adding"
    #print preparegenestowrite
        #print "writing genes"
        #print "MAT GENES"
        #print preparegenestowrite
        #print potrepeats
        writegenes(preparegenestowrite)
        pots = {}
        for ent in potrepeats:
            scafnam = findnam(ent[2], ordofscaf)
            st = ent[0]
            en = ent[1]
            if scafnam not in pots.keys():
                pots[scafnam] = [[st, en]]
            else:
                isitsmwhr = 0
                whichone = -1
                for rang in pots[scafnam]:
                    whichone = whichone + 1
                    if inrange(st, rang) or inrange(en, rang):
                        pots[scafnam][whichone] = [min(rang + [st, en]), max(rang + [st, en])]
                        isitsmwhr = 1
                if isitsmwhr == 0:
                    pots[scafnam].append([st, en])
        seqnum = 0
        repstowrite = []
        savedtypes = {}
        savedtypesdir = {}
        repnum = 0
        numnum = 0
        for chrom in pots.keys():
            for ent in pots[chrom]:
                #print ent
                #print chrom
                #print ent
                seqnum = seqnum + 1
                seq = scaffolds[chrom][ent[0]:ent[1]]
                numnum = numnum + 1
                with open("intrs.txt", "w") as que:
                    que.write(">seq_" + chrom[1:] + "_" + str(ent[0]) + "_" + str(ent[1]) + "\n" + seq)
                if blastingmode == 0:
                    p = subprocess.Popen("blastn -query {0} -db {1} -out {2} -outfmt {3}".format("intrs.txt", spec + ".fas", "blastpython.tab", "6").split())
                #p.communicate()
                ##p = subprocess.Popen("blastn -query {0} -db {1} -out {2}".format("intrs.txt", spec + ".fas", "blastpython_" + str(numnum) + ".tab").split())
                else:
                    p = subprocess.Popen("/opt/bifxapps/ncbi-blast-2.6.0+/bin/blastn -query {0} -db {1} -out {2} -outfmt {3}".format("intrs.txt", spec + ".fas", "blastpython.tab", "6").split())
                p.communicate()
                with open("blastpython.tab", "r") as res:
                    cntline = 0
                    for line in res:
                        #print line
                        cntline = cntline + 1
                        if cntline > 0:
                            #print "entered line"
                            spl = line.strip().split("\t")
                            perc = int(float(spl[2]))
                            #print perc
                            scafik = spl[1]
                            pocz = int(spl[8])
                            kon = int(spl[9])
                            qpocz = int(spl[6])
                            qkon = int(spl[7])
                            rev = 0
                            if kon < pocz:
                                rev = 1
                                zmiana = pocz
                                pocz = kon
                                kon = zmiana
                            mid = pocz + int((kon - pocz)/2)
                            #print "HERE IS MID"
                            #print mid
                            isitininterestingplace = 0
                            for chrom2 in pots.keys():
                                for ent2 in pots[chrom2]:
                                    #print pots
                                    #print "..."
                                    #print chrom2
                                    #print ent2
                                    #print "..."
                                    #print findnam(scafik, ordofscaf)
                                    #print findnam(chrom2, ordofscaf)
                                    #print inrange(mid, [ent2[0], ent2[1]])
                                    #print mid
                                    #print ent[0]
                                    #print ent[1]
                                    if findnam(scafik, ordofscaf) == findnam(chrom2, ordofscaf) and inrange(mid, [ent2[0], ent2[1]]):
                                        #print "WE WENT IN!!!"
                                        #print numnum
                                        #print line
                                        isitininterestingplace = 1
                            #print isitininterestingplace
                            #isitininterestingplace = 1
                            if isitininterestingplace == 1 and perc > 90:
                                #print "enter here"
                                chrtowr1 = chrom
                                scafnamtowr1 = findnam(chrtowr1, ordofscaf)
                                numtoadd = startofscaf[findindex(scafnamtowr1, ordofscaf)]
                                #print qpocz
                                #print numtoadd
                                #print ent[0]
                                #print qkon
                                #print numtoadd
                                #print ent[0]
                                sttowr1 = qpocz + numtoadd + ent[0]
                                entowr1 = qkon + numtoadd + ent[0]
                                midtowr1 = sttowr1 + int((qkon - qpocz) / 2)
                                spectowr1 = spec
                                modetowr1 = 0
                                writeintheend = 0
                                chrtowr2 = scafik
                                scafnamtowr2 = findnam(chrtowr2, ordofscaf)
                                numtoadd = startofscaf[findindex(scafnamtowr2, ordofscaf)]
                                #print pocz
                                #print numtoadd
                                #print kon
                                #print numtoadd
                                sttowr2 = pocz + numtoadd
                                entowr2 = kon + numtoadd
                                midtowr2 = sttowr2 + int((kon - pocz) / 2)
                                spectowr2 = spec
                                modetowr2 = rev
                                namtowr = "unset"
                                dontwrite1 = 0
                                dontwrite2 = 0
                                for keyy in savedtypes.keys():
                                    #print savedtypes
                                    lii = savedtypes[keyy]
                                    cnting = -1
                                    for ent4 in lii:
                                        cnting = cnting + 1
                                        if inrange(midtowr1, [ent4[0], ent4[1]]):
                                            namtowr = keyy
                                            #if namtowr == "repeat1":
                                                #print cnting
                                                #print "one"
                                            dontwrite1 = 1
                                            dirtowrite2 = (rev + savedtypesdir[keyy][cnting]) % 2
                                        if inrange(midtowr2, [ent4[0], ent4[1]]):
                                            namtowr = keyy
                                            #if namtowr == "repeat1":
                                                #print cnting
                                                #print "two"
                                            dontwrite2 = 1
                                            dirtowrite1 = (rev + savedtypesdir[keyy][cnting]) % 2
                                if namtowr == "unset":
                                    #print "FIRST TIME REPEat " + str(repnum + 1)
                                    if not (inrange(midtowr1, [sttowr2, entowr2]) or inrange(midtowr2, [sttowr1, entowr1])):
                                        repnum = repnum + 1
                                        namtowr = "repeat" + str(repnum)
                                        savedtypes[namtowr] = [[sttowr1, entowr1], [sttowr2, entowr2]]
                                        savedtypesdir[namtowr] = [0, rev]
                                else:
                                    if dontwrite1 == 0:
                                        savedtypes[namtowr].append([sttowr1, entowr1])
                                        savedtypesdir[namtowr].append(dirtowrite1)
                                    if dontwrite2 == 0:
                                        savedtypes[namtowr].append([sttowr2, entowr2])
                                        savedtypesdir[namtowr].append(dirtowrite2)
        for keyy in savedtypes.keys():
            cnt = -1
            for entt in savedtypes[keyy]:
                cnt = cnt + 1
                repstowrite.append([keyy, "irrelevant", entt[0], entt[1], spec, savedtypesdir[keyy][cnt]])
                                #if [namtowr, scafnamtowr1, sttowr1, entowr1, spectowr1, modetowr1] not in repstowrite:
                                #    repstowrite.append([namtowr, scafnamtowr1, sttowr1, entowr1, spectowr1, modetowr1])
                                #if [namtowr, scafnamtowr2, sttowr2, entowr2, spectowr2, modetowr2] not in repstowrite:
                                #    repstowrite.append([namtowr, scafnamtowr2, sttowr2, entowr2, spectowr2, modetowr2])
        #print "here1"
        #print repstowrite
        #print savedtypes
        if len(repstowrite) > 0:
            writemisc(repstowrite)




    #print preparegenestowrite
if __name__ == "__main__":
    sys.exit(main(sys.argv))
