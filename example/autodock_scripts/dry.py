#!/usr/bin/env python



""" 
Identify water molecules in close contact with a ligand in a PDB file

"""

# TODO add save best cluster option

# v.1.1 Adding support for normalizing the score by water molecules

from sys import argv
from numpy import sqrt, array
import os, getopt


#name  = argv[1].split(".")[0]

SOFTENING = .5 # softening value used in the AutoGrid maps
DEBUG = False
CUTOFF = 3.5 + SOFTENING# Angstrom. Maximum value for water interaction [ tested values : 3.15 3.51@1cbs, 3.25 1UY9; 3.5 2wi4 ]
CUTOFF = 4

DEFAULT_PENALTY = 0.3 # kcal/mol per W atom
OVERLAP_DISTANCE = 1. # angstroms to consider two W atoms overlapping for averaging
CLASH = 2.4
CLASH = 1.9
CLASH = 2.03
buffer = 5.
# TODO 
# TODO make the cutoffs agnostic of the absolute value,
# TODO but use the weight of the maps..
# TODO 
#
#
# TODO add a minimum-search protocol, for placing the water at the best spot
# TODO nearby to do a local search refinement... [ to do in AutoDock? ]

WAFFINITY_CUTOFF = -0.35 # map value to keep waters # 1efy:0.35
WAFFINITY_STRONG = -0.5 
WAFFINITY_WEAK   = -0.3

MATCH_DISTANCE = 2.0 # 1.5 + SOFTENING # to be the same as overlap_distance?

QUIET = False

lig_lig = []
lone = []
lig_rec = []

LOG=[]

cluster_number = 1
# accepted kwords for the input DLG
# INPUT:  docking.dlg [ receptor.pdbqt ]
#            if no receptor is provided, only water-water and water-ligand overlaps will be considered

# OUTPUT: cleaned docking.dlg (???)


def dist(firstline, secondline):  
    # INFO   : calculate the atomic distance between two PDB atom lines
    # INPUT  : two pdb lines
    # OUTPUT : a pdb line
    #coord1=[]
    #coord2=[]
    #temp=firstline[28:55]
    #coord1=temp.split()
    #temp=secondline[28:55]
    #coord2=temp.split()
    # floating the strings
    #for index in range(len(coord1)):
    #    coord1[index]=float(coord1[index])
    #    coord2[index]=float(coord2[index])
    coord1 = coord(firstline)
    coord2 = coord(secondline)
    measure=sqrt((coord1[0]-coord2[0])**2+(coord1[1]-coord2[1])**2+(coord1[2]-coord2[2])**2)
    return measure

def pc(v, tot):
    return ( (float(v)/float(tot)) * 100)

def get_lines(filename):
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()
    return lines

"""
def get_coords(filename = None, list = None, include_hydrogens = True): # TODO DISABLED
    " "" 
    designed to extract receptor coordinates from a PDB(QT)
    but it can be used with lists too.
    return [ text_pdb_atom_entries, numeric_array_coords, atype]
    " ""
    coord = []
    atoms = []
    atype = []
    if not filename and not list:
        return False
    if filename:
        source = get_lines(filename)
    if list:
        source = list
    for l in source:
        if l.startswith("ATOM") or l.startswith("HETATM"):
            at = l.rsplit(None, 1)[1]
            if not at == "HD" or include_hydrogens: # by default, HD are included
                coord.append([float(l[30:38]),float(l[38:46]),float(l[46:54])])
                atoms.append(l)
                atype.append(at) 
    return { 'text' : atoms, 'coord' : array( coord, 'f'), 'atype': atype }
"""

def coord(l):
    return [float(l[30:38]),float(l[38:46]),float(l[46:54])]


def parse_dlg (filename, raw_output = False):
    # INPUT : DLG file (or pdbqt?)
    # OUTPUT: pdb-compliant lines
    # OPTIONS: raw_output (True/False) get raw "DOCKED: " lines or the clustered lines

    file = open(filename, 'r')
    dlg = file.readlines()
    file.close()

    accepted_kw = [ "MODEL",
        "USER",
        "ATOM",
        "HETATM",
        "TER",
        "ENDMDL"
        ]
    raw_keywords = ["DOCKED:", "INPUT-PDBQ:"]

    output = []
    for line in dlg:
        line = line.strip()
        if line:
            try:
                kw = line.split()[0]
            except:
                pass
            if not raw_output:
                if kw in accepted_kw:
                    output.append(line)
            else:
                if kw in raw_keywords:
                    output.append( line.split(": ", 1)[1])
    return output

def parse_pdbqt (filename):
    # INPUT : DLG file (or pdbqt?)
    # OUTPUT: pdb-compliant lines
    # 
    # TODO add the flexible residues support

    fp = open(filename, 'r')
    dlg = fp.readlines()
    fp.close()

    accepted_kw = [ "MODEL",
        "USER",
        "ATOM",
        "HETATM",
        "TER",
        "ENDMDL",
        "REMARK",
        'ROOT',
        'ENDROOT',
        'BRANCH',
        'ENDBRANCH'
        'TORSDOF',
        ]

    output = []
    for line in dlg:
        try:
            kw = line.split()[0]
        except:
            pass #print line
        if kw in accepted_kw:
            output.append(line)
    return output

def get_lelc(filename, cluster_number = None):
    # get clustering information from a DLG
    # histogram graph.
    biggest = 0
    fp = open(filename, 'r')
    dlg = fp.readlines()
    fp.close()
    IN = False
    current_cluster = 0
    for line in dlg:
        if "_____|___________|_____|___________|_____|______________________________________" in line:
            IN = False

        if IN:
            tmp = line.split("|")
            cluster_size = int(tmp[4])
            current_cluster += 1
            if cluster_number and current_cluster == cluster_number:
                return [ int(tmp(0)), cluster_size, tmp[5] ] # TODO TO TEST!
            elif cluster_size > biggest:
                biggest = cluster_size
                cluster = int(tmp[0])
                run = int(tmp[2])
        else:
            if "RANKING" in line:
                tmp = line.split()
                if int(tmp[0]) == cluster and int(tmp[1]) == 1: # pick the subrank 1 of the most populated cluster
                    rmsd = tmp[5]
                    return [cluster, biggest, rmsd, run]

        if "_____|___________|_____|___________|_____|____:____|____:____|____:____|____:___" in line:
            IN = True




def parse_models(structure):
    models = []
    buffer = []

    for line in structure:
        if not "ENDMDL" in line:
            buffer.append(line)
        else:
            buffer.append(line)
            models.append(buffer)
            buffer = []

    if models:
        if not QUIET: print(" [ %d docking pose(s) found ] " % len(models), end=' ')

        return models
    else:
        if not QUIET: print("\t\t\t\t [ NO DOCKING POSES HAVE BEEN FOUND!!!! ] ")
        return False




def minmax(models):
    # INPUT : pdb-like line
    # OUTPUT: coordinates (x,y,z), (X,Y,Z)
    """ Return max/min coordinates of
    the box including the given atoms list"""

    x_min, y_min, z_min =  999999999,  9999999999,  9999999999
    x_max, y_max, z_max = -999999999, -9999999999, -9999999999
    for pose in models:
        for line in pose:
            if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                atom = coord(line)
                # min
                if atom[0] < x_min:
                    x_min = atom[0]
                if atom[1] < y_min:
                    y_min = atom[1]
                if atom[2] < z_min:
                    z_min = atom[2]
                # max
                if atom[0] > x_max:
                    x_max = atom[0]
                if atom[1] > y_max:
                    y_max = atom[1]
                if atom[2] > z_max:
                    z_max = atom[2]
    return [x_min, y_min, z_min], [x_max, y_max, z_max]


def load_receptor(filename, filter = True):
    file = open(filename, 'r')
    stack = file.readlines()
    file.close()
    if filter:
        output = []
        for line in stack:
            if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                output.append(line)
        return output
    else:
        return stack

def in_the_box(atom_list, MIN, MAX):
    # INPUT : pdb-like atom, MIN = [x,y,z], MAX = [x,y,z]
    # OUTPUT: True/False

    good = []
    for atom in atom_list:
        pos = coord(atom)
        if pos[0] < MAX[0]+buffer and pos[0] > MIN[0]-buffer:
            if pos[1] < MAX[1]+buffer and pos[1] > MIN[1]-buffer:
                if pos[2] < MAX[2]+buffer and pos[2] > MIN[2]-buffer:
                    good.append(atom)
    if good:
        return good
    else:
        return False

    
def atom_type(line):
    # make it to work with both PDB and PDBQT
    atype = ""
    line = line.strip()
    value = line.split()[2]
    value = line.split()[-1]
    for c in value:
        if c.isalpha():
            atype += c
    return value


def pdb_atom_type(line):
    atype = ""
    tmp = line[12:16].strip()
    for c in tmp:
        if c.isalpha():
            atype+= c
    return atype

def pdbqt_atom_type(line):
    return line.split()[-1].strip()


### TODO create a function to extract input ligand and corresponding atom type
def input_ligand_pdbqt(dlg):
    ligand = []
    atypes = []
    for line in dlg:
        if line[0:18] == "INPUT-LIGAND-PDBQT":
            line = line.split(None, 1)[1]
            ligand.append(line)
            splitting = line.split()
            if splitting[0] == "ATOM" or splitting[0] == "HETATM":
                atypes.append(splitting[-1])
    return ligand, atypes



def selfoverlap(model, cutoff = CLASH):
    global LOG
    # remove ligand-overlapping w-atoms
    # INPUT : pdb model
    # OUTPUT: cleaned pdb model
    count = 0
    h_cutoff = cutoff - 0.7 # 0.7 first value; 0.6
    for pose in model:
        waters = []
        atoms  = []
        lig_lig.append("MODEL\n")
        idx = model.index(pose)
        for entry in pose:
            if entry[0:4] == "ATOM" or entry[0:6] == "HETATM":
                if pdb_atom_type(entry) == "W":
                    waters.append(entry)
                else:
                    atoms.append(entry)
        if DEBUG:
            LOG.append("MODEL\nREMARK LIGAND_OVERLAP_WATERS [ pose %d ]\n###################\n" % (idx))
        for w in waters:
            try:
                for a in atoms:
                    if not pdb_atom_type(a) == "H" and not pdbqt_atom_type(a)[0] == "H":
                        distance = dist(w,a)
                        if distance < cutoff:
                            count += 1
                            lig_lig.append(w)
                            pose.remove(w)
                            if DEBUG:
                                LOG.append("REMARK water/ligand_atom")
                                LOG.append(w)
                                LOG.append(a)
                                LOG.append("REMARK dist: "+str(distance)+" [max : "+str(cutoff)+" ]"+atom_type(a)+"\n\n")
                            raise
            except:
                pass
        if DEBUG: LOG.append("ENDMDL")
        lig_lig.append("ENDMDL\n")
        model[idx] = pose
    if not QUIET: print(" [ %d water(s) removed ]" % count)
    return model


def get_waters(pdbqt_poses, type = "W"):
    waters = []
    for p in pdbqt_poses:
        for a in p:
            if a[0:4] == "ATOM" or a[0:6] == "HETATM":
                if pdbqt_atom_type(a) == type:
                    waters.append(a)
    if not QUIET: print(" [ filtered pose contains %d waters ]" % len(waters))
    return waters


def targetoverlap(model, target, cutoff = CLASH):
    global LOG
    # INPUT : pdb lig (multi)model and pdb target model
    # OUTPUT: cleaned pdb lig model
    h_cutoff = cutoff - 0.4 # TODO decrease to 0.4? keep 0.7?
    h_cutoff = 0 # TODO CHECK THIS EXCLUDING HD!!!
    count = 0
    for pose in model:
        lig_rec.append("MODEL\n")
        waters = []
        idx = model.index(pose)
        for entry in pose:
            if entry[0:4] == "ATOM" or entry[0:6] == "HETATM":
                if pdb_atom_type(entry) == "W":
                    waters.append(entry)
        if DEBUG:
            LOG.append("MODEL\nREMARK TARGET_OVERLAP_WATERS [ pose %d ]\n###################\n" % (idx))
        for w in waters:
            try:
                for a in target:
                    if not pdbqt_atom_type(a) == "HD":
                        max = cutoff
        
                        if dist(w, a) <= max:
                            
                            lig_rec.append(w)
                            lig_rec.append(a)
                            pose.remove(w)
                            count += 1
                            if DEBUG:
                                print("water kicked out:", w, a)
                                LOG.append("REMARK water/receptor_atom")
                                LOG.append(w)
                                LOG.append(a)
                                LOG.append("REMARK dist: "+str(distance)+" [max : "+str(max)+" ]"+atom_type(a)+"\n\n")
                            raise
            except:
                pass
        if DEBUG: LOG.append("ENDMDL")
        lig_rec.append("ENDMDL ")
        model[idx] = pose
    if not QUIET: print(" [ %d water(s) removed ]" % count)
    return model


def lonewaters(model, target, cutoff = CUTOFF): # tested values : 3.15 3.51@1cbs, 3.25 1UY9; 3.5 2wi4
    hcutoff = cutoff - 0.7
    count = 0
    for pose in model:
        waters = []
        mates = []
        idx = model.index(pose)
        for entry in pose:
            if entry[0:4] == "ATOM" or entry[0:6] == "HETATM":
                if atom_type(entry) == "W":
                    waters.append(entry)
        for entry in target:
            if entry[0:4] == "ATOM" or entry[0:6] == "HETATM":
                atype = entry.split()[-1]
                if atype == "NA" or atype == "OA" or atype == "HD":
                    mates.append(entry)
        
        lone.append("MODEL ")
        for w in waters:
            safe = False
            for a in mates:
                atype = a.split()[-1]
                if atype == "HD":
                    max = hcutoff
                else:
                    max = cutoff
                if dist(w, a) < max:
                    if DEBUG:
                        LOG.append("LONE_WATERS\n###[ pose"+str(idx)+"]######\nwater: "+w+"r_atm: "+a+"dist:"+str(dist(w,a)))
                    safe = True
                    break
            if not safe:
                count += 1
                lone.append(w)
                pose.remove(w)
        model[idx] = pose
        
        lone.append("ENDMDL ")
    
    if not QUIET: print("\t\t\t [ %d water(s) removed ]" % count)
    return model



def average(list):
    # INPUT : pdb atom/het lines
    # OUTPUT: average coordinates

    total = [ 0., 0., 0.]
    for atom in list:
        a_coords = coord(atom)
        total[0] += a_coords[0]
        total[1] += a_coords[1]
        total[2] += a_coords[2]

    avg = [ total[0]/len(list), total[1]/len(list), total[2]/len(list) ]
    return avg


def pdbatom(coords, kw = "ATOM  ", atype = "X", pcharge = 0.):
    # INPUT : three coordinates
    # OUTPUT: a pdb atom
    index = 99
    residue = "WAT"
    chain = 9
    line = "%s%5d %2s   %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00     %1.3f %1s" % (kw, index, atype, residue, chain, index, coords[0], coords[1], coords[2], pcharge, atype)
    return line


def cluster_waters(model, cutoff = OVERLAP_DISTANCE):
    # TODO decide if report the pop size as:
    #    W3+
    #   b-factor
    #   partial charge   <-

    count = 0
    for pose in model:
        waters = []
        candidates = []
        idx = model.index(pose)
        for entry in pose:
            if entry[0:4] == "ATOM" or entry[0:6] == "HETATM":
                if pdbqt_atom_type(entry) == "W":
                    waters.append(entry)
                    pose.remove(entry)
        print(" [cluster: found %d waters]" % len(waters))
    
        for w in waters:
            pop = [w]
            for mate in waters:
                if not mate == w:
                    if dist(mate, w) <= cutoff:
                        pop.append(mate)
                        waters.remove(mate)
                        waters.remove(w)
            if len(pop)> 0:
                candidates.append(pop)
    
        for pop in candidates:
            count = len(pop)
            avg = average(pop)
            avg_atom = pdbatom(avg, atype = "W", pcharge = count)
            pose.append(avg_atom)
        model[idx] = pose

    if count > 1:
        if not QUIET: print("\t\t\t [ %d clustered ]" % count)
    else:
        if not QUIET: print("\t\t\t [ none ]") 
    return model
            

def save_pdb(list, filename):
    file = open(filename, 'w')
    for i in list:
        print(i[:-1], file=file)
    file.close

    
def save_output_pdbqt(dictionary, output_filename, count = 0, nowaters = False, index = False, printrms = False):
    # NOTE output_filename is a pointer
    # count = 0 will process all the structures
    # any other value will keep going...
    output = []

    if index:
        index -=1
        model = dictionary[index]
        for item in model:
            if item[0:4] == "ATOM" or item[0:6] == "HETATM":
                atype = atom_type(item)
                if atype == "W" and nowaters:
                    pass
                else:
                    item = item[0:77]+atype
                    output.append(item)
            else:
                item = item[:-1]
                output.append(item)

    else:
        for model in dictionary:
            for item in model:
                if printrms and "RMSD from reference structure" in item:
                    print(item.split()[6])
                if item[0:4] == "ATOM" or item[0:6] == "HETATM":
                    atype = atom_type(item)
                    if atype == "W" and nowaters:
                        pass
                    else:
                        item = item[0:77]+atype
                        output.append(item)
                else:
                    item = item[:-1]
                    output.append(item)
            if not count == 0:
                break

    for line in output:
        print(line, file=output_filename)




def experimental_match( calc_list, exp_list, cutoff = MATCH_DISTANCE):
    real = len(exp_list)
    predicted = len(calc_list)
    if not QUIET: print("real/predicted", real,predicted)
    
    found = 0

    min = 9999999
    stack = []
    invented = 0
    if real > predicted:
        missed = real-predicted
    if real < predicted:
        invented = predicted-real
    for e in exp_list:
        for w in calc_list:
            if (not e in stack) and (not w in stack):
                measure = dist(w,e)
                if measure <= cutoff:
                    found += 1
                    stack.append(w)
                    stack.append(e)
                if measure < min:
                    min = measure
    missed = real-found
    success = (float(found)/float(real))*100
    return (found, success, missed, invented)



def confirm_predictions(calc_list, exp_list, cutoff = MATCH_DISTANCE):
    real = len(exp_list)
    predicted = len(calc_list)
    if DEBUG: print("confirm_predictions>real/predicted", real,predicted)
    
    found = 0
    distance_log = []
    for w in calc_list:
        for e in exp_list:
            d = dist(w,e)
            if d < cutoff:
                found += 1
                distance_log.append([w, e, d])
                break
    return [found, distance_log]

                


def getmap(data):

    spacing = float(data[3].split()[1])

    pts = data[4].split()[1:]
    for i in range(len(pts)):
        pts[i] = int(pts[i])+1

    center = data[5].split()[1:]
    for i in range(len(pts)):
        center[i] = float(center[i])

    data = data[6:]
    for i in range(len(data)):
        data[i] = float(data[i])

    step = []
    step.append(pts[0]/2 * spacing)
    step.append(pts[1]/2 * spacing)
    step.append(pts[2]/2 * spacing)

    min, max = [], []

    min.append(center[0]-step[0])
    min.append(center[1]-step[1])
    min.append(center[2]-step[2])

    max.append(center[0]+step[0])
    max.append(center[1]+step[1])
    max.append(center[2]+step[2])

    data = array(data, 'f').reshape(pts[2], pts[1], pts[0])

    data = data.reshape(pts[2], pts[1], pts[0])
    return { "values" : data, "spacing" : spacing, 'pts': pts, 'center' : center, 'min' : min, 'max' : max}
    


def getpoints(mapdata, coords=None, distance=None):
    data = mapdata['values']
    min = mapdata['min']
    max = mapdata['max']
    spacing = mapdata['spacing']
    pts = mapdata['pts']
    harvesting = []
    
    if not distance:
        distace = 1.0 # default 1A distance scanning

    if not coords and not distance:
        coords = [3., 4., 5.]
        distance = 1.0

    pt_scan = int( round(float(distance)/(spacing)) )
    if pt_scan == 0: pt_scan = 1 # at least one point around
    if DEBUG: print("Grid range: %2.2f [ +/- %d points : %2.2f ]" % (distance,pt_scan, pt_scan*spacing), end=' ') 
    best = 999

    if (coords[0] < max[0]) and (coords[0] > min[0]):
        if (coords[1] < max[1]) and (coords[1] > min[1]):
            if (coords[2] < max[2]) and (coords[2] > min[2]):
                z_pt = int(round((coords[2] - min[2])/spacing))
                y_pt = int(round((coords[1] - min[1])/spacing))
                x_pt = int(round((coords[0] - min[0])/spacing))

                for x_ofs in range(-pt_scan, pt_scan+1): 
                    x = x_pt + x_ofs
                    if x < 0 :
                        harvesting.append(0)
                        break
                    for y_ofs in range(-pt_scan, pt_scan+1): 
                        y = y_pt + y_ofs
                        if y < 0 :
                            harvesting.append(0)
                            break
                        for z_ofs in range(-pt_scan, pt_scan+1): 
                            z = z_pt + z_ofs
                            if z < 0 :
                                harvesting.append(0)
                                break
                            harvesting.append( data[z,y,x] )
                            if data[z,y,x] < best:
                                best = data[z,y,x]
            if DEBUG: print(" %d pts [ best: %2.2f ]" % (len(harvesting), best), end=' ')
            return harvesting, best
    return False




def scoreCorrection(atoms, atom_counts, penalty):
    """ ugly function, just temporary for helping Yang.
        To be rewritten in a more clean way...

        it works only with PDBQT+ structures!
    """
    #penalty = 0.3 # same as the entropy?
    heavy = atom_counts[0]
    waters = atom_counts[1]
    penalty = waters * penalty
    if not QUIET: print(" energy correction (penalty = %2.3f per %d W atom) " % (penalty,waters))
    if waters == 0: return atoms

    hpattern = "USER    AD_histogram> "
    epattern = "USER    AD_LE> "
    cpattern = "USER    AD_LC> "
    LONG=False
    for i in range(len(atoms)):
        l = atoms[i]
        #print l
        # histogram
        if hpattern in l:
            histo = l.split(hpattern)[1]
            histo = histo.split(",")
            for j in range(len(histo)):
                e,bulk = histo[j].split(":",1)
                histo[j] = "%2.3f:%s" % (float(e)-penalty, bulk)
            correct = hpattern
            for h in histo:
                correct += h
                correct += ","
            correct = correct[:-1]+"\n"
            atoms[i] = correct
        # le
        elif epattern in l:
            info = l.split(epattern)[1]
            # USER    AD_LE> -9.310,	-0.776,	1,	100.00

            info = info.split(",",2)
            e_org = float(info[0])
            leff_org = float(info[1])
            e = e_org + penalty
            leff = e/heavy
            correct = epattern
            correct += " %2.3f," % e
            correct += " %2.3f," % leff
            correct += info[2]
            atoms[i] = correct
        # lc
        elif cpattern in l:
            LONG=True

            info = l.split(cpattern)[1]
            # USER    AD_LC> -9.310,	-0.776,	1,	100.00
            info = info.split(",",2)
            ec_org = float(info[0])
            leffc_org = float(info[1])
            ec = ec_org + penalty
            leffc = ec/heavy
            correct = epattern
            correct += " %2.3f," % ec
            correct += " %2.3f," % leffc
            correct += info[2]
            atoms[i] = correct


    if not QUIET:
        print("\t(LE) energy            : %2.3f \t[org: %2.3f ]" % (e, e_org))
        print("\t     ligand efficiency : %2.3f \t[org: %2.3f ]" % (leff, leff_org))

        if LONG:
            print("\t(LC) energy            : %2.3f \t[org: %2.3f ]" % (ec, ec_org))
            print("\t     ligand efficiency : %2.3f \t[org: %2.3f ]" % (leffc, leffc_org))
    else:
        print("%s,%2.3f,%2.3f" % (input_file,e,leff), end=' ')
        if LONG:
            print(",%2.3f,%2.3f" % (ec,leffc))
        else:
            print("\n")
            
    return atoms



def countWatoms(madels,excludeH = True):
    w = 0
    a = 0
    for l in models[0]:
        if l.startswith("ATOM") or l.startswith("HETATM"):
            if l.split()[-1] == "W": w+=1
            else: 
                if not "HD" in l.split()[-1]:  a+=1
    return a,w
            




def usage():
    myname =  os.path.basename(argv[0])

    print("""\n
                                             ,---,  
                                          ,`--.' |  
                ,---,                     |   :  :  
              .'  .' `\                   '   '  ;  
            ,---.'     \   __  ,-.        |   |  |  
            |   |  .`\  |,' ,'/ /|        '   :  ;  
            :   : |  '  |'  | |' |   .--, |   |  '  
            |   ' '  ;  :|  |   ,' /_ ./| '   :  |  
            '   | ;  .  |'  :  /, ' , ' : ;   |  ;  
            |   | :  |  '|  | '/___/ \: | `---'. |  
            '   : | /  ; ;  : | .  \  ' |  `--..`;  
            |   | '` ,/  |  , ;  \  ;   : .--,_     
            ;   :  .'     ---'    \  \  ; |    |`.  
            |   ,.'                :  \  \`-- -`, ; 
            '---'                   \  ' ;  '---`"  
                                     `--`           

    """)

    print("\tUSAGE :\t%s [options] -i input.[dlg|pdbqt] \n" % myname)
    print("\t     -i\t\tinput file [ required ]. Accepted file formats are '.dlg', '.pdbqt' (including pdbqt+)\n")

    print("\tOptions")
    print("\t-------")
    print("\t    -F\t\tfull mode (write also the DRY,STRONG,WEAK files to cope with ADT limitations)")
    print("\t    -E\t\textract only the lowest energy structure, strip all waters (if present) and save the file")
    print("\t    -c\t\textract the lowest energy-largest cluster result ( DLG input only )")
    print("\t    -C\t\textract the lowest energy-largest cluster result, strip all waters (if present) and save the file (DLG input only)")
    print("\t    -p [integer]\t\textract the selected cluster number")

    print("\t    -n [float]\t\tnormalize the energy by using this value per W atoms [ default : %2.3f ]\n" % (DEFAULT_PENALTY))
    print("\tStructural input")
    print("\t----------------")
    print("\t    -r rec.pdbqt\tuse the receptor PDBQT file")
    print("\t    -w waters.pdb\tuse the experimental waters in the file to calculate distances and success rate")

    print("\n\tW-atoms scoring & cleaning")
    print("\t--------------------------")
    print("\t    Maps")
    print("\t    -m protein.X.map\tuse map file to decide wich waters are kept")
    print("\t    -x float\tmap affinity cutoff [default %2.2f ]" % WAFFINITY_CUTOFF)
#    print "\t    Distance"
#    print "\t    -M [float]\tmerge W-atoms that are closer than the distance [ default : %1.2f Angstrom ]" % OVERLAP_DISTANCE


    print("\t    -q\t\tquiet mode. Only dir name and the RMSD are printed out.")
    print("\t    -d\t\tsave the \"DEBUG_*\" files\n")

    print("\t\tIf the receptor.pdbqt is omitted, only self-overlapping ligand waters will be removed")
    print("\t\tBy default, the following filters are performed:")
    print("\t\t\t - ligand-atoms overlapping waters ")
    print("\t\t\t - hanging waters (not interacting with any protein atoms)")
    print("\t\t\t - off-waters (waters overlapping with receptor atoms [if provided])")
    print("\t\t\t - clustering of close waters (DISABLED!)")

    print("\n\tReference  ")
    print("\t----------")
    print("\t     Please cite:  Stefano Forli and Arthur J. Olson, J. Med. Chem., 2012, 55 (2), pp 623-638")
    print("\t                   DOI: 10.1021/jm2005145")
    
    print("\n\n")



# OPTIONS
try:
    options, structures = getopt.getopt(argv[1:], 'EhdcCqr:w:m:M:x:p:n:Fi:')
except getopt.GetoptError as err:
 usage()
 print(str(err), "\n")
 exit(2)

opts = {} # create the option dictonary with {option: argument} format

for o, a in options: # populate the options
 opts[o] = a

if '-i' in opts:
    input_file = opts['-i']
else:
    print("ERROR! Input file ('-i') is required!")
    usage()
    exit(0)


if '-h' in opts: # or '--help' in opts:
 usage()
 exit(0)

if '-q' in opts:
    QUIET = True

if '-F' in opts: FULL = True
else: FULL = False

if '-d' in opts:
  print("\n\n\t[ Debugging mode is *ON* ]\n\n")
  DEBUG = True

if not QUIET: 
    print("""
                  ____                      
                 /\  _`\                    
                 \ \ \/\ \  _ __  __  __    
                  \ \ \ \ \/\`'__\\\\ \/\ \   
                   \ \ \_\ \ \ \/\ \ \_\ \  
                    \ \____/\ \_\ \/`____ \ 
                     \/___/  \/_/  `/___/> \\
                                      /\___/
                                      \/__/ 

    """)
    print("========================== INPUT DATA =========================")



PDBQTplus = False
try:
    #if True:
    #name, ext = os.path.splitext(structures[0])
    name, ext = os.path.splitext(input_file)
    if not ext == ".dlg" and ('-c' in opts or '-C' in opts or '-E' in opts ):
        print("ERROR: clustering-related options can be used only with DLG files")
        exit(1)

    if ext == ".dlg":
        if not QUIET: print(" [ using DLG file %s ]" % input_file) #(structures[0])
        docking = parse_dlg(input_file, raw_output = True) # TODO it must become: docking, header, footer
        if not QUIET: print("    parsing docked MODELS...", end=' ')
        models = parse_models(docking)
        if not len(models):
            print("no docking results in the file")
            exit(1)
    else:
        # TODO put a function for loading simple pdbqt
        if not QUIET: print(" importing ATOMS from ", input_file) # structures[0]
        models = [parse_pdbqt(input_file) ] # TODO it must become: docking, header, footer
        # CHECK HERE IF THERE'S A PDBQT+
        tmp = get_lines(input_file)
        for l in tmp:
            if "ADVS_result" in l:
                PDBQTplus=True
                atom_counts = countWatoms(models)
                if not QUIET: print(" [found PDBQT+ : %d atoms ; %d waters ]" % (atom_counts))
                break
        if not len(models):
            print("no docking results in the file")
            exit(1)
        pass
    #print "\n",
except:
    #else:
    if not QUIET:
        print("ERROR: input file is missing")
        usage()
        #exit(1)
    #else:
    #    print "[missing %s]" % structures[0]
    exit(1)


if '-n' in opts:
    try:
        norm = float(opts['-n'])
        print(" normalize energy %2.3f Kcal/mol x W atom" % norm)

    except:
        print("Error in the normalization value!")
        usage()
        exit(1)

else: norm = DEFAULT_PENALTY 


if '-M' in opts:
    # cutoff for overlapping distance between W atoms
    CLUSTER = True
    try:
        OVERLAP_DISTANCE = float(opts['-M'])
    except:
        pass
else:
    CLUSTER = False
    

if '-m' in opts:
    try:
        if not QUIET: print("\n [ using map file %s ]" % (opts['-m']))
        grid = get_lines(opts['-m'])
        if '-x' in opts:
            try:
                WAFFINITY_CUTOFF = -float(opts['-x'])
            except:
                print("Error in maps cutoff value: %s" % opts['-x'])
                exit(1)
            if not QUIET: print("Affinity cutoff = ", WAFFINITY_CUTOFF)
    except:
        print("ERROR! Impossible to open map file : %s" % opts['-m'])
        exit(1)
else:
    grigrid = False



if '-w' in opts:
    if not QUIET: print("\n [ getting waters from %s ]" % opts['-w'])
    waters = []
    tmp = get_lines(opts['-w'])
    for a in tmp:
        if atom_type(a) == "O":
            waters.append(a)
else:
    waters = False


if not QUIET:
    print("===============================================================\n")



if '-E' in opts:
    if not QUIET: print("\n[ saving the lowest energy result only ]")
    output_filename = name+"_RINSE_LE.pdbqt"
    out = open(output_filename, 'w')
    pdb = os.getcwd()
    pdb = os.path.basename(pdb)
    print(pdb+","+name+",", end=' ')
    save_output_pdbqt(models, out, count = 1, nowaters = True, printrms = True)
    exit(0)



if '-C' in opts:
    out = open(name+"_RINSE_LELC.pdbqt", 'w')
    lelc = get_lelc(structures[0])
    if not len(lelc):
        print("no docking results in the file")
        exit(1)
    if not QUIET:
        print("\n[ saving the rinsed LELC result only ]")
    else:
        pdb = os.getcwd()
        pdb = os.path.basename(pdb)
        print(pdb+","+name+","+lelc[2])
    save_output_pdbqt(models, out, count = 0, nowaters = True, index=lelc[0] )
    exit(0)



# Check if we're in receptor mode
if '-r' in opts:
    try:
        rec = load_receptor(opts['-r'])
        if not QUIET: print("\n receptor structure loaded\t", end=' ')
        if not QUIET: print("\t\t [ %d atoms ]" % len(rec))
        MIN, MAX = minmax (models)
        # include only receptor atoms in the max/min ligand coordinates to speed up the computation
        rec_box = in_the_box(rec, MIN, MAX)
        if DEBUG: save_pdb(rec_box, "DEBUG_box.pdb")
        if not QUIET: print(" receptor 5A shell extracted ", end=' ')
        if not rec_box:
            print("ERROR: no receptor atoms are close to the ligand")
            exit(1)
        else:
            if not QUIET: print("\t\t\t [ %d atoms in 5 A shell ] " % len(rec_box))
        REC = True
    except:
        print("ERROR: something bad happened when reading the receptor file", rec)
        exit(1)
else:
    if not QUIET: print(" [ no receptor ]")
    REC = False


if not models:
    print("ERROR: no docking poses have been found")
    exit(1)


suffix = "" #name



#print "\n-PROCESSING THE DATA --------------------------"

### Pose extraction
if '-c' in opts:
    lelc = get_lelc(input_file)
    #lelc = get_lelc(structures[0])
    lelc_idx = lelc[3]
    if not QUIET:
        print("\n LELC docked pose\t\t\t\t [ biggest cluster: ", lelc[1], " AutoDock RMSD :", lelc[2], " Run #", lelc[3], "] \n")

    models = [ models[ lelc_idx-1 ] ]
    suffix = "_LELC"
elif '-e' in opts:
    models = [ models[0]]
    print("WARNING! NEVER TESTED!")
    suffix = "_LE"

elif '-p' in opts:
    try:
        cluster_number = int(opts['-p'])
        if not QUIET: print(" [extracting cluster # %d ]" % (cluster_number))
        models = [ models[cluster_number-1]]
        suffix = ("_pose_%d" % int(opts['-p']))

    except:
        print("ERROR! please specify a number")
        exit(1)
### Pose extraction

if not QUIET: print(" removing ligand/ligand overlapping waters\t", end=' ')
clean = selfoverlap(models)
suffix += "_DRY"

if REC:
    if not QUIET: print(" removing ligand/receptor overlapping waters\t", end=' ') 
    clean = targetoverlap(clean, rec_box)

if CLUSTER:
    if not QUIET: print(" clustering remaining waters", end=' ')
    clean = cluster_waters(clean)
    #print "CLEAN1",len(clean)

def write_this(filename, list):
    fp = open(filename, 'w')
    for a in list:
        print(a.strip(), file=fp)
    fp.close()

if grid :
    if not QUIET: print("\n scanning grid map for conserved waters...\t", end=' ') 
    grid = getmap(grid)

    watoms = get_waters(clean)
    #print "len", len(watoms)
    #print "CLEAN2", len(clean)
    conserved_waters_s = {}
    conserved_waters_w = {}

    if not QUIET: print("\n water grid score results [ map: %s ] " % (opts['-m']))
    for w in watoms:
        harvest, affinity = getpoints(grid, coord(w), distance = 1.0) # angstroms
        #WAFFINITY_STRONG = -0.45 
        #WAFFINITY_WEAK   = -0.3

        # WAFFINITY_CUTOFF = -0.35 # map value to keep waters # 1efy:0.35
        # WAFFINITY_STRONG = -0.5 
        # WAFFINITY_WEAK   = -0.3

        # strong water
        if affinity < WAFFINITY_STRONG and not (w.strip() in conserved_waters_s):
                if not QUIET: 
                    print("\t [ Water STRONG ( %2.2f ) +++ ]" % affinity)
                conserved_waters_s[w] = affinity

        # weak water
        #elif (affinity < WAFFINITY_WEAK) and not (affinity > WAFFINITY_CUTOFF) and not (w.strip() in conserved_waters_w):
        elif (affinity < WAFFINITY_WEAK) and not (w.strip() in conserved_waters_w):
                if not QUIET: print("\t [ Water  WEAK  ( %2.2f )  +  ]" % affinity)
                conserved_waters_w[w] = affinity
        # displaced water
        else:
            if not QUIET: print("\t [ Water DISPLC ( %2.2f )  D  ]" % affinity)

    if DEBUG:
        confirmed_water = open("test_maps_confirmed.pdb", 'w')
        if conserved_waters_s:
            print("MODEL\nREMARK Strong waters", file=confirmed_water)
            for w in conserved_waters_s:
                print(w.strip(), file=confirmed_water)
                print(("REMARK map affinity : %2.3f" % conserved_waters_s[w]), file=confirmed_water)
            print("ENDMDL", file=confirmed_water)
        if conserved_waters_w:
            print("MODEL\nREMARK Weak waters", file=confirmed_water)
            for w in conserved_waters_w:
                print(w.strip(), file=confirmed_water)
                print(("REMARK map affinity : %2.3f" % conserved_waters_w[w]), file=confirmed_water)
            print("ENDMDL", file=confirmed_water)
        confirmed_water.close()

    suffix += "_SCORED"

    atoms = []
    dry_version = []
    weak_w = []
    strong_w = []
    #print "NAME", name
    #print "NAME NAME NAME = ", name+suffix+'.pdbqt'
    for a in models[0]:
        if pdbqt_atom_type(a) == "W":
            if a in list(conserved_waters_s.keys()):
               atoms.append(("REMARK  STRONG water ( score: %2.2f )" % conserved_waters_s[a]))
               strong_w.append(a)
               a = a[0:55]+" 1.00 30.00           "+a[77:]
               atoms.append(a)

            elif a in list(conserved_waters_w.keys()):
               atoms.append(("REMARK  WEAK   water ( score: %2.2f )" % conserved_waters_w[a]))
               weak_w.append(a)
               a = a[0:55]+" 1.00 80.00           "+a[77:]
               atoms.append(a)
        else:
            dry_version.append(a)
            if a.startswith("ATOM") or a.startswith("HETATM"):
                a = a[0:55]+" 1.00 10.00           "+a[77:]
            atoms.append(a)
            
    output = open(name+suffix+".pdbqt", 'w')
    if not PDBQTplus:
        for a in atoms:
            output.write(a+'\n')
        output.close()
    else:
        atoms = scoreCorrection(atoms, atom_counts, norm)
        for a in atoms:
            output.write(a)
        output.close()

    if FULL:
        # option to handle ADT bad limitations...
        write_this(name+"_ADT_DRY.pdb", dry_version)
        write_this(name+"_ADT_STRONG.pdb", strong_w)
        write_this(name+"_ADT_WEAK.pdb", weak_w)
else:
    print("[ WARNING: maps are ignored in multi-poses results]")
#exit()



if waters : #and ('-c' in opts or '-E' in opts):
    print("\n Predicted VS experimental \t2.0 A    \t   2.5 A")
    print(" ========================= \t-----    \t   -----")
    if map:
        if conserved_waters_s:
            calc_list = list(conserved_waters_s.keys())
            found_2, distance_log = confirm_predictions(calc_list, waters, cutoff = MATCH_DISTANCE)
            found_2_atoms = []
            for x in distance_log:
                found_2_atoms.append(x[1])
            found_25, distance_log = confirm_predictions(calc_list, waters, cutoff = MATCH_DISTANCE+.5)
            percentage_2 =  pc(found_2, len(calc_list))
            percentage_25 = pc(found_25, len(calc_list))
            print(" Predicted STRONG  [ %d ] :\t%02.2f%% (%d)\t %2.2f%% (%d)" % (len(calc_list), percentage_2, found_2, percentage_25, found_25))
            write_this(name+"_matching_waters_STRONG.pdb", found_2_atoms)

            if DEBUG: 
                found_3, distance_log = confirm_predictions(calc_list, waters, cutoff = 3.51)
                experimental_in_3A_range_STRONG = []
                for x in distance_log:
                    experimental_in_3A_range_STRONG.append(x[1])
                write_this(name+"_DEBUG_matching_waters_3Amax_range.pdb", experimental_in_3A_range_STRONG)

        if conserved_waters_w:
            calc_list = list(conserved_waters_w.keys())
            found_2, distance_log = confirm_predictions(calc_list, waters, cutoff = MATCH_DISTANCE)
            found_2_atoms = []
            for x in distance_log:
                found_2_atoms.append(x[1])
            found_25, distance_log = confirm_predictions(calc_list, waters, cutoff = MATCH_DISTANCE+.5)
            percentage_2 =  pc(found_2, len(calc_list))
            percentage_25 = pc(found_25, len(calc_list))
            print(" Predicted  WEAK   [ %d ] :\t%02.2f%% (%d)\t %2.2f%% (%d)" % (len(calc_list), percentage_2, found_2, percentage_25, found_25))
            write_this(name+"_matching_waters_WEAK.pdb", found_2_atoms)

            if DEBUG: 
                found_3, distance_log = confirm_predictions(calc_list, waters, cutoff = 3.51)
                experimental_in_3A_range_WEAK = []
                for x in distance_log:
                    experimental_in_3A_range_WEAK.append(x[1])
                write_this(name+"_DEBUG_matching_waters_3Amax_range_WEAK.pdb", experimental_in_3A_range_STRONG)
        print("\n", end=' ')

    else:
        calc_list = get_waters(get_lines(output_filename))
        for a in get_lines(output_filename):
            if a[0:4] == "ATOM" or a[0:6] == "HETATM":
                if pdbqt_atom_type(a) == "W":
                    calc_list.append(a)
    if False:
        found_2, distance_log = confirm_predictions(calc_list, waters, cutoff = MATCH_DISTANCE)
        found_25, distance_log = confirm_predictions(calc_list, waters, cutoff = MATCH_DISTANCE+.5)
        percentage_2 =  pc(found_2, len(calc_list))
        percentage_25 = pc(found_25, len(calc_list))
        if not QUIET: 
            print("Reality percentage [2.0 A ]: %2.2f [ %d ]" % (percentage_2, found_2))
            print("Reality percentage [2.5 A ]: %2.2f [ %d ] " % (percentage_25, found_25))
        else:
            print("found_2.0A:%d(%2.2f)," % (found_2, percentage_2), end=' ')
            print("found_2.5A:%d(%2.2f)"  % (found_25, percentage_25))


if DEBUG:
    rec = open("DEBUG_rec_overlap.pdb", 'w')
    lig = open("DEBUG_lig_overlap.pdb", 'w')
    lone_lone = open("DEBUG_lone.pdb", 'w')
    log_file = open("DEBUG_log_file.log", "w")
#    print "LIG_LIG"
    for i in lig_lig:
        print(i[:-1], file=lig)

#    print "LONE"
    for i in lone:
        print(i[:-1], file=lone_lone)

#    print "REC"
    for i in lig_rec:
        print(i[:-1], file=rec)

    for i in LOG:
        print(i[:-1], file=log_file)


exit(0)



### USELESS (REMOVE)
# waters overlapping ligand atoms
# waters overlapping receptor atoms
# waters hanging in the space

### USEFUL

# waters close to any receptor acceptor atoms (NA, OA, SA)

# clustering algorithm

# waters overlapping waters ( counter for how many water are overlapping
   # average water position?

