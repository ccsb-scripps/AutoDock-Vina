#!/usr/bin/env python
#
# WET
#
# v.0.9  Stefano Forli
#
# Copyright 2011, Molecular Graphics Lab
#     The Scripps Research Institute
#        _  
#       (,)  T  h e
#      _/
#     (.)    S  c r i p p s
#      '\_
#       (,)  R  e s e a r c h
#      ./'
#     ( )    I  n s t i t u t e
#      "
#
#
# add water "atoms" in potentially hydrated
# regions of a given ligand
#
# 
# distance between water and ligand atom
space = 3.0 # Angstrom
bond_dist = 1.85 
pcharge = 0.000
residue = "WAT"
ATYPE = "W"

from numpy import *
from sys import argv, exit
import os 
import getopt

DEBUG = False
quiet = False

# TODO add multiple ligands input


# XXX Known bugs XXX
# - methoxy of 1BXO ligand recognized as SP2 (instead of SP3)
# - ether-O angle is not very precise


water_mates = []
GPF = False
numbering_stuff = []
FORCE=False
EXTENDED_ATOMS=False # phosphate_sulphate



def dist(f, s):  
    return sqrt((float(f[30:38])-float(s[30:38]))**2+(float(f[38:46])-float(s[38:46]))**2+ (float(f[46:54])-float(s[46:54]))**2)


def mean_pdb(firstline, secondline):
    # INFO   : calculate the mean point between two PDB atom lines
    # INPUT  : two pdb lines, an number for residue and atom numbering
    # OUTPUT : a pdb line

    coord1 = atom_coord(firstline)
    coord2 = atom_coord(secondline)

    x=(coord1[0]+coord2[0])/2
    y=(coord1[1]+coord2[1])/2
    z=(coord1[2]+coord2[2])/2
    atype=firstline[12:16]
    residue="MEA"
    chain="Y"
    count =1
    index = 1
    mean_atom="ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00          %1s" % (count, ATYPE, residue, chain, index, x, y, z, ATYPE)
    return mean_atom

def closest (first_atom, atom_list, cutoff=99999):
    # INFO   : find the closest atom 
    # INPUT  : a PDB atom, a list of PDB atoms [, cutoff distance]
    # OUTPUT : the closest atom [=null, if cutoff not satisfied] and short distance found
    # EXTRA  : dist function required
    best_distance=999999
    best_candidate=None
    for second_atom in atom_list:
            distance=dist(first_atom, second_atom)
            if distance < best_distance:
                best_distance=distance
                if best_distance < cutoff:
                    best_candidate=second_atom
    if best_candidate != "":
        return best_candidate, best_distance
    else:
        return best_candidate, best_distance

def rotatePoint(pt,m,ax):
    # From Ludo
    x=pt[0]
    y=pt[1]
    z=pt[2]
    u=ax[0]
    v=ax[1]
    w=ax[2]
    ux=u*x
    uy=u*y
    uz=u*z
    vx=v*x
    vy=v*y
    vz=v*z
    wx=w*x
    wy=w*y
    wz=w*z
    sa=sin(ax[3])
    ca=cos(ax[3])
    pt[0]=(u*(ux+vy+wz)+(x*(v*v+w*w)-u*(vy+wz))*ca+(-wy+vz)*sa)+ m[0]
    pt[1]=(v*(ux+vy+wz)+(y*(u*u+w*w)-v*(ux+wz))*ca+(wx-uz)*sa)+ m[1]
    pt[2]=(w*(ux+vy+wz)+(z*(u*u+v*v)-w*(ux+vy))*ca+(-vx+uy)*sa)+ m[2]
    return pt


def atom_coord(atom):
    coord = atom[28:56].split()
    for i in range(len(coord)):
        coord[i] = float(coord[i])
    return coord
    return [ float(f[30:38]), float(f[38:46]), float(f[46:54])]
     
def mean3(firstline, secondline, thirdline):
    # INFO   : calculate the mean point between two PDB atom lines
    # INPUT  : two pdb lines, an number for residue and atom numbering
    # OUTPUT : a pdb line
    coord1 = atom_coord(firstline)
    coord2 = atom_coord(secondline)
    coord3 = atom_coord(thirdline)

    x=(coord1[0]+coord2[0]+coord3[0])/3
    y=(coord1[1]+coord2[1]+coord3[1])/3
    z=(coord1[2]+coord2[2]+coord3[2])/3
    atype=firstline[12:16]
    residue="MEA"
    chain="Y"
    count =1
    index = 1
    mean_atom="ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00          %1s" % (count, ATYPE, residue, chain, index, x, y, z, ATYPE)
    return mean_atom


def hydro (atom1, atom2, atom3 = None, spacing = 3):
    # INFO   : place a water atom W at given distance from an atom
    # INPUT  : (a) two pdb lines for N-H; (b) three pdb lines for -[N|O]- acceptor, -O-H donors; the spacing distance between the water and the atom
    # OUTPUT : one pdb line 
    #
    # Synopsys:
    #
    # -C-NH-C  => atom1 = N, atom2 = H
    # -C-N=C   => atom1 = C, atom2 = N, atom3 = C
    # -C-O-H   => atom1 = C, atom2 = O, atom3 = H

    if atom2.split()[-1] == "HD":
        spacing -= 1 # to put the W at 3A from the N-H
    index = 99 
    coord2 = atom_coord(atom2)
    x=(coord2[0])
    y=(coord2[1])
    z=(coord2[2])
    chain="X"


    if atom3:
        atom4 = mean_pdb(atom1, atom3)
        vec_module = dist(atom2, atom4)
        coord1 = atom_coord(atom4) 
    else:
        coord1 = atom_coord(atom1)
        vec_module = dist(atom1, atom2)

    alpha = math.acos((coord2[0]-coord1[0])/vec_module)    # x-axis angle
    beta  = math.acos((coord2[1]-coord1[1])/vec_module)    # y-axis angle
    gamma = math.acos((coord2[2]-coord1[2])/vec_module)    # z-axis angle
    wat_x = spacing*math.cos(alpha)+x
    wat_y = spacing*math.cos(beta)+y
    wat_z = spacing*math.cos(gamma)+z

    wet = "%s%5d %2s   %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00     %1.3f %1s\n" % (keyw, index, ATYPE, residue, chain, index, wat_x, wat_y, wat_z, pcharge, ATYPE)
    return wet



def hydroH (atom1, atom2, atom3):
    
    print("middle")
    middle = hydro(atom1, atom2, spacing = 0)
    print("MIDDLE\n", middle)
    print("avg")
    avg = mean_pdb(middle, atom3)
    print("AVERAGE\n", avg)
    print("last")
    last = hydro(avg,atom2, spacing = space)
    print("LAST\n", last)



def vector(p1 , p2 = None):
    # accept both atoms and vectors as input
    #
    if type(p1) == type(str()):
        p1 = atom_coord(p1)
    x1 = p1[0]
    y1 = p1[1]
    z1 = p1[2]

    if type(p2) == type(str()):
        p2 = atom_coord(p2)
    
    if not p2 == None:
        x2 = p2[0]
        y2 = p2[1]
        z2 = p2[2]
        
        vec_x = x2-x1
        vec_y = y2-y1
        vec_z = z2-z1

        # it must be an array
        vec = array([vec_x, vec_y, vec_z], 'f')
        #print "REAL VECTOR", vec
    else:

        vec = array([p1[0], p1[1], p1[2] ], 'f' )
        #print "ATOM VECTOR", vec
    return vec

def norm(A):
        "Return vector norm"
        return sqrt(sum(A*A))

def normalize(A):
        "Normalize the Vector"
        return A/norm(A)



def bound(atom, structure, bond_dist = bond_dist, exclude = None):
    """
    identify all the atoms in "structure" that are @bond_dist from "atom"

    NOTE: this should be made with a lookup table!
    """
    bound_list = []
    tolerance = 0

    DEBUG = 0
    if DEBUG:
        print("Finding mates for ", atom)
    atype = atom.split()[-1]
    if atype == "HD":
        bond_dist = 1.15 #5 # previous bond dist of 1.1 cound't be enough for -S-H
    elif atype == "S" or atype == "SA":
        bond_dist = 1.95
    for candidate in structure:
        if candidate == atom or candidate == exclude:
            pass
        else:
            if candidate[0:4] == "ATOM" or candidate[0:6] == "HETATM":
                c_atype = candidate.split()[-1]
                if c_atype == "SA" or c_atype == "S" or c_atype  == "P":
                    if not atype == "HD":
                        tolerance = .20
                    else:
                        tolerance = .30
                elif c_atype == "HD":
                    tolerance = -.5
                else:
                    tolerance = 0

                #print candidate.strip(), tolerance+bond_dist
                if dist(atom, candidate) <= bond_dist + tolerance:
                    if not candidate in bound_list:
                        bound_list.append(candidate)
                        #print candidate,
                        #print tolerance+bond_dist
                else:
                    pass
    if len(bound_list) > 0:
        # Extra clean-up required for planar structures
        # where the HD lies between too many interested atoms...
        if atype == "HD":
            min = 999
            for b in bound_list:
                d = dist(atom,b)
                if d < min:
                    min = d
                    closest = b
            bound_list = [closest]
        return bound_list
    else:
        if not quiet:
            print("ERROR: this atom seems to be disconnected:", atom)
            print(atom)
        else:
            print("%s : error in atoms connections : %s" % (pdbqt, atom.strip()))
        exit(1)


def bound2(atom, structure, bond_dist = bond_dist, exclude = None):
    # from AD4_parameter.dat Rii/2 values
    cov_radii = { 'H': 1.00, 'HD': 1.00, 'HS': 1.00, 'C': 2.00,
            'A': 2.00, 'N': 1.75, 'NA': 1.75, 'NS': 1.75, 'OA': 1.60,
            'OS': 1.60, 'F': 1.54, 'Mg': 0.65, 'MG': 0.65, 'P': 2.10,
            'SA': 2.00, 'S': 2.00, 'Cl': 2.04, 'CL': 2.04, 'Ca': 0.99,
            'CA': 0.99, 'Mn': 0.65, 'MN': 0.65, 'Fe': 0.65, 'FE': 0.65,
            'Zn': 0.74, 'ZN': 0.74, 'Br': 2.165, 'BR':2.165, 'I':2.36,
            'Z' : 2.00, 'G' : 2.00, 'GA': 2.00, 'J' :2.00, 'Q' :2.00,
            'X': 2 } # default vdW for unknown atom
    bound_list = []
    atype = atom.split()[-1]
    r1 = cov_radii[atype]
    for candidate in structure:
        if candidate == atom or candidate == exclude:
            pass
        else:
            if candidate[0:4] == "ATOM" or candidate[0:6] == "HETATM":
                c_atype = candidate.split()[-1]
                r2 = cov_radii[atype]
                if dist(atom, candidate) <= r1+r2:
                    if not candidate in bound_list:
                        bound_list.append(candidate)
    if len(bound_list) > 0:
        # Extra clean-up required for planar structures
        # where the HD lies between too many interested atoms...
        if atype == "HD":
            min = 999
            for b in bound_list:
                d = dist(atom,b)
                if d < min:
                    min = d
                    closest = b
            bound_list = [closest]
        return bound_list
    else:
        if not quiet:
            print("ERROR: this atom seems to be disconnected:", atom)
            print(atom)
        else:
            print("%s : error in atoms connections : %s" % (pdbqt, atom.strip()))
        exit(1)

def calc_plane(atom1, atom2, atom3):
    DEBUG = 0
    # weird but it works...
    v12 = vector(atom1, atom2)
    v13 = vector(atom3, atom2)
    plane = cross(v12, v13)
    plane = normalize(plane)
    if DEBUG:
        print(atom1)
        print(atom2)
        print(atom3)
        print("PLANE FREE> coords: ", plane)
        print("PLANE FREE> type: ", type(plane))
        print("PLANE FREE> atoms:")
        keyw, index, ATYPE, residue, chain, pcharge = "ATOM  ", 1, "C", "RES", 1, 0.0000
        print("%s%5d %2s   %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00     %1.3f %1s" % (keyw, index, ATYPE, residue, chain, index, plane[0], plane[1], plane[2], pcharge, ATYPE))
        print("%s%5d %2s   %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00     %1.3f %1s\n" % (keyw, index, "X" , residue, chain, index, 0,0,0 , pcharge, "X"))

        centroid = mean3(atom1, atom2, atom3)
        arrow = vec_sum(vector(atom1, centroid), plane)
        print("%s%5d %2s   %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00     %1.3f %1s\n" % ("ATOM  ", 1, "A", 1, 1, 1, arrow[0], arrow[1], arrow[2], 0.000, "A"))
    return plane

def vec_sum(vec1, vec2):
    return array([vec1[0]+vec2[0], vec1[1]+vec2[1], vec1[2]+vec2[2] ], 'f')

def coplanar(plane, structure, reference, tolerance = .2):
    # identify all the atoms in "structure" that are co-planar with "plane"
    # the dot product between the point and the plane must be ~= 0
    
    coplane_list = []
    for atom in structure:
        position = vector(reference, atom)
        if dot(plane, position) <= tolerance:
            coplane_list.append(atom)
    return coplane_list


def dot(vector1, vector2):
    dot_product = 0.
    for i in range(0, len(vector1)):
        dot_product += (vector1[i] * vector2[i])
    return dot_product


def furanbolic(atom, structure, max = 2.35):
    # it should walk with pre-filtered co-planar atoms
    # HD's are automatically excluded
    the_ring = [atom]
    if atom.split()[-1] == "SA":
        max = 2.7
    for item in structure:
        if not item == atom:
            if not item.split()[-1] == "HD":
                if dist(atom, item) < max:
                    the_ring.append(item)
    if len(the_ring) == 5:
        if not quiet: print(" - possible furan/oxazole found...")
        return True
    if len(the_ring) > 6:
        if not quiet: print("WARNING: multiple atoms match the furan/oxazole check...")
        return True
    else:
        return False

def Osp2(oxygen, atom1, atom2):
    waters = []
    # hydroxyl/ether mode
    angles = [120, -120]
    #angles = range(0, 360, 10)
    oxyvector = vector(oxygen, atom1)
    oxyvector = normalize(oxyvector)
    for a in angles:
        roto = [oxyvector[0], oxyvector[1], oxyvector[2],  radians(a)]
        lone_pair_vector = vector(atom2, oxygen)
        lone_pair_vector = normalize(lone_pair_vector)
        water =  rotatePoint(-lone_pair_vector*space, atom_coord(oxygen), roto)  
        residue = "99"
        chain = "1"
        wet = "%s%5d %2s   %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00     %1.3f %1s\n"%(keyw, 1, ATYPE, residue, chain, 1, water[0], water[1], water[2], pcharge, ATYPE)
        waters.append(wet)
    return waters

def Osp2_NEW(oxygen, atom1, atom2):
    waters = []
    # hydroxyl/ether mode
    #angles = [120, -120]
    angles = list(range(0, 360, 10))
    oxyvector = vector(oxygen, atom1)
    oxyvector = vector(atom1, atom2)
    oxyvector = normalize(oxyvector)
    mid = mean_pdb(atom1, atom2)
    for a in angles:
        roto = [oxyvector[0], oxyvector[1], oxyvector[2],  radians(a)]
        lone_pair_vector = vector(mid, oxygen)
        lone_pair_vector = normalize(lone_pair_vector)
        water =  rotatePoint(+lone_pair_vector*space, atom_coord(oxygen), roto)  
        residue = "99"
        chain = "1"
        wet = "%s%5d %2s   %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00     %1.3f %1s\n"%(keyw, 1, ATYPE, residue, chain, 1, water[0], water[1], water[2], pcharge, ATYPE)
        waters.append(wet)
    return waters



def gpfminmax(gpf):
    # INPUT : gpf_file
    # OUTPUT: box coordinates (x,y,z), (X,Y,Z)
    """ Return max/min coordinates of
    the box described in the GPF"""
    if not GPF:
        return True
    #print gpf    
    file = open(gpf, 'r')
    lines = file.readlines()
    file.close()

    for line in lines:
        tmp=line.split()
        if tmp[0] == "gridcenter":
            center_x = float(tmp[1])
            center_y = float(tmp[2])
            center_z = float(tmp[3])
        if tmp[0] == "npts":
            pts_x = float(tmp[1])
            pts_y = float(tmp[2])
            pts_z = float(tmp[3])
        if tmp[0] == "spacing":
            res = float(tmp[1])

    step_x = pts_x/2 * res
    step_y = pts_y/2 * res
    step_z = pts_z/2 * res

    x_min = center_x - step_x
    x_max = center_x + step_x

    y_min = center_y - step_y
    y_max = center_y + step_y

    z_min = center_z - step_z
    z_max = center_z + step_z

    print(" - using the GPF box filter [ %s ]"% gpf)

    return [x_min, y_min, z_min], [x_max, y_max, z_max]


def in_the_box(atom_list, MIN, MAX):
    # INPUT : pdb-like atom, MIN = [x,y,z], MAX = [x,y,z]
    # OUTPUT: True/False

    good = []
    for atom in atom_list:
        pos = atom_coord(atom)
        if pos[0] < MAX[0] and pos[0] > MIN[0]:
            if pos[1] < MAX[1] and pos[1] > MIN[1]:
                if pos[2] < MAX[2] and pos[2] > MIN[2]:
                    good.append(atom)
    if good:
        return good
    else:
        return False




    
def usage():
    #print "\n WET\n\n"
    myname =  os.path.basename(argv[0])

    print("""\n
                                                      ,---,  
                                                   ,`--.' |  
                         .---.              ___    |   :  :  
                        /. ./|            ,--.'|_  '   '  ;  
                    .--'.  ' ;            |  | :,' |   |  |  
                   /__./ \ : |            :  : ' : '   :  ;  
               .--'.  '   \\' .   ,---.  .;__,'  /  |   |  '  
              /___/ \ |    ' '  /     \ |  |   |   '   :  |  
              ;   \  \;      : /    /  |:__,'| :   ;   |  ;  
               \   ;  `      |.    ' / |  '  : |__ `---'. |  
                .   \    .\  ;'   ;   /|  |  | '.'| `--..`;  
                 \   \   ' \ |'   |  / |  ;  :    ;.--,_     
                  :   '  |--" |   :    |  |  ,   / |    |`.  
                   \   \ ;     \   \  /    ---`-'  `-- -`, ; 
                    '---"       `----'               '---`"  
                                              

    """)

    print("\tUSAGE\n\t\t%s  [ -o output_filename | output_directory ] [ -g gpf.gpf ] [ -F ] -i input.pdbqt\n\n" % myname)

    print("\tINPUT\n\t\tThe \"input.pdbqt\" filename is required.\n")
    print("\tOUTPUT\n\t\tWet pdbqt file with W atoms. All atoms and tree-items (BRANCH,ENDBRANCH) are re-numbered accordingly.")
    print("\t\tBy default the program saves the output by adding the suffix \"_HYDRO\" to the input PDBQT.\n")

    print("\tOPTIONS")
    print("\t\t-o\tIf \"output_filename\" is speficied, the output will be saved with this filename.\n\t\t\tFull path names can be used.")

    print("\t\t\tIf \"output_directory\" is provided, the output will be saved in the specified path by\n\t\t\tusing the default filename.\n\n")
    print("\t\t-g\thydrate only the input atoms comprised in the volume specified by the GPF.\n\n")
    print("\t\t-F\tThe output will be saved as \"input_HYDRO.pdbqt\" even if no waters are added.")
    print("\t\t\t(the default is to not to write the output if no waters have been added)")
    
    print("\n\tReference  ")
    print("\t----------")
    print("\t     Please cite:  Stefano Forli and Arthur J. Olson, J. Med. Chem., 2012, 55 (2), pp 623-638")
    print("\t                   DOI: 10.1021/jm2005145")

    print("\n\n")

#####################################################################################
#    Beginning
#####################################################################################

try:
    options, extra = getopt.getopt(argv[1:], 'Fi:pho:g:q')
    # i = input ligand
    # d = output directory
    # o = output filename
    # g = reference GPF
    # p = do phosphate/sulphate
    # q = quiet
    # 
except getopt.GetoptError as err:
 usage()
 print(str(err), "\n")
 exit(2)

opts = {} # create the option dictonary with {option: argument} format

#print options
for o, a in options: # populate the options
 opts[o] = a

if '-h' in opts: # or '--help' in opts:
 usage()
 exit(0)

if '-i' in opts:
    pdbqt = opts['-i']

if '-q' in opts:
    quiet = True



try:
    input = open(pdbqt, 'r').readlines() # ugly, fix
    #name = pdbqt.rsplit(".")[0]
    name = os.path.splitext(pdbqt)[0]
    print("PDBQT", pdbqt)
    print("NAME is ",name)
    if not quiet:
        print("\n====================================")
        print("  ________         __   ")
        print(" |  |  |  |.-----.|  |_ ")
        print(" |  |  |  ||  -__||   _|")
        print(" |________||_____||____|\n")

        print("   processing %s"% pdbqt)
except:
    # print "\n ...I'm only asking for an input file... #
    print("\n\n\t# ERROR #\n\t Input filename is required.")
    usage()
    exit(1)


if '-o' in opts: # or '--type' in opts:
    output = opts['-o']
    #print "SPECIAL OPTION REQUIRED"
    #try:
    if os.path.isdir(output):
        if not quiet: print(" - saving the file in the path => %s" % output)
        name = output+os.path.sep+name+"_HYDRO.pdbqt"
    else:
        if not quiet: print(" - saving the file => %s" % output)
        name = output
#    except:
else:
    name += "_HYDRO.pdbqt"
    

if '-g' in opts:
    GPF = True
    MIN, MAX = gpfminmax(opts['-g'])


if '-F' in opts:
    # if no waters are added, return the input as "xxxx_HYDRO.pdbqt"
    # used to process multiple files
    FORCE = True


if '-p' in opts:
    if not quiet: print(" - include Phosphate/Sulphate groups")
    EXTENDED_ATOMS=True


####

if not EXTENDED_ATOMS:
    if not quiet: print(" - ignoring phosphate/sulphate groups")

hydrate_list = []
atoms_list = []

for line in input:
    if line[0:4] == "ATOM" or line[0:6] == "HETATM":
        atype = line.split()[-1]
        atoms_list.append(line)
        if atype == "OA" or atype == "NA" or atype == "HD" or atype == "SA":
            hydrate_list.append(line)
    if line[0:4] == "ATOM":
        keyw = "ATOM  "
    if line[0:6] == "HETATM":
        keyw = "HETATM"
    
if len(hydrate_list):
    if GPF:
        hydrate_list = in_the_box(hydrate_list, MIN, MAX)
    if not quiet: print(" - hydratable atoms : %d / %d " % (len(hydrate_list), len(atoms_list)))
else:
    if not FORCE:
        if not quiet: print(" [ No atoms to hydratate ]")
        exit(0)
    else:
        if not quiet: print(" [ No atoms to hydratate... FORCING TO SAVE OUTPUT...  ]")



numbering_stuff = []

# Scan the list to add waters
for atom in hydrate_list:
    atype = atom.split()[-1]
    HYDROXYL = False
    # ordinal position in the original file
    position = int(atom.split()[1])
    # add the atom to be hydrated to the buffer list
    waters_generated = [ atom ]
    # find the master atom(s)
    master = bound(atom, atoms_list)
    if len(master) == 0:
        print("\n\nERROR: this atom is disconnected:\n", atom)
        exit(1)
    # HYDROGENS #####################################################################################
    if atype == "HD":
        # check for errors
        if len(master) > 1:
            if not quiet:
                print("\n\nERROR (HD) : there is a proximity error and the following hydrogen is in close contact with more than one atom")
                print(atom)
                print("Bound mates:")
                for m in master:
                    print(m[:-1]," ==>", dist(m,atom))
            else:
                print("%s : HD proximity error" % (pdbqt))
            exit(1)
        else:
            # calculate the Water vector
            wet = hydro(master[0], atom)
            waters_generated.append(wet)
            water_mates.append(waters_generated)
        numbering_stuff.append([position, (len(waters_generated)-1)] )


    # OXYGENS #####################################################################################
    if atype == "OA" or atype == "SA":
        # OA includes the following options:
        #
        ## Two mates
        #  ---------
        #
        # -X-OA-X-    (ethers)
        #
        # -X-OA-HD    (hydroxyls)
        #
        # [-X-OA-X-]  (furan/pyran like)
        #
        #
        ## One mate
        # ----------
        #
        # -X=O        (carbonyl, carboxylate, nitro)   [ DONE]
        #
        #
        #  |
        # -X=O        (sulphonyl, phosphonyl) 
        #  |
        #

        if len(master) == 1:
            # # identify the mates of the master
            mates = bound(master[0], atoms_list, exclude = atom)
            # Phosphates check
            if len(mates) <=2:
                if DEBUG: print("[ carbonyl found ]")
                # 1. calculate the plane
                v12 = vector(atom, master[0])
                v23 = vector(mates[0], master[0])
                plane = cross(v12, v23)
                plane = normalize(plane)
                chain = 1 # TODO read this from the atom
                # O-lone pair 1
                roto = [ plane[0], plane[1], plane[2], radians(50) ]
                wat = rotatePoint(normalize(-v12)*space, atom_coord(atom),roto)
                wet = "%s%5d %2s   %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00     %1.3f %1s\n" %\
                            (keyw, 1, ATYPE, residue, chain, 1, wat[0], wat[1], wat[2], pcharge, ATYPE)
                waters_generated.append(wet)
                # O-lone pair 2
                roto = [ plane[0], plane[1], plane[2], radians(-50) ]
                wat = rotatePoint(normalize(-v12)*space, atom_coord(atom),roto)
                wet = "%s%5d %2s   %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00     %1.3f %1s\n" %\
                            (keyw, 1, ATYPE, residue, chain, 1, wat[0], wat[1], wat[2], pcharge, ATYPE)
                waters_generated.append(wet)
                water_mates.append(waters_generated)
                numbering_stuff.append([int(position), len(waters_generated)-1])
            else:
                if EXTENDED_ATOMS:
                    chain = 1
                    residue = 1 
                    directive = vector(master[0],atom)
                    directive = normalize(directive)*space
                    for q in mates:
                        position_q = vector(q) 
                        position_q = vec_sum(position_q, directive)
                        push = normalize(vector(atom, position_q))
                        start = vector(atom)
                        lpair = vec_sum(start, push*space)
                        wet = "%s%5d %2s   %3s %1s%4d    %8.3f%8.3f%8.3f  1.00 10.00     %1.3f %1s\n" %\
                            (keyw, 1, ATYPE, residue, chain, 1, lpair[0], lpair[1], lpair[2], pcharge, ATYPE)
                        waters_generated.append(wet)
                    water_mates.append(waters_generated)
                    numbering_stuff.append([int(position), len(waters_generated)-1])

        if len(master) == 2:
            for m in master:
                if m.split()[-1] == "HD":
                    HYDROXYL = True
                    if DEBUG: print(" [ >>> Found hydroxyl ]")
            if not HYDROXYL:
                # calculate the plane formed by oxygen and the two masters
                O_plane = calc_plane(atom, master[0], master[1])
                # get the list of all the co-planar atoms in the structure
                coplanar_mates = coplanar(O_plane, atoms_list, atom)
                # check if there are at least 4 coplanar atoms to make a ring with the OA
                if len(coplanar_mates) >=4 and furanbolic(atom, coplanar_mates):    
                    wet = hydro(master[0], atom, master[1])
                    waters_generated.append(wet)
                else:
                    lp_waters = Osp2(atom, master[0], master[1])
                    for w in lp_waters:
                        waters_generated.append(w)
            else:
                lp_waters = Osp2(atom, master[0], master[1])
                for w in lp_waters:
                    waters_generated.append(w)
            water_mates.append(waters_generated)
            numbering_stuff.append([int(position), len(waters_generated)-1])

    # NITROGEN #####################################################################################
    if atype == "NA":
        if len(master) == 1:
            # calculate the Water vector
            wet = hydro(master[0], atom)
            waters_generated.append(wet)
            water_mates.append(waters_generated)
            numbering_stuff.append([int(position), len(waters_generated)-1])

        # nitrile mode
        if len(master) == 2:
            wet = hydro(master[0], atom, master[1])
            waters_generated.append(wet)
            water_mates.append(waters_generated)
            numbering_stuff.append([int(position), len(waters_generated)-1])

        # tertiary amine HB receptor
        if len(master) == 3: 
            master_center = mean3(master[0], master[1], master[2])
            wet = hydro(master_center, atom)
            waters_generated.append(wet)
            water_mates.append(waters_generated)
            numbering_stuff.append([int(position), len(waters_generated)-1])

for mates in water_mates:
    index = input.index(mates[0])
    line = ""
    for atom in mates:
        line += atom
    input[index] = line
    

count = 1

# nicely split and format the lines that need to be splitted from
# the previous addition of water molecules
final = []
for line in input:
    line = line.split("\n")
    for item in line:
        if not item == "": final.append(item)

# renumbering the PDBQT atoms
for line in final:
    if line[0:4] == "ATOM" or line[0:6] == "HETATM":
        value = "%4s" % count
        idx = final.index(line)
        final[idx] = line[0:7]+value+line[11:]
        count += 1

# process the BRANCH numbers
for line in final:
    idx = final.index(line)
    if "BRANCH" in line:
        #print "Processing line => ", line
        line = line.split()
        value1, value2 = int(line[1]), int(line[2])
        #print "BEFORE =>", value1, value2 
        addendum1 = 0
        addendum2 = 0
        for mark in numbering_stuff:
            #print "Mark is :\t POSITION:", mark[0], " COUNT", mark[1]
            if value1 > mark[0]:
                addendum1 += mark[1]
                #print "value 1 is bigger than mark[0]:", value1, " | ", mark[0]
            if value2 > mark[0]:
                addendum2 += mark[1]
                #print "value 2 is bigger than mark[0]:", value2, " | ", mark[0]
        value1 += addendum1
        value2 += addendum2
        #print "AFTER  =>", value1, value2
        #print "- - - - - "
        final[idx] = line[0]+" "+str(value1)+" "+str(value2)

# Writing the output

try:
    hyd_ligand = open(name, 'w')
except:
    print("%s : # Error in saving the file %s #" % (pdbqt, name))
    exit(1)

count_waters = 0
for line in final:
    if line.split()[-1] == "W":
        count_waters += 1
    print(line, file=hyd_ligand)
if not quiet: print(" - %d waters added" % count_waters)

exit()
