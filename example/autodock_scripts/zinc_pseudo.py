#!/usr/bin/env python

import math
from copy import deepcopy
import sys

def dist(a,b):
    return math.sqrt(sum([(a[i]-b[i])**2 for i in range(min(len(a),len(b)))]))

def angle(a,b,c):
    """Calculate the angle between three points. First point in the middle"""
    d12 = dist(a, b)
    d13 = dist(a, c)
    d23 = dist(b, c) 
    #round. To avoid things like 1.000000001
    angle = math.acos(round((d12**2 + d13**2 - d23**2)/(2*d12*d13),7))
    return angle

def angled(a,b,c):
    return angle(a,b,c)*180/math.pi

def dihedral(a,b,c,d):
    """ Calculate dihedral considering a in the beggining"""
    v1 = [b[i] - a[i] for i in range(3)]
    v2 = [c[i] - b[i] for i in range(3)]
    v3 = [d[i] - c[i] for i in range(3)]
    temp = [dist((0,0,0),v2) * v1[i] for i in range(3)]
    y = dotProd(temp ,crossProd(v2,v3))
    x = dotProd(crossProd(v1,v2),crossProd(v2,v3))
    rad = math.atan2(y,x)
    return rad*(180/math.pi) 

def dotProd(a,b):
    N = min(len(a),len(b))
    return sum([a[i] * b[i] for i in range(N)])

def crossProd(a,b):
    """Pretty self-explanatory, this function bakes cookies"""
    normal_vect = [
    a[1]*b[2] - a[2]*b[1],
    a[2]*b[0] - a[0]*b[2],
    a[0]*b[1] - a[1]*b[0]]
    return normal_vect

def rawVec(a,b):
    N = min(len(a),len(b))
    return [b[i]-a[i] for i in range(N)]

class PDBQT():
    def __init__(self, line):
        self._parse_common(line)     # used here (PDB) and in PDBQT
        self._parse_specific(line)   # differs in PDBQT

    def getline(self):
        txt = self._print_common()    # no \n; PDB + PDBQT
        txt += self._print_specific() # differs in PDBQT
        return txt

    def dist(self, a):
        return dist(self.getcoords(), a.getcoords())

    def angle(self, a, b):
        return angled(self.getcoords(), a.getcoords(), b.getcoords())

    def getcoords(self):
        return (self.x, self.y, self.z)

    def setcoords(self, coords):
        self.x, self.y, self.z = coords
    # blanks: [11:12], [20:21], [27:30], [66:76](pdb) or [66:68](pdbqt)

    def _parse_common(self, line):
        """Common to PDB and PDBQT formats"""
        self.keyword     = line      [ 0: 6]     # ATOM or HETATM
        self.serial      = int(line  [ 6:11])    # atom id
        #                            [11:12]
        self.name        = line      [12:16]     # atom name
        self.altLoc      = line      [16:17]     # Alternate location
        self.resName     = line      [17:20]     # Residue name
        #                            [20:21] 
        self.chain       = line      [21:22]     # chain
        self.resNum      = int(line  [22:26])    # Residue number
        self.icode       = line      [26:27]     # ???
        #                            [27:30]
        self.x           = float(line[30:38])    # X
        self.y           = float(line[38:46])    # Y
        self.z           = float(line[46:54])    # Z
        self.occupancy   = float(line[54:60])    # Occupancy
        self.bfact       = float(line[60:66])    # Temperature factor

    def _parse_specific(self, line):
        """ PDBQT characters [68:79] """
        self.charge      = float(line[68:76])   # Charge
        self.atype       = line      [77:79]    # Atom type
        self.atype = self.atype.strip().upper()
        self.atomnr = self.atype_to_atomnr(self.atype)

    def _print_common(self):
        """ Characters [0:68]"""
        linestr = ''
        linestr += '%6s' % (self.keyword)
        linestr += '%5d' % (self.serial)
        linestr += ' ' 
        linestr += '%4s' % (self.name)
        linestr += '%1s' % (self.altLoc) 
        linestr += '%3s' % (self.resName)
        linestr += ' ' 
        linestr += '%1s' % (self.chain)
        linestr += '%4d' % (self.resNum)
        linestr += '%1s' % (self.icode)
        linestr += ' ' * 3 
        linestr += '%8.3f' % (self.x)
        linestr += '%8.3f' % (self.y)
        linestr += '%8.3f' % (self.z)
        linestr += '%6.2f' % (self.occupancy)
        linestr += '%6.2f' % (self.bfact)
        return linestr

    def _print_specific(self):
        """ PDBQT characters [68:79] """
        linestr =  ' ' * 2                      # [66:68]
        linestr += '%8.3f' % (self.charge)      # [68:76]
        linestr += ' ' * 1                      # [76:77]
        linestr += '%2s' % (self.atype)       # [77:79]
        linestr += '\n'
        return linestr


    def dist(self, other):
        s = [(a-b)**2 for (a, b) in zip(self.getcoords(), other.getcoords())]
        return math.sqrt(sum(s))

    def isbound(self, other_atom, cut_off_percent = .5):
        """ Depends on atomic number """
        threshold = cut_off_percent * (
            .5 * self.atomnr_vdw(self.atomnr) +
            .5 * self.atomnr_vdw(other_atom.atomnr))
        is_bound = self.dist(other_atom) < threshold
        return is_bound

    def atomnr_vdw(self, atomnr):
        D = {1:2.0, 6:4.0, 7:3.5, 8:3.2, 9:3.1, 12:1.3, 15:4.2, 16:4.0, 17:4.1, 
        20:2.0, 25:1.3, 26:1.3, 30:1.5, 35:4.3, 53:4.7} 
        return D[atomnr]

    def atype_to_atomnr(self, atype):
        D = {'H':1, 'HD':1, 'HS':1, 'C':6, 'A':6, 'N':7, 'NA':7, 'NS':7, 
             'OA':8,'OS':8, 'F':9, 'MG':12, 'S':16, 'SA':16, 'CL':17, 
             'CA':20, 'MN':25, 'FE':26, 'ZN':30, 'BR':35, 'I':53, 'G':6, 
             'J':6, 'P':15, 'Z':6, 'GA':6, 'Q':6, 'TZ':-1}
        try:
            return D[atype.strip().upper()]
        except:
            sys.stderr.write(
    'unexpected atom type: %s (not in standard forcefield)\n' % atype)
            return -1 

def load_pdbqt(filename):
    """Creates a list of PDBQT atom objects"""
    atoms_list = []
    max_id = 0
    num_tz = 0
    non_atom_text = {}
    f = open(filename)
    counter = 0
    for line in f:
        if line.startswith('ATOM  ') or line.startswith('HETATM'):
            atom = PDBQT(line)
            if atom.atype == 'TZ':
                num_tz += 1
            else:
                counter += 1
                if atom.atype.upper() == 'ZN':
                    # set zinc charge to zero
                    atom.charge = 0.0
                atoms_list.append(atom)
                max_id = max(max_id, atom.serial)
        else:
            if counter not in non_atom_text:
                non_atom_text[counter] = []
            non_atom_text[counter].append(line)
    f.close()
    if counter in non_atom_text: # text after all atoms
        non_atom_text['last'] = non_atom_text.pop(counter)
    return atoms_list, num_tz, max_id, non_atom_text

def bruteNearbyAtoms(atomsList, atype = 'ZN', cutOff = 4.5):
    """Find atoms close to given atype in pdb/pdbqt"""
    metalsList = [a for a in atomsList if a.atype.upper() == atype]
    allNearbyLists = []
    for metal in metalsList:
        bht_indx = []
        bht_dist = []
        for (i, atom) in enumerate(atomsList):
            if atom.dist(metal) < cutOff:
                bht_indx.append(i)
                bht_dist.append(atom.dist(metal))
        idx = [i[0] for i in sorted(enumerate(bht_dist),
                 key = lambda x:x[1])]
        nearbyAtoms = []
        for i in idx:
            nearbyAtoms.append(atomsList[bht_indx[i]])
        allNearbyLists.append(nearbyAtoms)
    return allNearbyLists

class znShell():
    """Process closest atoms to get coordination sphere
        Init with parameters, and call PROCESS to the current init
"""
    def __init__(self, ZnAtom, cutOffDist = 2.5, carboxyExp = 0.5):
        self.zn = ZnAtom
        self.c = cutOffDist
        self.e = carboxyExp
        self.lig = []
        self.rec = []

    def buildShell(self, atoms):
        """start the analysis process"""
        # Connectivity stuff
        connect, n_conn = self._getBonds(atoms)

        # Find carboxy 
        carboxyIndx = self._getCarboxyOxyIndx(atoms, connect, n_conn) 
        Os = [(atoms[o], atoms[O])for i,o,O in carboxyIndx] # Oxygen atoms
        atoms = self._rmCarboxy(atoms, carboxyIndx)         # remove carboxy
        connect, n_conn = self._getBonds(atoms)             # re-do connect

        # Remove 1-3 connections and invalid elements
        reach3 = self._build13(connect)             # carboxyls disconnected
        atoms = self._selectBinders(atoms, reach3)

        # Finally add carboxy pseudos
        n = len(atoms)                              # n atoms before carboxy
        atoms += self._avgCarboxy(Os)               # append pseudos
        atoms = self._rm_too_far(atoms)             # rm if dist > self.c
        coop_indx = [n + i for i in range(len(atoms) - n)] # indx of coo
        coop_indx.reverse()

        return atoms, Os, coop_indx

    def set_carboxyExp(self, e):
        self.e = e
        [self.rec.pop(i) for i in self.recOmask]
        [self.lig.pop(i) for i in self.ligOmask]
        r = len(self.rec)
        l = len(self.lig)
        self.rec += self._avgCarboxy(self.recO)
        self.lig += self._avgCarboxy(self.ligO)
        self.rec = self._rm_too_far(self.rec)
        self.lig = self._rm_too_far(self.lig)
        self.recOmask = [r + i for i in range(len(self.rec) - r)]
        self.ligOmask = [l + i for i in range(len(self.lig) - l)]
        self.recOmask.reverse()
        self.ligOmask.reverse()

    def proc_rec(self, nearAtoms):
        self.start_rec = nearAtoms # copy object
        self.rec, self.recO, self.recOmask = self.buildShell(nearAtoms)


    def proc_lig(self, nearAtoms):
        self.start_lig = nearAtoms # copy object
        self.lig, self.ligO, self.ligOmask = self.buildShell(nearAtoms)

    def _filter_NOS(self, atoms):
        valid = [7, 8, 16]
        to_pop = [i for i, a in enumerate(atoms) if a.atomnr not in valid]
        to_pop.reverse()
        [atoms.pop(i) for i in to_pop]
        return atoms


    def _getBonds(self, atoms):
        n = len(atoms)
        n_conn = [0 for i in range(n)] # connections by atom
        connect = []
        for i in range(n):
            for j in range(i+1,n):
                if atoms[i].isbound(atoms[j]):
                    connect.append((i,j))
                    n_conn[i] += 1
                    n_conn[j] += 1
        return (connect, n_conn)

    def _getCarboxyOxyIndx(self,atoms,connect,n_conn):
        
        # Find C=O (forcing n_conn 1 for O should eliminate hydroxyls)
        co = []
        for i,j in connect:
            if (atoms[i].atype in ['O','OA'] and 
                atoms[j].atype in ['C','A'] and n_conn[i] == 1):
                co.append((j,i))
            if (atoms[j].atype in ['O','OA'] and 
                    atoms[i].atype in ['C','A'] and n_conn[j] == 1):
                co.append((i,j))
        # Finally find carboxyl 
        coo = []
        for i in range(len(co)):
            for j in range(i+1,len(co)):
                if co[i][0] == co[j][0]: 
                    twoOxygens = (co[i][1], co[j][1])
                    coo.append((co[i][0], min(twoOxygens), max(twoOxygens)))
        return coo

    def _rmCarboxy(self, nearAtoms, indxPairs):
        indx2rm = []
        for c,i,j in indxPairs:
            [indx2rm.append(x) for x in (c,i,j)]
        indx2rm = sorted(indx2rm)
        indx2rm.reverse()
        for i in indx2rm:
            nearAtoms.pop(i)
        return nearAtoms

    def _avgCarboxy(self, oxys):
        avgs = []
        if self.e < 0:
            out = []
            for oxypair in oxys:
                out.append(oxypair[0])
                out.append(oxypair[1])
            return out
        for oxypair in oxys:
            A = deepcopy(oxypair[0])
            A.x, A.y, A.z = self.avgCarboxyInner(
                    oxypair[0].getcoords(), 
                    self.zn.getcoords(), 
                    oxypair[1].getcoords())
            A.name = 'AVG'
            A.atype = 'OC'
            avgs.append(A)
        return avgs

    def avgCarboxyInner(self, o1, zn, o2):
        """Find a pseudo point that is representative of the two Oxygen atoms
            with a dist based approach. This will be used to calculate the
            angle between the carboxylate and other coordinating atoms"""
        d1 = dist(o1,zn)
        d2 = dist(o2,zn) 
        oo = dist(o1,o2)
        ratio = ((d2-d1)/oo)**self.e # 0 to 1, 1 beign o2 more distant from zn
        weight = (1-ratio)/2 # if ratio==0 w12 should be half way between 01-o2
        v12 = rawVec(o1,o2)
        w12 = [v12[i]*weight for i in range(3)]
        return tuple([o1[i]+w12[i] for i in range(3)])

    def _build13(self, connect):
        # Build 1-3 indx pairs
        connect13 = []
        for i,j in connect:
            for x,y in connect:
                if j == x:
                    connect13.append((min(i,y), max(i,y)))
                if j == y:
                    connect13.append((min(i,x), max(i,x)))
                if i == x:
                    connect13.append((min(j,y), max(j,y)))
                if i == y:
                    connect13.append((min(j,x), max(j,x)))
        # Prune
        connect13clean = []
        for i,j in connect13:
            if i != j:
                pair = (i,j)
                if pair not in connect13clean:
                    connect13clean.append(pair)
        connect13 = connect13clean

        # Build an All connections list to delete faster
        # This contains, for each atom, the indices of the
        # 1-2 and 1-3 connected atoms.
        # Then, if an atom is chosen to be bound (carbons/hydrogens wont)
        # the connected atoms will be removed.
        # This is robust for the cases a C would be closer then a 
        # 1-2 connected NA - although the C has a lower index (is closer)
        # that fact does not remove the NA, which is the actual binder
        reach3 = {}
        for i,j in connect+connect13:
            if i in reach3:
                reach3[i].append(j)
            else:
                reach3[i] = [j]
            if j in reach3:
                reach3[j].append(i)
            else:
                reach3[j] = [i]
        #print '***', reach3
        return reach3

    def _rm_too_far(self, atoms):
        # Wrap me  ------ remove far away...
        too_far = []
        for i,a in enumerate(atoms):    
            if a.dist(self.zn) > self.c:
                too_far.append(i)
        too_far.reverse()
        for i in too_far:
            atoms.pop(i)
        return atoms

    def _selectBinders(self, nearAtoms, reach3):
        """This function here is pretty cool.
        It goes up atom-index by atom-index,
        starting with the closest atoms to Zn,
        (they are sorted by distances), and if it is a
        valid zinc binder, like NA, it removes the 
        1-2 and 1-3 bound atoms for the list.
        The remaining atoms are the binders."""
        ValidAtoms = ['O','OA','NA','N','S','SA','MG','MN','ZN','CA','FE','CU']
        remove = []
        for i,atom in enumerate(nearAtoms):
            invalid = (atom.atype not in ValidAtoms)
            too_far = (atom.dist(self.zn) > self.c)
            not_removed = (i not in remove)
            connected = (i in reach3)   # reach3 has no lonely atoms
            #print i+1, atom.dist(self.zn),invalid, too_far, not_removed, connected
            if invalid or too_far:
                remove.append(i) 
            elif not_removed and connected: 
                [remove.append(j) for j in reach3[i] if j not in remove]
        # new atoms list with leftovers form remove
        binderAtoms = [atom for i,atom in enumerate(nearAtoms) 
                            if i not in remove] 
        return binderAtoms

    def tetrahedral_pseudo(self, d = 2.0):
        n = len(self.rec)
        zn = self.zn.getcoords()
        # some self.tetrhdrlAngDev() test in the future
        tetrhdrlAngDev = 3.23 # random value just for fun
        wise_limit = 12.0 # another randm value for fun :D
        if n==3 and tetrhdrlAngDev < wise_limit:
            pseudo = deepcopy(self.zn)
            pseudo.name = 'TZ'
            pseudo.atype = 'TZ' 
            a,b,c = [a.getcoords() for a in self.rec] # 
            w = 2*(int(dihedral(a,b,c,zn)) > 0)-1       # 1 or -1
            # normal vector
            nv = crossProd(rawVec(a,b), rawVec(b,c))    
            # Normalize length
            nvn = [w*d*nv[i]/dist(nv,(0,0,0)) for i in range(3)]
            x,y,z = [zn[i]+nvn[i] for i in range(3)]     # origin on zn
            pseudo.x = x
            pseudo.y = y
            pseudo.z = z
        else:
            pseudo = None
        return pseudo

    def getAngles(self, atoms):
        out = []
        N = len(atoms)
        #print N
        for i in range(N):
            for j in range(i+1, N):
                out.append(self.zn.angle(atoms[i], atoms[j]))
        return out

    #def getCarboxy_Atoms_others(self):
    #coo_oxys = [znobj.rec[i] for i in znobj.recOmask]
    #coo_oxys += [znobj.lig[i] for i in znobj.ligOmask]
    #others = [a for a in znobj.rec if a not in coo_oxys]
    #others += [a for a in znobj.lig if a not in coo_oxys]

    def recTetraDev(self):
        """ Fast function - could be easily wrapped outside"""
        angs = self.getAngles(self.rec)
        devs = [(ang - 109.5)**2 for ang in angs]
        return math.sqrt( (sum(devs)) / len(devs))

    def ligTZrmsd(self, d = 2):
        tz = self.tetrahedral_pseudo(d)
        if tz == None:
            return tz
        devs = [(a.dist(tz))**2 for a in self.lig]
        return math.sqrt( float(sum(devs)) / len(devs))

def about(): 
    print('Description:')
    print('   places tetrahedral zinc pseudo atoms (atom type = TZ)')
    print('   required by the forcefield "AutoDock4Zn"')
    print('   http://pubs.acs.org/doi/abs/10.1021/ci500209e')
    print('')
    print('   TZ pseudo atoms are placed next to Zn atoms that may become')
    print('     tetrahedraly coordinated after ligand binding.')
    print('   Also, the charge of Zn atoms is set to zero.')
    print('   Only ATOM and HETATM records are kept in the output')

def usage():
    print('Typical usage:')
    print('   python zinc_pseudo.py -r receptor.pdbqt')
    print('Options:')
    print('   -r   input receptor filename')
    print('   -o   output receptor filename | default = "input"_TZ.pdbqt')
    print('   -h   print this help message and exit') 
    print('   -a   print "about" message and exit')

def main():

    import getopt
    from os.path import splitext

    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], 'r:o:ah', ['help'])
    except getopt.GetoptError as err:
        print((str(err))) # will print something like "option -a not recognized"
        usage()
        sys.exit(2) 

    # Parse arguments
    for o, a in opts:
        if o in ['--help', '-h']:
            usage()
            sys.exit()
        if o == '-a':
            about()
            sys.exit()
        if o == '-r':
            input_name = a
        if o == '-o':
            output_name = a

    if not 'input_name' in locals():
        usage()
        sys.stderr.write('Error:\n    missing input receptor\n')
        sys.exit(2)

    if not 'output_name' in locals():
        (stem_name, extension) = splitext(input_name)
        output_name = stem_name + '_TZ' + extension

    # Fixed values
    distance = 2.0  # zinc to TZ pseudo atom
    carboxy = 0.5   # carboxylate averaging parameter (Supp. Info of paper)
    cutoff = 2.5    # to consider receptor atoms as zinc-coordinated

    # Load molecules 
    r, num_tz, lastserial, non_atom_text = load_pdbqt(input_name) 
    if num_tz:
        sys.stderr.write(
            'WARNING: ignoring TZ pseudo-atoms in %s\n' % input_name)

    # Get list of lists of atoms nearby each Zn atom
    atoms_by_zn = bruteNearbyAtoms(r)

    # Process Zn spheres
    znobjs = []
    tz_list = []
    for alist in atoms_by_zn: # for each list of atoms nearby a zn atom
        zn = alist[0]
        znobj = znShell(zn, cutoff, carboxy) 
        znobj.proc_rec(alist[1:]) # alist[0] is zn
        znobjs.append(znobj)
        tz_list.append(znobj.tetrahedral_pseudo(distance))

    # generate TZ pseudo atoms text
    tz_list = [tz for tz in tz_list if tz] # Removes tz == None
    for tz in tz_list:
        lastserial += 1
        tz.serial = lastserial
        r.append(tz)
    print(('Wrote %d TZ atoms on %s.' % (len(tz_list), output_name)))

    # write file
    recfile = open(output_name, 'w')
    for (counter, atom) in enumerate(r):
        if counter in non_atom_text:
            for line in non_atom_text[counter]:
                recfile.write(line)
        recfile.write(atom.getline())
    recfile.close()

main()

