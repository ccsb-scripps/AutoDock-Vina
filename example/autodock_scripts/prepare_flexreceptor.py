#!/usr/bin/env pythonsh
#
# 
#
# $Header: /mnt/raid/services/cvs/python/packages/share1.5/AutoDockTools/Utilities24/prepare_flexreceptor4.py,v 1.7.6.1 2015/08/26 22:45:31 sanner Exp $
#
import os 

from MolKit import Read
from MolKit.protein import ProteinSet, ResidueSet, AtomSet
from MolKit.molecule import BondSet
from MolKit.stringSelector import CompoundStringSelector

from AutoDockTools.MoleculePreparation import AD4FlexibleReceptorPreparation



if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        "Print helpful, accurate usage statement to stdout."
        print "Usage: prepare_flexreceptor4.py -r receptor_filename -s list_of_names_of_residues_to_move"
        print "    Description of command..."
        print "         -r     receptor_filename (.pdbqt)"
        print "         -s     specification for flex residues" 
        print "                Use underscores to separate residue names:"
        print "                  ARG8_ILE84  "
        print "                Use commas to separate 'full names' which uniquely identify residues:"
        print "                  hsg1:A:ARG8_ILE84,hsg1:B:THR4 "
        print "                [syntax is molname:chainid:resname]"
        print "    Optional parameters:"
        print "        [-v]    verbose output"
        print "        [-N]    type(s) of bonds to disallow: "
        print "        [-M]    interactive mode (automatic is default)"
        print "        [-P]    pairs of atom names bonds between which to disallow: hsg1:A:ARG8:CA_CB,CB_CG;hsg1:B:ARG8:CA_CB"
        print "        [-g pdbqt_filename] (rigid output filename)"
        print "        [-x pdbqt_filename] (flexible output filename)"

    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'r:vs:N:P:g:x:Mh')
    except getopt.GetoptError, msg:
        print 'prepare_flexreceptor4.py: %s' %msg
        usage()
        sys.exit(2)

    # initialize required parameters
    #-r: receptor
    receptor_filename =  None
    #-s: residues_to_move
    residues_to_move =  None
    # optional parameters
    verbose = None
    #-N: type of bonds to  disallow
    disallow = ""
    #-P: pairs of atom names bonds between which to disallow; residues separated by ';'
    all_disallowed_pairs = ""
    #-g  : rigid output filename
    rigid_filename = None
    #-x  : flexible output filename
    flexres_filename = None
    #-M  : mode
    mode = 'automatic'

    #'r:vs:N:g:x:Mh'
    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('-r', '--r'):
            receptor_filename = a
            if verbose: print 'set receptor_filename to ', a
        if o in ('-v', '--v'):
            verbose = True
            if verbose: print 'set verbose to', True
        if o in ('-s', '--s'):
            residues_to_move = a
            if verbose: print 'set residues_to_move to ', a
        if o in ('-N', '--N'):
            disallow = a
            if verbose: print 'set disallow to ', a
        if o in ('-P', '--P'):
            all_disallowed_pairs = a
            if verbose: print 'set all_disallowed_pairs to ', a
        if o in ('-g', '--g'):
            rigid_filename = a
            if verbose: print 'set rigid_filename to ', a
        if o in ('-x', '--x'):
            flexres_filename = a
            if verbose: print 'set flexres_filename to ', a
        if o in ('-M', '--M'):
            mode = 'interactive'
        if o in ('-h', '--'):
            usage()
            sys.exit()


    if not  receptor_filename:
        print 'prepare_flexreceptor4: receptor filename must be specified!\n'
        usage()
        sys.exit()

    if not  residues_to_move:
        print 'prepare_flexreceptor4: residues to move must be specified!\n'
        usage()
        sys.exit()

    extension = os.path.splitext(receptor_filename)[1]
    if extension!=".pdbqt":
        print 'prepare_flexreceptor4: receptor file must be in .pdbqt format\n'
        usage()
        sys.exit()

    rec = Read(receptor_filename)[0]
    bnds = rec.buildBondsByDistance()
    if verbose: print 'read ', receptor_filename

    all_res = ResidueSet()
    # hsg1:A:ARG8_ILE82;hsg1:B:THR4 
    # ARG8_ILE84 
    names_by_chain = residues_to_move.split(',')
    if verbose: 
        print "Specified flexres selection strings are:"
        for strN in names_by_chain: 
            print "   %s" %strN
    #1. ['hsg1:A:ARG8_ILE82','hsg1:B:THR4'] 
    # OR 
    #2. ARG8_ILE84
    for res_names_by_chain in names_by_chain:
        #check for ':'
        if verbose: print "selecting ", res_names_by_chain
        if res_names_by_chain.find(':')==-1:
            # no ':' in selection string = simple case eg: ARG8_ILE84 
            res_names = res_names_by_chain.split('_')
            if verbose: print "res_names=", res_names
            res = ResidueSet()
            for n in res_names:
                if verbose: print "looking for  ",n
                matched = rec.chains.residues.get(lambda x: x.name==n)
                #if matched.__class__ == AtomSet:
                if matched.__class__ == ResidueSet:
                    #res = matched.parent.uniq() #@@
                    res = matched.uniq() #@@
                if len(res):
                    for residue in res:
                        all_res += res
                        if verbose: print " added ", res.name, " to ", all_res
                else:
                    print "WARNING: no residue named " + n 
        else:
            # '1JFF_protein:A:THR179_GLU183,1JFF_protein:B:TRP21_ARG64_LEU275'
            chainStart = res_names_by_chain.index(':') + 1
            molName = res_names_by_chain[:chainStart-1]
            if verbose: print "molName = ", molName, " chainStart=", chainStart
            n = res_names_by_chain[chainStart:].replace("_", ",")
            if verbose: print "after comma replaced '_', n is ", n, '\n'
            selStr = molName + ":" + n
            if verbose: print "now selStr", selStr
            #eg: 'hsg1:A:ARG8,ILE82'
            result, msg = CompoundStringSelector().select(ProteinSet([rec]), selStr)
            if verbose: print "prepare_flexreceptor4 line 157: result=", result
            if result.__class__!= ResidueSet:
                print residues_to_move, " is not a ResidueSet instead it is a ", result.__class__
            #if verbose: print selStr, " selection =", res, " msg=", msg, '\n'
            resS = ResidueSet()
            if result.__class__ == AtomSet:
                resS = result.parent.uniq()
            else:
                resS = result
            if len(resS):
                all_res += resS
            else:
                print "no residue found using string ", selStr
    #if verbose: print "built all_res=", all_res.full_name()
    #check for duplicates
    d = {}
    for res in all_res: d[res] = 1
    all_res = d.keys()
    all_res = ResidueSet(all_res).uniq()
    all_res.sort()
    if verbose: 
        print "located ", len(all_res),  " residues to format:"
        for z in all_res: print "   %s" %z.full_name()

    all_bnds = BondSet()
    #inactivate bonds between specified atoms:
    #all_disallowed_pairs eg "1g9v_rec:A:ARG532:CA_CB,CB_CG,C_CA;1g9v_rec:B:ARG532:CA_CB,CB_CG"
    if len(all_disallowed_pairs):
        if verbose: print "line 183: all_disallowed_pairs = ", all_disallowed_pairs
        disallowed_pairs = all_disallowed_pairs.split(';') #';' between different residues
        for dp in disallowed_pairs: #1g9v_rec:A:ARG532:CA_CB,CB_CG,C_CA 
            #find the residue
            recN, chainN, resN, res_bnds = dp.split(':') #1g9v_rec, A, ARG532, CA_CB,CB_CB,C_CA
            if recN!=rec.name:
                print dp, " does not match ", rec.name
                exit()
            #result, msg = CompoundStringSelector().select(ProteinSet([rec]), dp)
            if verbose: print "line 199: res_bnds = ", res_bnds
            #r_bnds = res_bnds.split(',')
            #for bnd_pair in res_bnds: # find the residue then process the rest
            #    resK = bnd_pair
            #1g9v_rec:A:ARG532:CA_CB,CB_CG,C_CA   
            #CA_CB, CB_CG          
            #print "r_bnds = ", r_bnds
            for pair in res_bnds.split(','):
                if verbose: print "processing pair ", pair
                names = pair.split('_')
                if verbose: print "names = ", names
                bnds = all_res.atoms.bonds[0].get(lambda x: x.atom1.name in names and x.atom2.name in names)
                if verbose: print "bnds = ", bnds
                if len(bnds):
                    all_bnds += bnds
                    if verbose: 
                        print "line 207: added ",bnds, " to all_bnds"
                        print "all_bnds=", all_bnds
    if verbose:  print "beginning formatting ..."
    fdp = AD4FlexibleReceptorPreparation(rec, residues=all_res, rigid_filename=rigid_filename, 
                                            flexres_filename=flexres_filename,
                                            non_rotatable_bonds=all_bnds, mode = mode)

    if verbose:  print "finished!"

# To execute this command type:
# prepare_flexreceptor4.py -r receptor_filename -s list_of_names_of_residues_to_move"
# eg: in MGLToolsPckgs
#../bin/pythonsh AutoDockTools/Utilities24/prepare_flexreceptor4.py -r hsg1.pdbqt -s hsg1:A:ARG8_ILE84,hsg1:B:THR4 -P "hsg1:A:ARG8:CA_CB" -g hsg1_rigid.pdbqt -x hsg1_flex.pdbqt





