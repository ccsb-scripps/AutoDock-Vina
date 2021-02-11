#!/usr/bin/env python

# Combine OA and HD maps for providing the W map
# for water molecules affinity
# The first map should be OA, the second HD

import getopt
from sys import argv, exit
import os
#from numpy import array

# VARIABLES

# weight of the affinity (to scale the interaction)
default_weight = .6
mode = "BEST"
ENTROPY=-0.2
#O_weight = .33
#hd_boost = 1.7
O_weight = 1.
hd_boost = 1.

#mode = "AVG"
#mode = "COOP"
#mode = "BOOST"

#if mode == "BOOST":
#    default_weight /= .05

#mode = "BEST"
#mode = "BOOST"
#mode = "BEST_W"

#ENTROPY=0.0 #TODO  <============= TEST THIS ON 2CCS to check if the penalty will lead to submerge the ligand inside the binding site
#ENTROPY=0.5 

O=0
H=0

#--OXYGEN
#O_weight = .33
#O_weight = .5
#O_weight = 1

#--HYDROGEN
#hd_boost = 1.8
#hd_boost = 1.7
#hd_boost = 1.0
#hd_boost = 0.6


def getLines(file, mode = "r"):
    fp = open(file, 'r')
    l = fp.readlines()
    fp.close()
    return l

def usage():
    #print "\n WET\n\n"
    myname =  os.path.basename(argv[0])

    print("""\n
                                        /\/\   __ _ _ __  
                                       /    \ / _` | '_ \ 
                                      / /\/\ \ (_| | |_) |           .-.
                                      \/    \/\__,_| .__/           (   )_
                                                   |_|              _/-'(_)
                                   __    __      _                 (_)
                                  / / /\ \ \__ _| |_ ___ _ __ 
                                  \ \/  \/ / _` | __/ _ \ '__|
                                   \  /\  / (_| | ||  __/ |   
                                    \/  \/ \__,_|\__\___|_|   
                                  
    """)

    print("\tUSAGE\n\t\t%s  -r | -o -h [options] \n\n" % myname)

    print("\tINPUT\n\t\tEither the receptor filename or the map files.\n")
    print("\t\t\t-r protein.pdbqt    :   the receptor filename is used to guess the OA and HD filenames")
    print("\n\t\t\t\t *OR*\n")
    print("\t\t\t-o protein.OA.map   :   OA map filename (omit if -r is used)")
    print("\t\t\t-h protein.HD.map   :   HD map filename (omit if -r is used)")

    print("\tOUTPUT\n\t\tW map.\n")
    #print "\t\tBy default the program saves the output by adding the suffix \"_HYDRO\" to the input PDBQT.\n"

    #print "\tOPTIONS"
    #print "\t\t-o\tIf \"output_filename\" is speficied, the output will be saved with this filename.\n\t\t\tFull path names can be used."

    #print "\t\t\tIf \"output_directory\" is provided, the output will be saved in the specified path by\n\t\t\tusing the default filename.\n\n"
    #print "\t\t-g\thydrate only the input atoms comprised in the volume specified by the GPF.\n\n"
    #print "\t\t-F\tThe output will be saved as \"input_HYDRO.pdbqt\" even if no waters are added."
    #print "\t\t\t(the default is to not to write the output if no waters have been added)"
    #print "\n\n"

    print("\tOPTIONS\n")
    print("\t\t-m (best,avg,coop..):   mix method (default 'best')")
    print("\t\t-e float            :   entropy penalty value (default: 0.0)")
    print("\t\t-O float            :   oxygen map weight (default: 1.0)")
    print("\t\t-H float            :   hydrogen map weight (default: 1.0)")
    print("\t\t-s protein.W.map    :   W map output filename (default: protein.W.map.[mix method].w[eight].o[xygen weight].h[ydrogen weight].E[ntropy correction])")
    print("\n\n")
    print("\n\tReference  ")
    print("\t----------")
    print("\t     Please cite:  Stefano Forli and Arthur J. Olson, J. Med. Chem., 2012, 55 (2), pp 623-638")
    print("\t                   DOI: 10.1021/jm2005145")
    print("\n\n") 

try:
    options, pdbqt = getopt.getopt(argv[1:], 'o:h:m:e:O:H:s:w:r:')
    # -r protein.pdbqt
    # -o protein.OA.map
    # -h protein.HD.map
    # -s protein.W.map
    # -w 0.5
    # -m best,avg,coop...
    # -e entropy
    # -O oxygen_weight
    # -H hydrogen_weight
except getopt.GetoptError as err:
 usage()
 print(str(err), "\n")
 exit(2)

opts = {}
for o, a in options: # populate the options
 opts[o] = a

try:
    opts
except:
    usage()
    exit(1)

print("ADD PWD AND FILE SUMMARY")


if '-r' in opts:
        protein =  os.path.basename(opts['-r'])
        name = protein.rsplit(".")[0]
        print("  receptor : ", protein)
        firstMap = name+".OA.map"
        secondMap = name+".HD.map"
        print("      OA map ->", firstMap)
        print("      HD map ->", secondMap)
        

else:
    if not "-o" in opts or not "-h" in opts:
        print("missing input1,2 or output")
        usage()
        exit(1)

    firstMap = opts['-o']
    secondMap = opts['-h']



if '-w' in opts:
    try:
        weight = float(opts['-w'])
    except:
        print("wrong weight")
        exit(1)
else:
    weight = default_weight
    print(" => Water map weight : DEFAULT [ %1.2f ]" % weight)


if '-m' in opts:
    if opts['-m'] == "best":
        mode = "BEST"    
    if opts['-m'] == "best_new":
        mode = "BEST_NEW"    
    if opts['-m'] == "optimal":
        mode = "OPTIMAL"    
    if opts['-m'] == "optimal_spherical":
        mode = "OPTIMAL_SPHERICAL"    
    if opts['-m'] == "avg":
        mode = "AVG"    
    if opts['-m'] == "coop":
        mode = "COOP"    
    if opts['-m'] == "boost":
        mode = "BOOST"    
    if opts['-m'] == "best_w":
        mode = "BEST_W"    


#else:
#    mode = "BEST"    

if mode == "BOOST":
    default_weight /= .05

if "-e" in opts:
    ENTROPY = float(opts['-e'])



if '-O' in opts:
    O_weight = float(opts['-O'])
if '-H' in opts:
    hd_boost = float(opts['-H'])


if '-s' in opts:
    outputFile = opts['-s']
else:
    # OLD for compatibility sake
    # outputFile = "protein.W.map.%s.w%1.2f.O%1.1f.H%1.1f.E%1.1f" % (mode, weight, O_weight, hd_boost, ENTROPY ) 
    outputFile = "protein.W.%s.w%1.2f.O%1.2f.H%1.2f.E%1.2f.map" % (mode, weight, O_weight, hd_boost, ENTROPY )
    # protein.W.map.BEST.w0.5.O0.3.H1.3.E0


print("\n  MapWater generator\n =====================")
print("  mode      : ", mode)
print("  weight    :  ", weight)
print("  HD_weight :  ", hd_boost)
print("  OA_weight :  ", O_weight)
print("  entropy   :  ", ENTROPY)



##################
# mixing modes


def best(first, second, positive = False ): #, hd_boost = hd_boost):
    global O, H
    """
    - Return the best value between two map values
    - If any positive values are found, 0 is returned by default
    - by default it should be first=OA, second=HD
    """
    
    first *= O_weight
    second *= hd_boost

    # TODO TODO TODO TODO
    # TODO POTENTIAL ERROR???
    # TODO - positive values for oxygen are misleading
    # TODO TODO TODO TODO

    if not positive and (first > 0 or second > 0): #
            # TODO add a smoothing function for positive values?
            return ENTROPY
    else:
        if first < second :
            O += 1
            return first * weight
        else:
            H += 1
            #print "First: %2.6f\tSecond: %2.6f" % (first, second)
            return second * weight


def BEST(oxygen, hydrogen):
    global O, H
    """
    - Return the best value between two map values
    - If any positive values are found, 0 is returned by default
    - by default it should be first=OA, second=HD
    """
    oxygen *= O_weight
    hydrogen *= hd_boost

    # TODO TODO TODO TODO
    # TODO POTENTIAL ERROR!!!
    # TODO - positive values for oxygen are misleading
    # ATTEMPT TO FIX...
    # TODO TODO TODO TODO

    if hydrogen > 0:  # HD is the smaller probe: if it doesn't fit, neither does OA
            return ENTROPY
    else:
        if oxygen < hydrogen :
            O += 1
            return oxygen * weight
        else:
            H += 1
            return hydrogen * weight



def optimal(oxygen, hydrogen):
    global O, H
    # NEVER TESTED
    oxygen *= O_weight
    hydrogen *= hd_boost

    #if (hydrogen < 0 and hydrogen > -0.2) and ():
    #    return 0


    if hydrogen >= 0:  # HD is the smaller probe: if it doesn't fit, neither does OA
            return ENTROPY
    else:
        if oxygen < 0 and hydrogen < 0 : # both values are negative: COMBINE
            O += 1
            H += 1
            return (oxygen + hydrogen) * weight
                                         # ..... or take the best
        elif oxygen < hydrogen:   
            O += 1
            return oxygen * weight
        elif hydrogen < oxygen:
            H += 1
            return hydrogen * weight


def optimal_spherical(oxygen, hydrogen):
    global O, H
    # NEVER TESTED
    oxygen *= O_weight
    hydrogen *= hd_boost

    oh_balance = 1

    #if (hydrogen < 0 and hydrogen > -0.2) and ():
    #    return 0


    if hydrogen >= 0:  # HD is the smaller probe: if it doesn't fit, neither does OA
            return ENTROPY
    else:
        if oxygen < 0 and hydrogen < 0 : # both values are negative: COMBINE
            O += 1
            H += 1
            #print "sphere!"
            return ((oxygen + hydrogen)*oh_balance) * weight
                                         # ..... or take the best
        elif oxygen < hydrogen:   
            O += 1
            return oxygen * weight
        elif hydrogen < oxygen:
            H += 1
            return hydrogen * weight


def cooperate(first, second, positive = False):
    global O, H
    """
    - Return the sum of value between two map values
    - If any positive values are found, 0 is returned by default
    """

    first *= O_weight
    second *= hd_boost

    if not positive and (first > 0 or second > 0):
            # TODO add a smoothing function for positive values?
            return ENTROPY
    else:
        value = (first + second) * weight
        O += 1
        H += 1
        # check
        if first > 0 or second >0:
            print("################ERROR!!!! ###################")

        #

        return value


def coop_blend():
    global O, H
    """
    - Return the sum of value between two map values
    - If any positive values are found, 0 is returned by default
    """

    first *= O_weight
    second *= hd_boost



    if not positive and (first > 0 or second > 0):
            # TODO add a smoothing function for positive values?
            return ENTROPY
    else:
        # TODO 
        # average of the total two contributions?
        value = ((first + second)/2) * weight
        O += 1
        H += 1
        # check
        if first > 0 or second >0:
            print("################ERROR!!!! ###################")

        #



def best_weighted(first, second, positive = False):
    global O, H
    # TODO test different parameters to check if there's
    # a real advantage
    """
    - Return the best value between two map values
    - If any positive values are found, 0 is returned by default
    - WEIGHTED: the OA map is scaled by weight to account for the higher
                intensity of OA maps versus HD
    """
    
    first *= O_weight
    second *= hd_boost

    if not positive and (first > 0 or second > 0):
            return 0 + ENTROPY
    else:
        if first < second :
            O += 1
            return first*weight
        else:
            H += 1
            return second*weight



def avg(first, second, positive = False):
    """
    - Return the average value between two map values
    - If any positive values are found, 0 is returned by default
    """
    if not positive and (first > 0 or second > 0):
            return ENTROPY
    else:
        value = ( (first + second)/2 ) * weight
        return value





def boost(first, second, positive = False):
    if not positive and (first > 0 or second > 0):
        return ENTROPY
    else:
        value =  ((first*7+second*7)/2) * weight
        return value





#try:
#    weight = float(argv[4])
#    print " => Water map weight :", weight
#except:
#    weight = default_weight
#    print " => Water map weight : DEFAULT [ %1.2f ]" % weight

#weight = 1
#weight = 4

#if not hd_boost == 1:
#    print " => HD boost ", hd_boost, "!"

#if not ENTROPY == 0:
#    print " => Entropy factor   :", ENTROPY
#else:
#    print "[ Entropy is ignored ]"

#print " => Calculated value :", mode

#print "\n\n"
#firstMap = argv [1]
#secondMap = argv [2]
#outputFile = argv [3]

data=[]
data2=[]
total=[]
HEADER=[]

try:
    #MAP1 = open( firstMap, 'r' )
    MAP1 = getLines(firstMap)
except:
    print("ERROR: the map %s can't be open" % firstMap)
    exit(1)

try:
    #MAP2 = open( secondMap, 'r' )
    MAP2 = getLines(secondMap)
except:
    print("ERROR: the map %s can't be open" % secondMap)
    exit(1)

try:
    DIFFERENCE = open ( outputFile, 'w' )
except:
    print("ERROR: impossible to open the output map", outputFile)
    exit(1)



HEADER = MAP1[0:6]
data = list(map(float, MAP1[6:]))
data2 = list(map(float, MAP2[6:]))



"""
for skip in "123456":
        HEADER.append(MAP1.next())


for line in MAP1:
     data.append(float(line))


for skip in "123456":
        MAP2.next()
for line in MAP2:
     data2.append(float(line))

"""


"""
MAP1 = open( firstMap, 'r' ).readlines()
MAP2 = open( secondMap, 'r' ).readlines()
DIFFERENCE = open ( outputFile, 'w' )


HEADER = MAP1[0:6]
data = array(MAP1[6:], float)
data2 = array(MAP2[6:], float)
"""

highest = 9999
lowest = -9999

if mode == "BEST":
    function = best
elif mode == "AVG":
    funcion = avg
elif mode == "COOP":
    function = cooperate
elif mode == "BOOST":
    function = boost
elif mode == "BEST_W":
    function = best_weighted
elif mode == "BEST_NEW":
    function = BEST
elif mode == "OPTIMAL":
    function = optimal
elif mode == "OPTIMAL_SPHERICAL":
    function = optimal_spherical

for i in range(len(data)):
    first=data[i]
    second=data2[i]
    value = function(first, second)

    if value > lowest:
        lowest = value
    if value < highest:
        highest = value
    total.append(value)

    """
    if mode == "BEST":
        value = best(first, second)
    elif mode == "AVG":
        value = avg(first, second)
    elif mode == "COOP":
        value = cooperate(first, second)
    elif mode == "BOOST":
        value = boost(first, second)
    elif mode == "BEST_W":
        value = best_weighted(first, second)
    elif mode == "BEST_NEW":
        value = BEST(first, second)
    elif mode == "OPTIMAL":
        value = optimal(first, second)
    elif mode == "OPTIMAL_SPHERICAL":
        value = optimal_spherical(first, second)
    """
        

#total = map(function, data, data2)

# output the results
for index in range(len(HEADER)):
    DIFFERENCE.write(str(HEADER[index]))

for index in range(len(total)):
    try:
        value = "%1.4f" % (total[index])
    except:
        print("problem")
    DIFFERENCE.write(value+'\n')
DIFFERENCE.close

O = float(O)
H = float(H)
tot = H+O
o_pc = (O/tot)*100
h_pc = (H/tot)*100

print("\n     Output info  ")
print("  --------------------")
print("  filename  :", outputFile)
print("  OA points : %3.2f%s" % ( o_pc , "%"))
print("  HD points : %3.2f%s\n" % ( h_pc, "%"))
#print O, H
print("  lowest  map value : %3.2f" % highest)
print("  highest map value : %3.2f" % lowest)
print("\n")

