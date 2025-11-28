load ../nextcloud/ProximityReactions/final_Ecoli/AF2_outputs/P0AES2_P23522/ranked_1.pdb
set transparency, 0.5, All
hide everything, All
show surface, All
show wire, All
util.color_chains("(6r2g)",_self=cmd)
pseudoatom 1, pos=[30.865, 43.314, 30.424999999999997]
show spheres, 1
color pink, 1
pseudoatom 2, pos=[30.865, 43.314, 10.424999999999997]
show spheres, 2
color pink, 2
distance 0-1, 0, 1
color green, 0-1
distance 1-2, 1, 2
color green, 1-2

distance eucl,1,2
color red, eucl
