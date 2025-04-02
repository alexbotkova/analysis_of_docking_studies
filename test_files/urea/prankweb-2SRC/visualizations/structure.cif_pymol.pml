from pymol import cmd,stored

set depth_cue, 1
set fog_start, 0.4

set_color b_col, [36,36,85]
set_color t_col, [10,10,10]
set bg_rgb_bottom, b_col
set bg_rgb_top, t_col      
set bg_gradient

set  spec_power  =  200
set  spec_refl   =  0

load "data/structure.cif", protein
create ligands, protein and organic
select xlig, protein and organic
delete xlig

hide everything, all

color white, elem c
color bluewhite, protein
#show_as cartoon, protein
show surface, protein
#set transparency, 0.15

show sticks, ligands
set stick_color, magenta




# SAS points

load "data/structure.cif_points.pdb.gz", points
hide nonbonded, points
show nb_spheres, points
set sphere_scale, 0.2, points
cmd.spectrum("b", "green_red", selection="points", minimum=0, maximum=0.7)


stored.list=[]
cmd.iterate("(resn STP)","stored.list.append(resi)")    # read info about residues STP
lastSTP=stored.list[-1] # get the index of the last residue
hide lines, resn STP

cmd.select("rest", "resn STP and resi 0")

for my_index in range(1,int(lastSTP)+1): cmd.select("pocket"+str(my_index), "resn STP and resi "+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.show("spheres","pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_scale","0.4","pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_transparency","0.1","pocket"+str(my_index))



set_color pcol1 = [0.361,0.576,0.902]
select surf_pocket1, protein and id [1697,2452,2460,2552,2553,2554,2555,1696,1584,1689,1693,1567,2453,1535,1537,1538,1544,1545,1568,2066,2086,2476,2477,2087,2074,1583,1680,2057,1529,1525,1526,1527,1530,1531,1532,1533,1582,2090,2094,2095,2119,1917,1919,2041,2040,1936,2026,2025,1918,2564,2565,2566,2547,2548,2549,1548,2424,1554,1558,1561,1563,1565,1711,2423,2441,1552,1839,1840,2551,2575,2576,2577,2578,1783,1785,2579,2580,2608,2610,1747,2583,2655,2701,2654,2719,2720] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.365,0.278,0.702]
select surf_pocket2, protein and id [1331,1339,1358,1333,1335,61,101,74,81,98,100,78,82,1351,1948,1352,1357,431,67,72,1947,1663,1657,1661,1658,1662,2043,2051] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.792,0.361,0.902]
select surf_pocket3, protein and id [1885,1886,1887,2304,2328,3443,3444,2357,2358,1870,1879,1888,1892,576] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.702,0.278,0.533]
select surf_pocket4, protein and id [1319,543,545,547,1320,548,546,539,540,541,483,529,506,497,54,56,416,44,45,1324] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.902,0.361,0.361]
select surf_pocket5, protein and id [2110,2111,3008,2449,2954,2453,3007,2750,2749,2439,2718,2441] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.702,0.533,0.278]
select surf_pocket6, protein and id [2900,2741,2863,2864,2757,2832,2834,2407,2408,2409,2713,2719,2742,2418,2423] 
set surface_color,  pcol6, surf_pocket6 
   

deselect

orient
