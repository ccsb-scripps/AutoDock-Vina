#!/usr/bin/env python
import sys
def oddify(i):
    return i + int(i % 2 == 0)
with open(sys.argv[1]) as f:
    parameters_read_count = 0 
    for line in f:
        if line.startswith("SPACING"):
            spacing = float(line.split()[1])
            parameters_read_count += 1
        if line.startswith("NELEMENTS"):
            nx, ny, nz = line.split()[1:]
            nx = oddify(int(nx))
            ny = oddify(int(ny))
            nz = oddify(int(nz))
            parameters_read_count += 1
        if line.startswith("CENTER"):
            cx, cy, cz = line.split()[1:]
            cx = float(cx)
            cy = float(cy)
            cz = float(cz)
            parameters_read_count += 1
        if parameters_read_count == 3:
            break
# vina uses 0.375 spacing and std::ceil to get the number of grid points
# here we subtract 0.1*spacing to guarantee that vina does not exceed
# the number of grid points due to floatint point resolution
print("size_x = %.3f" % (nx * spacing - spacing * 0.1))
print("size_y = %.3f" % (ny * spacing - spacing * 0.1))
print("size_z = %.3f" % (nz * spacing - spacing * 0.1))
print("center_x = %.3f" % cx)
print("center_y = %.3f" % cy)
print("center_z = %.3f" % cz)

