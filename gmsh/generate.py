import csv

f = open('t.csv', 'r', encoding='utf-8')
c = csv.reader(f)
rows = list(c)
city = rows[-1]
rows = rows[:-1]
with open('script.geo', 'w', encoding='utf-8') as file:
    file.write("SetFactory(\"OpenCASCADE\");\n")
    for r in range(1, len(rows) + 1):
        file.write(f'Point({r}) = {{{rows[r-1][0]}, {rows[r-1][1]}, 0, 1.0}};\n')
    for r in range(1, len(rows)):
        file.write(f'Line({r}) = {{{r}, {r+1}}};\n')    
    file.write(f'Line({len(rows)}) = {{{len(rows)}, 1}};\n')        
        
    file.write(f"Curve Loop(1) = { set([i for i in range(1, len(rows) + 1)]) };\n")
    file.write("Plane Surface(1) = {1};\n")  
    file.write(f'Point({len(rows) + 1}) = {{{city[0]}, {city[1]}, 0, 0.001}};\n')
    file.write(f'Point{{{len(rows) + 1}}} In Surface{{1}};\n')
