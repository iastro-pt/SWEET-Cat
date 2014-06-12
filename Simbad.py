#!/usr/bin/python

# Build Simbad query to get info


result = open('tempSimbadQuery.txt', 'w')
result.write('echo Simbad script for planet hosts\n')

ttt = '|%IDLIST(HD|Gl|GJ|BD|HIP|CoroT|WASP|Kepler|KOI|KIC|HAT|NGC|XO|Qatar|'
ttt += 'TrES|OGLE|1)|%COO(A)|%COO(D)|%SP(S)|%FLUXLIST(V;F)|%FLUXLIST(V;E)|'
ttt += '%PLX(V)|%PLX(E)\"\n'


# #with open('Planethosts.rdb') as f:
# #   f.readline()
with open('NEWNEW1') as f:
    for line in f.readlines():
        line2 = line.replace('\n', '')
        line3 = 'format object form1 \"' + line2 + ttt
        result.write(line3.replace('\r', ''))
        result.write('query id '+line2.replace('-', ' ') + '\n')

result.write('format display')
result.close()
