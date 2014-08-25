#!/usr/bin/python


def simbad(stars):
    """Function to make script for Simbad. Takes a list of stars

    :stars: A list of stars
    :returns: An output file called SimbadQuery.txt for Simbad to read

    """

    t = ''
    for i, star in enumerate(stars):
        if i == 0:
            t += star
        else:
            t += '\n' + star
    t += '\n'
    f = open('NEWNEW1', 'w')
    f.write(t)
    f.close()

    result = open('SimbadQuery.txt', 'w')
    result.write('echo Simbad script for planet hosts\n')

    ttt = '|%IDLIST(HD|Gl|GJ|BD|HIP|CoroT|WASP|Kepler|KOI|KIC|HAT|NGC|XO'
    ttt += '|Qatar|TrES|OGLE|1)|%COO(A)|%COO(D)|%SP(S)|%FLUXLIST(V;F)|%'
    ttt += 'FLUXLIST(V;E)|%PLX(V)|%PLX(E)\"\n'

    with open('NEWNEW1') as f:
        for line in f.readlines():
            line2 = line.replace('\n', '')
            line3 = 'format object form1 \"' + line2 + ttt
            result.write(line3.replace('\r', ''))
            result.write('query id '+line2.replace('-', ' ') + '\n')

    result.write('format display')
    result.close()
