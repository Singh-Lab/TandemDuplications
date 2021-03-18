from ete3 import Tree
import numpy as np

def createTreeRepresentation(t):
    i = 0
    for node in [n for n in t.traverse()][::-1]:
        node.add_feature('label', str(i))
        #node.name = str(i)
        i += 1

    g = np.zeros((i,i))
    for node in [n for n in t.traverse()][::-1]:
        if node.children == []:
            pos = int(node.label)
            g[pos][pos] = 1
        else:
            #lchild = g[int(node.children[0].label)]
            #rchild = g[int(node.children[1].label)]
            children = [g[int(child.label)] for child in node.children]
            me = int(node.label)
            for k in range(len(g[0])):
                #g[me][k] = max(lchild[k], rchild[k])
                g[me][k] = max([child[k] for child in children])
            g[me][me] = 1

    return g

def createDistMatrix(host):
    size = len([node for node in host.traverse()])
    d = np.zeros((size, size), dtype=int)
    for a in host.traverse():
        for b in host.traverse():
            i, j = int(a.label), int(b.label)
            if i < j:
                d[i][j] = a.get_distance(b, topology_only=True)
            elif i > j:
                d[i][j] = b.get_distance(a, topology_only=True)
            else:
                d[i][j] = 0
    return d

def combine_dict(a,b):
    out = {key : a[key] for key in a}
    for key in b:
        if key in out:
            out[key] = min(out[key], b[key])
        else:
            out[key] = b[key]
    return out

def interleaved(i,j,k,l):
    """
    Assumes i, j are children or node a and k,l are children of node b
    """
    ikjl = i < k and i < j and i < l and k < j and k < l and j < l
    ilkj = i < l and i < k and i < j and l < k and l < j and k < j
    jkil = j < k and j < i and j < l and k < i and k < l and i < l
    jlik = j < l and j < i and j < k and l < i and l < k and i < k
    kilj = k < i and k < l and k < j and i < l and i < j and l < j
    kjli = k < j and k < l and k < i and j < l and j < i and l < i
    likj = l < i and l < k and l < j and i < k and i < j and k < j
    ljki = l < j and l < k and l < i and j < k and j < i and k < i

    return ikjl or ilkj or jkil or jlik or kilj or kjli or likj or ljki

def createEMatrix(guest):

    #Generate position dictionaries
    for node in [z for z in guest.traverse()][::-1]:
        if node.children == []:
            hname = 'h' + node.name.split("_")[0][1:]
            node.add_feature('posdict', {hname : node.pos})
        else:
            a,b = node.children[0].posdict, node.children[1].posdict
            node.add_feature('posdict', combine_dict(a,b))

    #Generate TD Eligibility Matrix
    size = len([node for node in guest.traverse()])
    e = np.zeros((size, size))

    for node1 in guest.traverse():
        for node2 in guest.traverse():
            if node1.children == [] or node2.children == []:
                continue
            if guest.get_common_ancestor(node1, node2) in [node1, node2]:
                continue

            i,j = node1.label, node2.label

            keys = set(node1.posdict.keys()).intersection(set(node2.posdict.keys()))
            child1, child2 = node1.children[0].posdict, node1.children[1].posdict
            child3, child4 = node2.children[0].posdict, node2.children[1].posdict
            keys = set(child1.keys()).intersection(child2.keys())
            keys = keys.intersection(set(child3.keys()))
            keys = keys.intersection(set(child4.keys()))

            if len(keys) == 0:
                e[int(i), int(j)] = 1
                continue

            for key in keys:
                x,y,z,w = child1[key], child2[key], child3[key], child4[key]
                if interleaved(x,y,z,w):
                    e[int(i), int(j)] = 1
                    continue

    return e

def createMapping(mapping):
    newmap = {}
    for node in mapping:
        newmap[int(node.label)] = int(mapping[node].label)
    return newmap

def addcolumn(var, coeff, eqname, coldict):
    entry = "   " + var + "\t" + eqname + '\t' + str(coeff)
    if var in coldict:
        coldict[var].append(entry)
    else:
        coldict[var] = [entry]

def addrhs(eqname, val, rhs):
    rhs.append("   RHS1\t" + eqname + '\t' + str(val))

def namevar(varname, ids):
    out = varname
    for i in ids:
        out += "_" + str(i)
    return out

#Generates ILP Equations for the full CDL Problem
def createEqns(h, g, e, hrep, grep, mapping, distmat, eqns='ALL'):
    eqnames = []
    rhs = []
    coldict = {}

    #COST
    eqnames.append(" N COST")

    for u in g.traverse():
        if u.children == []:
            addcolumn(namevar("X", [u.label, u.label]), 2, "COST", coldict)
        else:
            addcolumn(namevar("X", [u.label, u.label]), 4, "COST", coldict)
            v, w = u.children
            for i in range(len(hrep)):
                for j in range(len(hrep)):
                    addcolumn(namevar("Y", [u.label, v.label, i, j]), distmat[i][j]-1, "COST", coldict)
                    addcolumn(namevar("Y", [u.label, w.label, i, j]), distmat[i][j]-1, "COST", coldict)
        for k in range(1,7):
            #c = -1.5
            c = -2
            coeff = c * ((k - 1.) / k)
            addcolumn(namevar("T", [u.label, k]), coeff, "COST", coldict)

    #(0.0) - Mapping Equalities
    if eqns == "ALL" or 0 in eqns:
        eqcounter = 0
        for u in mapping:
            for i in range(len(hrep)):
                varname = namevar("M", [u, i])
                eqname = "INVAR" + str(eqcounter)
                eqnames.append(" E " + eqname)
                eqcounter += 1
                if mapping[u] == i:
                    addrhs(eqname, 1, rhs)
                else:
                    addrhs(eqname, 0, rhs)
                addcolumn(varname, 1, eqname, coldict)

    #(0.1) - Duplication Invariants
    if eqns == "ALL" or 1 in eqns:
        eqcounter = 0
        for u in g:
            varname = namevar("X", [u.label, u.label])
            eqname = "DUP_INVAR" + str(eqcounter)
            eqnames.append(" E " + eqname)
            eqcounter += 1
            addrhs(eqname, 0, rhs)
            addcolumn(varname, 1, eqname, coldict)

    #(2)
    if eqns == "ALL" or 2 in eqns:
        eqcounter = 0
        for u in g.traverse():
            if u.children == []:
                continue
            v, w = u.children
            for i in range(len(hrep)):
                for j in range(len(hrep)):
                    if hrep[i][j] == 0:

                        parent = namevar("M", [u.label, i])
                        child1 = namevar("M", [v.label, j])

                        eqname = "MAP" + str(eqcounter)
                        eqnames.append(" L " + eqname)
                        eqcounter += 1
                        addrhs(eqname, 1, rhs)
                        addcolumn(parent, 1, eqname, coldict)
                        addcolumn(child1, 1, eqname, coldict)

                        child2 = namevar("M", [w.label, j])

                        eqname = "MAP" + str(eqcounter)
                        eqnames.append(" L " + eqname)
                        eqcounter += 1
                        addrhs(eqname, 1, rhs)
                        addcolumn(parent, 1, eqname, coldict)
                        addcolumn(child2, 1, eqname, coldict)

    #(3)
    if eqns == "ALL" or 3 in eqns:
        eqcounter = 0
        for u in range(len(grep)):
            for v in range(len(grep)):
                if u == v:
                    continue
                eqname = "SMAP" + str(eqcounter)
                eqnames.append(" L " + eqname)
                eqcounter += 1
                addrhs(eqname, 0, rhs)
                addcolumn(namevar("X", [u, v]), 1, eqname, coldict)
                for i in range(len(hrep)):
                    addcolumn(namevar("Z", [u, v, i]), -1, eqname, coldict)

    #(4)
    if eqns == "ALL" or 4 in eqns:
        eqcounter = 0
        for u in range(len(grep)):
            for v in range(len(grep)):
                for i in range(len(hrep)):

                    #4.1
                    eqname = "SAMEMAP" + str(eqcounter)
                    eqnames.append(" G " + eqname)

                    addrhs(eqname, -1, rhs)
                    addcolumn(namevar("Z", [u, v, i]), 1, eqname, coldict)
                    addcolumn(namevar("M", [u, i]), -1, eqname, coldict)
                    addcolumn(namevar("M", [v, i]), -1, eqname, coldict)

                    #4.2
                    eqname = "SAMEUI" + str(eqcounter)
                    eqnames.append(" L " + eqname)

                    addrhs(eqname, 0, rhs)
                    addcolumn(namevar("Z", [u, v, i]), 1, eqname, coldict)
                    addcolumn(namevar("M", [u,i]), -1, eqname, coldict)

                    #4.3
                    eqname = "SAMEVI" + str(eqcounter)
                    eqnames.append(" L " + eqname)

                    addrhs(eqname, 0, rhs)
                    addcolumn(namevar("Z", [u, v, i]), 1, eqname, coldict)
                    addcolumn(namevar("M", [v, i]), -1, eqname, coldict)

                    eqcounter += 1

    #(5)
    if eqns == "ALL" or 5 in eqns:
        eqcounter = 0
        for u in range(len(grep)):
            for v in range(len(grep)):
                if u == v:
                    continue
                eqname = "NOANC" + str(eqcounter)
                eqnames.append(" L " + eqname)
                addrhs(eqname, 1 - grep[u][v], rhs)
                addcolumn(namevar("X", [u, v]), 1, eqname, coldict)

                eqcounter += 1

    #(6)
    if eqns == "ALL" or 6 in eqns:
        eqcounter = 0
        for u in range(len(grep)):
            for v1 in range(len(grep)):
                if u == v1:
                    continue
                for v2 in range(v1 + 1, len(grep)):
                    if grep[v1][v2] + grep[v2][v1] == 0:
                        continue
                    eqname = "LARGEDUP" + str(eqcounter)
                    eqnames.append(" L " + eqname)
                    addrhs(eqname, 1, rhs)
                    addcolumn(namevar("X", [u, v1]), 1, eqname, coldict)
                    addcolumn(namevar("X", [u, v2]), 1, eqname, coldict)

                    eqcounter += 1

    #(7)
    if eqns == "ALL" or 7 in eqns:
        eqcounter = 0

        def do6(u, v, a, b, eqcounter):
            if grep[u][a] * grep[v][b] == 0 and grep[a][u] * grep[v][b] == 0:
                return eqcounter
            eqname = "DUPCROSS" + str(eqcounter)
            eqnames.append(" L " + eqname)
            addrhs(eqname, 1, rhs)
            addcolumn(namevar("X", [u, v]), 1, eqname, coldict)
            addcolumn(namevar("X", [a, b]), 1, eqname, coldict)

            return eqcounter + 1

        for u in range(len(grep)):
            for v in range(u+1, len(grep)):
                for a in range(v+1, len(grep)):
                    for b in range(a+1, len(grep)):
                        eqcounter = do6(u,v,a,b,eqcounter)
                        eqcounter = do6(u,a,v,b,eqcounter)
                        eqcounter = do6(u,b,a,v,eqcounter)


    #(9)
    if eqns == "ALL" or 9 in eqns:
        eqcounter = 0
        for u in g.traverse():
            if u.children == []:
                continue
            v,w = u.children
            for i in range(len(hrep)):
                for j in range(len(hrep)):
                    eqname = "DUP" + str(eqcounter)
                    eqnames.append(" G " + eqname)
                    eqcounter += 1
                    addrhs(eqname, 0, rhs)
                    addcolumn(namevar("X", [u.label, u.label]), 1, eqname, coldict)
                    coeff = (hrep[i][j] + hrep[j][i]) / -2.0
                    varname = namevar("Y", [v.label, w.label, i, j])
                    addcolumn(varname, coeff, eqname, coldict)

    #(10)
    if eqns == "ALL" or 10 in eqns:

        def do9(u, v, i, j, eqcounter):

            #9.1
            eqname = "DUPAND" + str(eqcounter)
            eqnames.append(" G " + eqname)

            addrhs(eqname, -1, rhs)
            addcolumn(namevar("Y", [u.label, v.label, i, j]), 1, eqname, coldict)
            addcolumn(namevar("M", [u.label, i]), -1, eqname, coldict)
            addcolumn(namevar("M", [v.label, j]), -1, eqname, coldict)

            #9.2
            eqname = "DUPUI" + str(eqcounter)
            eqnames.append(" L " + eqname)

            addrhs(eqname, 0, rhs)
            addcolumn(namevar("Y", [u.label, v.label, i, j]), 1, eqname, coldict)
            addcolumn(namevar("M", [u.label, i]), -1, eqname, coldict)

            #9.3
            eqname = "DUPVJ" + str(eqcounter)
            eqnames.append(" L " + eqname)

            addrhs(eqname, 0, rhs)
            addcolumn(namevar("Y", [u.label, v.label, i, j]), 1, eqname, coldict)
            addcolumn(namevar("M", [v.label, j]), -1, eqname, coldict)

            return eqcounter + 1

        eqcounter = 0
        for u in g.traverse():
            if u.children == []:
                continue
            v, w = u.children
            for i in range(len(hrep)):
                for j in range(len(hrep)):
                    eqcounter = do9(u, v, i, j, eqcounter)
                    eqcounter = do9(u, w, i, j, eqcounter)
                    eqcounter = do9(v, w, i, j, eqcounter)
                    eqcounter = do9(w, v, i, j, eqcounter)

    #(11)
    if eqns == "ALL" or 11 in eqns:
        eqcounter = 0
        for i in range(len(grep)):
            eqname = "UMAP" + str(eqcounter)
            eqnames.append(" E " + eqname)
            eqcounter += 1
            addrhs(eqname, 1, rhs)
            for j in range(len(hrep)):
                addcolumn(namevar("M", [i, j]), 1, eqname, coldict)

    #(12)
    if eqns == "ALL" or 12 in eqns:
        eqcounter = 0
        for u in range(len(grep)):
            for v in range(len(grep)):
                if u == v:
                    continue
                eqname = "TDMAX" + str(eqcounter)
                eqnames.append(" L " + eqname)
                addrhs(eqname, 0, rhs)
                addcolumn(namevar("X", [u, v]), 1, eqname, coldict)
                addcolumn(namevar("X", [u, u]), -1, eqname, coldict)

                eqname = "TDEQUAL" + str(eqcounter)
                eqnames.append(" E " + eqname)
                addrhs(eqname, 0, rhs)
                addcolumn(namevar("X", [u, v]), 1, eqname, coldict)
                addcolumn(namevar("X", [v, u]), -1, eqname, coldict)

                eqcounter += 1

    #(13)
    if eqns == "ALL" or 13 in eqns:

        eqcounter = 0
        #(13.1)
        for u in range(len(grep)):
            eqname = "TDSIZE" + str(eqcounter)
            eqnames.append(" L " + eqname)
            addrhs(eqname, 1, rhs)
            for k in range(1,7):
                addcolumn(namevar("T", [u, k]), 1, eqname, coldict)
            eqcounter += 1

        #(13.2)
        eqcounter = 0
        for u in range(len(grep)):
            for k in range(1,7):
                eqname = "TDLEQ" + str(eqcounter)
                eqnames.append(" L " + eqname)
                addrhs(eqname, 0, rhs)
                addcolumn(namevar("T", [u, k]), 1, eqname, coldict)
                for v in range(len(grep)):
                    addcolumn(namevar("X", [u, v]), -1./k, eqname, coldict)
                eqcounter += 1

    #(14 - eligibility matrix) 
    if eqns == "ALL" or 14 in eqns:

        eqcounter = 0
        for u in range(len(grep)):
            for v in range(len(grep)):
                if u == v:
                    continue
                eqname = "ELIGIBLE" + str(eqcounter)
                eqnames.append(" L " + eqname)
                addrhs(eqname, int(e[u][v]), rhs)
                addcolumn(namevar("X", [u, v]), 1, eqname, coldict)
                eqcounter += 1

    return eqnames, rhs, coldict

#Generates ILP for a subtree within a single node
#Uses a subset of the equations in createEqns (removes mapping and host variables, and loss cost)
def createEqnsSingle(g, e, grep):
    eqnames = []
    rhs = []
    coldict = {}

    #COST
    eqnames.append(" N COST")

    for u in g.traverse():
        if u.children != []:
            addcolumn(namevar("X", [u.label, u.label]), 2, "COST", coldict)
        for k in range(1,8):
            c = -2
            coeff = c * ((k - 1.) / k)
            addcolumn(namevar("T", [u.label, k]), coeff, "COST", coldict)

    """
    #(0.1) - Duplication Invariants
    eqcounter = 0
    for u in g:
        varname = namevar("X", [u.label, u.label])
        eqname = "DUP_INVAR" + str(eqcounter)
        eqnames.append(" E " + eqname)
        eqcounter += 1
        addrhs(eqname, 0, rhs)
        addcolumn(varname, 1, eqname, coldict)
    """

    #(5)
    eqcounter = 0
    for u in range(len(grep)):
        for v in range(len(grep)):
            if u == v:
                continue
            eqname = "NOANC" + str(eqcounter)
            eqnames.append(" L " + eqname)
            addrhs(eqname, 1 - grep[u][v], rhs)
            addcolumn(namevar("X", [u, v]), 1, eqname, coldict)

            eqcounter += 1

    #(6)
    eqcounter = 0
    for u in range(len(grep)):
        for v1 in range(len(grep)):
            if u == v1:
                continue
            for v2 in range(v1 + 1, len(grep)):
                if grep[v1][v2] + grep[v2][v1] == 0:
                    continue
                eqname = "LARGEDUP" + str(eqcounter)
                eqnames.append(" L " + eqname)
                addrhs(eqname, 1, rhs)
                addcolumn(namevar("X", [u, v1]), 1, eqname, coldict)
                addcolumn(namevar("X", [u, v2]), 1, eqname, coldict)

                eqcounter += 1

    #(Triangle Inequality, see 13 in BCB Paper)
    eqcounter = 0
    for u in range(len(grep)):
        for v in range(len(grep)):
            for w in range(len(grep)):
                if u == v or u == w or v == w:
                    continue
                eqname = "TREQ" + str(eqcounter)
                eqnames.append(" G " + eqname)
                addrhs(eqname, -1, rhs)
                addcolumn(namevar("X", [u, v]), 1, eqname, coldict)
                addcolumn(namevar("X", [u, w]), -1, eqname, coldict)
                addcolumn(namevar("X", [v, w]), -1, eqname, coldict)

                eqcounter += 1

    #(7)
    eqcounter = 0

    def do6(u, v, a, b, eqcounter):
        if grep[u][a] * grep[v][b] == 0 and grep[a][u] * grep[v][b] == 0:
            return eqcounter
        eqname = "DUPCROSS" + str(eqcounter)
        eqnames.append(" L " + eqname)
        addrhs(eqname, 1, rhs)
        addcolumn(namevar("X", [u, v]), 1, eqname, coldict)
        addcolumn(namevar("X", [a, b]), 1, eqname, coldict)

        return eqcounter + 1

    for u in range(len(grep)):
        for v in range(u+1, len(grep)):
            for a in range(v+1, len(grep)):
                for b in range(a+1, len(grep)):
                    eqcounter = do6(u,v,a,b,eqcounter)
                    eqcounter = do6(u,a,v,b,eqcounter)
                    eqcounter = do6(u,b,a,v,eqcounter)

    #(Modified 8)
    eqcounter = 0
    for u in g.traverse():
        eqname = "DUPREQ" + str(eqcounter)
        eqnames.append(" EQ " + eqname)
        if u.name == 'FAKE_ROOT' or u.event == 'SPECIATION':
            addrhs(eqname, 0, rhs)
        else:
            addrhs(eqname, 1, rhs)
        addcolumn(namevar("X", [u.label, u.label]), 1, eqname, coldict)

        eqcounter += 1

    #(12)
    eqcounter = 0
    for u in range(len(grep)):
        for v in range(len(grep)):
            if u == v:
                continue
            eqname = "TDMAX" + str(eqcounter)
            eqnames.append(" L " + eqname)
            addrhs(eqname, 0, rhs)
            addcolumn(namevar("X", [u, v]), 1, eqname, coldict)
            addcolumn(namevar("X", [u, u]), -1, eqname, coldict)

            eqname = "TDEQUAL" + str(eqcounter)
            eqnames.append(" E " + eqname)
            addrhs(eqname, 0, rhs)
            addcolumn(namevar("X", [u, v]), 1, eqname, coldict)
            addcolumn(namevar("X", [v, u]), -1, eqname, coldict)

            eqcounter += 1

    #(13)
    eqcounter = 0
    #(13.1)
    for u in range(len(grep)):
        eqname = "TDSIZE" + str(eqcounter)
        eqnames.append(" L " + eqname)
        addrhs(eqname, 1, rhs)
        for k in range(1,8):
            addcolumn(namevar("T", [u, k]), 1, eqname, coldict)
        eqcounter += 1

    #(13.2)
    eqcounter = 0
    for u in range(len(grep)):
        for k in range(1,8):
            eqname = "TDLEQ" + str(eqcounter)
            eqnames.append(" L " + eqname)
            addrhs(eqname, 0, rhs)
            addcolumn(namevar("T", [u, k]), 1, eqname, coldict)
            for v in range(len(grep)):
                addcolumn(namevar("X", [u, v]), -1./k, eqname, coldict)
            eqcounter += 1

    #(14 - eligibility matrix) 
    eqcounter = 0
    for u in range(len(grep)):
        for v in range(len(grep)):
            if u == v:
                continue
            eqname = "ELIGIBLE" + str(eqcounter)
            eqnames.append(" L " + eqname)
            addrhs(eqname, int(e[u][v]), rhs)
            addcolumn(namevar("X", [u, v]), 1, eqname, coldict)
            eqcounter += 1

    return eqnames, rhs, coldict

def write(filename, eqnames, rhs, coldict):
    f = open(filename, 'w')
    f.write("NAME TREESOLVE\nROWS\n")

    for line in eqnames:
        f.write(line + '\n')

    f.write("COLUMNS\n   MARK0000\t'MARKER'\t'INTORG'\n")
    for var in sorted(coldict.keys()):
        for stmt in coldict[var]:
            f.write(stmt +'\n')
    f.write("   MARK0000\t'MARKER'\t'INTEND'\n")

    f.write("RHS\n")
    for line in rhs:
        f.write(line + '\n')

    f.write("ENDATA")
    f.close()

def extractMapping(m, host, guest):
    """
    Extracts a mapping from model m. Can only be called after m.optimize().
    Assumes mapping is stored in M_X_Y variables, finds those that are set to 1.

    Args:
        m (gurobipy model): The model object containin the ILP solution
    
    Outputs:
        mapping (dict): guest -> host mapping of all nodes in the guest tree
    """
    labelMapping = {}
    for item in m.getVars():
        if item.x == 1 and item.varName[0] == 'M':
            vals = item.varName.split("_")
            labelMapping[vals[1]] = vals[2]

    labelToHost = {}
    for node in host.traverse():
        labelToHost[node.label] = node

    mapping = {}
    for node in guest.traverse():
        hnode = labelToHost[labelMapping[node.label]]
        mapping[node] = hnode

    return mapping
