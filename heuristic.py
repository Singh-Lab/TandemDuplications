from ete3 import Tree
from TreeUtils import run_multrec, parse_multrec, convert_trees, map_isometric_nodes
from ilp_generator import createDistMatrix, createTreeRepresentation, createEMatrix, \
                          createEqnsSingle, write
from TreeSearch import reconcileDL
import argparse
from gurobipy import *
import numpy as np
import sys
import copy

"""
Heuristic:
1) Run multrec
2) Run per node ILP to get dup sets
3) Move all dup sets to their LCA
4) Repeat?
"""

DEBUG = True

#Purely debugging stuff
def checkAbove(host, guest, realmap, infermap):
    above = 0
    for node in guest.traverse():
        real, infer = realmap[node], infermap[node]
        if host.get_common_ancestor(real, infer) == real and real != infer:
            above += 1
            print host.get_ascii()
            print real.name, infer.name
            print [realmap[i].name for i in node.children]

    print "Number of nodes above (should be 0):", above

#Util
def read_mapping(host, guest, mapfile):
	fullmap = {}
	f = list(open(mapfile))
	for line in f:
		hname, gname = line.strip().split('\t')
		try:
			hnode, gnode = host&hname, guest&gname
		except:
			continue
		fullmap[gnode] = hnode
	return fullmap

def read_dupfile(dupfile):
    dups = []
    f = list(open(dupfile))
    for line in f:
        dups.append(set([i.strip() for i in line.split('\t')]))
    return dups

def accuracy(real_tds, id_tds):
    if len(real_tds) == 0:
        return 1. if len(id_tds) == 0 else 0.
    correct = sum([i in real_tds for i in id_tds])
    return correct / (1. * len(real_tds))

#Util
def gen_map(guest, gtree):
	g2g = {}
	for leaf in gtree:
		gname = 'g' + leaf.name[1:]
		g2g[leaf] = guest&gname
	return g2g

#Util
def isolate_subtree(guest, hnode, mapping):
    nodes = [node for node in guest.traverse() if mapping[node] == hnode]
    roots = [node for node in nodes if (node.up == None or node.up not in nodes)]
    leaves = [node for node in nodes if (node.children == [] or node.children[0] not in nodes or node.children[1] not in nodes)]

    #Disconnect Everything, add fake root
    if len(roots) == 0:
        return None
    if len(roots) > 1:
        fakeroot = Tree()
        fakeroot.name = 'FAKE_ROOT'
        fakeroot.label = roots[0].up.label
        fakeroot.children = roots
        rootpars = [root.up for root in roots]
        for root in roots:
            root.up = fakeroot 
    else:
        fakeroot = roots[0]
        rootpars = [fakeroot.up]
        fakeroot.up = None
    leafchildren = [leaf.children for leaf in leaves]
    for leaf in leaves:
        if leaf.children == [] or leaf.children[0] not in nodes and leaf.children[1] not in nodes:
            leaf.children = []
        elif leaf.children[0] not in nodes:
            leaf.children = [leaf.children[1]] 
        else:
            leaf.children = [leaf.children[0]] 

    #Copy subtree
    subtree = fakeroot.copy()

    #Reconnect Everything
    for (i, root) in enumerate(roots):
        root.up = rootpars[i]
    for (i, leaf) in enumerate(leaves):
        leaf.children = leafchildren[i]

    return subtree

#Step 2: Extract tandem duplications from ILP
def identified_tds(m, guest):
    tds = []
    for item in m.getVars():
        if item.varName[0] == "X" and item.x == 1:
            a,b = item.varName[2:].split("_")
            a,b = guest.search_nodes(label=str(a))[0].name, guest.search_nodes(label=str(b))[0].name
            if a == 'FAKE_ROOT' or b == 'FAKE_ROOT':
                continue
            found = False
            for td in tds:
                if a in td or b in td:
                    td.add(a)
                    td.add(b)
                    found = True
            if not found:
                x = set()
                x.add(a)
                x.add(b)
                tds.append(x)

    return tds

#Step 2: Use td info on isomorphic copy tree to get tds ing guest
def build_dupsets(subtree, guest, tds):
    dups = []
    for td in tds:
        iso_td = [guest&node for node in td]
        dups.append(iso_td)

    return dups

#Step 2: Run ILP on guest subtree within a single host node
def run_ilp_single(subtree, ematrix):

    old_ordering = {node.name : int(node.label) for node in subtree.traverse()}

    #DEBUG CODE
    if DEBUG:
        valid_pairs = set()
        for a in subtree:
            for b in subtree:
                if ematrix[int(a.label)][int(b.label)] == 1:
                    valid_pairs.add((a.name, b.name))
                    valid_pairs.add((b.name, a.name))

    g = createTreeRepresentation(subtree)
    new_ordering = [x[1] for x in sorted([(int(i.label), i.name) for i in subtree.traverse()])]
    indices = [old_ordering[i] for i in new_ordering]
    new_ematrix = np.zeros((len(indices), len(indices)))
    
    for (i, io) in enumerate(indices):
        for (j, jo) in enumerate(indices):
            new_ematrix[i][j] = ematrix[io][jo]

    #DEBUG CODE
    if DEBUG:
        new_pairs = set()
        for a in subtree:
            for b in subtree:
                if new_ematrix[int(a.label)][int(b.label)] == 1:
                    new_pairs.add((a.name, b.name))
                    new_pairs.add((b.name, a.name))

        try:
            assert valid_pairs == new_pairs

        except:
            print 'valid pairs', valid_pairs
            print 'new pairs', new_pairs
            print subtree.get_ascii()
            print
            print subtree.get_ascii(attributes=['label'])
            print
            print 'Should be 1:', [a for a in valid_pairs if a not in new_pairs]
            print 'Should not be 1:', [a for a in new_pairs if a not in valid_pairs]
            sys.exit(1)

    eqnames, rhs, coldict = createEqnsSingle(subtree, new_ematrix, g)
    write('treesolve.mps', eqnames, rhs, coldict)
    #print "running ilp on tree with", len(subtree), 'nodes'
    m = read('treesolve.mps') # pylint: disable=undefined-variable
    m.Params.outputflag = 0
    m.optimize()
    return m, identified_tds(m, subtree)

def add_event_info(host, guest, mapping):

    for node in guest.traverse():
        if node.children == []:
            node.event = 'SPECIATION'
        else:
            ga, gb = node.children
            ha, hb = mapping[ga], mapping[gb]
            hnode = host.get_common_ancestor(ha, hb)
            if hnode in [ha, hb]:
                node.event = 'DUPLICATION'
            else:
                node.event = 'SPECIATION'

#Given a set of tds from the same host node, order them top to bottom 
def order_tds(tds):

    """
    def pairwise_td_order(td1, td2):
        #Returns 0,1,2 for td1 > td2, td2 > td1, incomparable
        thing = set(td2)
        for node in td1:
            x = node
            while x.up != None:
                if x.up in thing:
                    return 0

        thing = set(td1)
        for node in td2:
            x = node
            while x.up != None:
                if x.up in thing:
                    return 1

        return 2
    """

    def no_deps(td, unused):
        for node in td:
            if node.up in unused:
                return False

        return True

    all_dups = set()
    for td in tds:
        for d in td:
            all_dups.add(d)

    working_set = [i for i in tds]
    output = []

    while len(working_set) != 0:
        for (i, item) in enumerate(working_set):
            if no_deps(item, all_dups):
                output.append(item)
                for node in item:
                    all_dups.discard(node)
                working_set.pop(i)
                break

    return output

#Find the tds with no children in other tds in nodes
def bottom_tds(tds, nodes):
    bottoms = []

    for td in tds:
        bottom = True
        for node in td:
            if node.child in nodes:
                bottom = False
                break
        if bottom:
            bottoms.append(td)

    return bottom

#Find the tds with no parents in other tds in the input
def top_tds(tds, nodes):
    tops = []

    for td in tds:
        top = True
        for node in td:
            if node.up in nodes:
                top = False
                break
        if top:
            tops.append(td)

    return top

#Step 3: Identify nodes that should be moved down based on tandem duplication sets
def push_down(host, guest, mapping, td_dict):
    for hnode in [z for z in host.traverse()][::-1]:
        if hnode.children == []:
            continue

        tds = td_dict[hnode]
        all_dups = []
        for td in tds:
            all_dups += td

        all_dups = set(all_dups)

        for key in mapping:
            if mapping[key] == hnode and key not in all_dups:
                tds.append([key])

        for td in tds: #order_tds(tds)[::-1]:
            #Get LCA host node of each guest node in a td, then move the entire td to
            #the lca of those nodes.
            lcas = []

            for gnode in td:
                ga, gb = gnode.children
                ha, hb = mapping[ga], mapping[gb]
                lcas.append(host.get_common_ancestor(ha, hb))

            lca = host.get_common_ancestor(*lcas) if len(lcas) > 1 else lcas[0]

            for gnode in td:
                if mapping[gnode] != host&"h0" and lca == host&"h0":
                    print 'fuck', gnode.name
                mapping[gnode] = lca

    return mapping

def eligible(tdA, tdB, guest, ematrix):
    for nodeA in tdA:
        for nodeB in tdB:
            if ematrix[int(nodeA.label)][int(nodeB.label)] != 1:
                return False
            if guest.get_common_ancestor(nodeA, nodeB) in [nodeA, nodeB]:
                return False

    return True

def push_up(host, guest, mapping, td_dict, ematrix):
    for hnode in host.traverse():
        if hnode.up == None:
            continue

        tds = td_dict[hnode]

        full_mapping = set()
        for key in mapping:
            if mapping[key] == hnode:
                full_mapping.add(key)

        top = top_tds(tds, full_mapping)
        for td in top:
            parent_hnode = hnode.up
            parent_tds = td_dict[parent_hnode]
            for (i, ptd) in enumerate(parent_tds):
                if eligible(td, ptd, guest, ematrix):
                    tds.remove(td)
                    #add to new td
                    parent_tds[i] += td

def map_with_multrec(host, guest):
    s, g = convert_trees(host, guest)
    cost, gtree = run_multrec(s, g)

    host = Tree(host, format=1)
    guest = Tree(guest, format=1)
    #print 'guest immediate features', guest.features

    gtree, m_to_host = parse_multrec(gtree, host)
    m_to_g = map_isometric_nodes(guest, gtree, gen_map(guest, gtree))

    inferred_map = {}

    for key in m_to_host:
        gnode = m_to_g[key]
        hnode = m_to_host[key]
        inferred_map[gnode] = hnode

    return host, guest, inferred_map

def map_with_lca(host, guest):
    host = Tree(host, format=1)
    guest = Tree(guest, format=1)
    #print 'full tree has', len(host), 'host leaves and', len(guest), 'guest leaves'

    leafmap = {}
    for gnode in guest:
        hname = 'h' + gnode.name[1:].split("_")[0]
        hnode = host&hname
        leafmap[gnode] = hnode

    cost, fullmap = reconcileDL(host, guest, leafmap)
    return host, guest, fullmap

#host, guest are file paths, not Tree objects
def pipeline(host, guest):
    #Step 1: Run/Parse Multrec

    #host, guest, inferred_map = map_with_multrec(host, guest)
    host, guest, inferred_map = map_with_lca(host, guest)

    #Add in dup/spec information for the ILP
    add_event_info(host, guest, inferred_map)

    #Step 2: Run Per Node ILP
    grep = createTreeRepresentation(guest)
    ematrix = createEMatrix(guest)
    td_dict = {}
    for hnode in host.traverse():
        subtree = isolate_subtree(guest, hnode, inferred_map)
        if subtree == None:
            td_dict[hnode] = []
            continue
        m, iso_tds = run_ilp_single(subtree, ematrix)
        tds = build_dupsets(subtree, guest, iso_tds)
        td_dict[hnode] = tds

    intermediate_tds = []
    for key in td_dict.keys():
        item = td_dict[key]
        for td in item:
            intermediate_tds.append(td)

    
    #Step 3: Move Things to LCA
    mapping = push_down(host, guest, {key : inferred_map[key] for key in inferred_map}, td_dict)

    grep = createTreeRepresentation(guest)
    ematrix = createEMatrix(guest)
    td_dict = {}
    for hnode in host.traverse():
        subtree = isolate_subtree(guest, hnode, mapping)
        if subtree == None:
            td_dict[hnode] = []
            continue
        m, iso_tds = run_ilp_single(subtree, ematrix)
        tds = build_dupsets(subtree, guest, iso_tds)
        td_dict[hnode] = tds

    final_tds = []
    for key in td_dict.keys():
        item = td_dict[key]
        for td in item:
            final_tds.append(td)

    return host, guest, inferred_map, mapping, intermediate_tds, final_tds
    
    #return host, guest, inferred_map, intermediate_tds

if __name__ == '__main__':
    #Step 0: Input Arguments
    parser = argparse.ArgumentParser()

    parser.add_argument('host', type=str, help='The input host tree in newick format')
    parser.add_argument('guest', type=str, help='The input guest tree in newick format')
    parser.add_argument('mapping', type=str, help='Path to txt file containing host->guest leaf mapping')

    args = parser.parse_args()

    host = Tree(args.host, format=1)
    guest = Tree(args.guest, format=1)
    g_to_host = read_mapping(host, guest, args.mapping) #Ground Truth Mapping DO NOT USE IN INFERENCE
    dupfile = read_dupfile(args.host[:-8] + 'dupfile.txt')

    #Step 1: Run/Parse Multrec
    s, g = convert_trees(args.host, args.guest)
    cost, gtree = run_multrec(s, g)

    gtree, m_to_host = parse_multrec(gtree, host)
    m_to_g = map_isometric_nodes(guest, gtree, gen_map(guest, gtree))

    inferred_map = {}

    for key in m_to_host:
        gnode = m_to_g[key]
        hnode = m_to_host[key]
        inferred_map[gnode] = hnode

    if DEBUG:
        checkAbove(host, guest, g_to_host, inferred_map)

    #Add in dup/spec information for the ILP
    add_event_info(host, guest, inferred_map)

    #Step 2: Run Per Node ILP
    if DEBUG and False:
        print "\nFIRST RUN\n"
    grep = createTreeRepresentation(guest)
    ematrix = createEMatrix(guest)
    td_dict = {}
    for hnode in host.traverse():
        subtree = isolate_subtree(guest, hnode, inferred_map)
        m, iso_tds = run_ilp_single(subtree, ematrix)
        tds = build_dupsets(subtree, guest, iso_tds)
        td_dict[hnode] = tds

    if DEBUG and False:
        for key in td_dict:
            print key.name, [[j.name for j in i] for i in td_dict[key]]

    #Step 3: Move Things to LCA
    mapping = push_down(host, guest, {key : inferred_map[key] for key in inferred_map}, {key : td_dict[key] for key in td_dict})

    #Step 4: Reinfer tandem duplications with new mapping
    #TODO: After testing step 2, copy paste here
    if DEBUG and False:
        print "\nSECOND RUN\n"
    grep = createTreeRepresentation(guest)
    ematrix = createEMatrix(guest)
    td_dict = {}
    for hnode in host.traverse():
        subtree = isolate_subtree(guest, hnode, mapping)
        m, iso_tds = run_ilp_single(subtree, ematrix)
        tds = build_dupsets(subtree, guest, iso_tds)
        td_dict[hnode] = tds

    if DEBUG and False:
        for key in td_dict:
            print key.name, [[j.name for j in i] for i in td_dict[key]]
