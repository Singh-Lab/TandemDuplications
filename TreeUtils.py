#Provides several utility functions useful for reading/writing trees, etc
from ete3 import Tree
from math import exp
from Utils import sortBy, printProgressBar
import numpy as np
import os
from sys import platform
#from ConfigParser import ConfigParser as CP

alphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 
                'T', 'V', 'W', 'Y']

#Tree I/O 
def readTree(filename):
    """
    Reads a tree in newick or NHX format to an ete3 object. To read trees not created by 
    this package, see readHostTree in HostTreeGen.py
    """
    return Tree(filename, format=1)

def writeReconciliation(host, guest, mapping):
    """Writes out a basic output version of the input trees and mapping without any
    information from the trees in nwk format, and writes a guest -> host mapping file.
    Automatically names all nodes that do not have a name.
    """
    os.system('mkdir output')

    i=0
    for node in host.traverse():
        if node.name == '':
            node.name = str(i)
            i += 1

    i=0
    for node in guest.traverse():
        if node.name == '':
            node.name = str(i)
            i += 1

    f = open('output/host.nwk','w')
    f.write(host.write(format=1)[:-1] + host.name + ';')
    f.close()

    g = open('output/guest.nwk','w')
    g.write(guest.write(format=1)[:-1] + guest.name + ';')
    g.close()

    mapfile = open('output/mapping.txt','w')
    for key in mapping:
        out = key.name + '\t' + mapping[key].name
        mapfile.write(out + '\n')
    mapfile.close()

def writeTree(tree, filename, features=[]):
    """Writes an ete3 tree to the given filename in NHX format"""
    output = tree.write(format=1, features=features)[:-1]

    #For whatever reason ete3 doesn't include root name + features, so this adds it on
    output += tree.name + "[&&NHX"
    for feature in tree.features:
        output += ":" + str(feature) + '=' + str(getattr(tree, feature))
    output += '];'

    outputHandle = open(filename, 'w')
    outputHandle.write(output)
    outputHandle.close()

def readMapping(host, guest, mapfile=None):
    """
    Requires that the host and guest tree that the mapping refers to have already 
    been read into memory. See writeMapping for mapfile format

    Args:
        host    (Tree): The host tree
        guest   (Tree): The guest tree
        mapfile (str ): Name of the file containing mapping between host and guest nodes
                        If None, mapping is inferred from host and guest node names

    Output:
        nodeMap (dict): A mapping of host -> [guest] nodes
    """

    if mapfile is not None:
        nodemap = {}
        for node in host.traverse():
            nodemap[node] = []

        nodemapFile = list(open(mapfile))

        for line in nodemapFile:
            line = line.split()
            hostNode = host&line[0]
            mapped = []
            for guestName in line[1:]:
                guestNode = guest&guestName
                mapped.append(guestNode)
            nodemap[hostNode] = mapped

        return nodemap

    else:
        pass

def genMap(host, guest, names=False):
    """
    Works only in the case that node gx_y maps to node hx
    """
    #{guest -> host}
    nodemap = {}
    for leaf in guest:
        gname = leaf.name
        hname = 'h' + gname.split("_")[0][1:]
        if names:
            nodemap[gname] = hname
        else:
            nodemap[leaf] = host&hname
    return nodemap

def writeMapping(nodemap, filename):
    """
    Writes out a mapping between host nodes and guest nodes. Each line of the output
    consists of a host node name followed by the names of each of the guest nodes 
    that maps to it (all separated by spaces).
    """
    outputHandle = open(filename, 'w')

    """
    for node in nodemap:
        out = node.name + '\t' + nodemap[node].name + '\n'
        #out = node.name + ' ' + ' '.join([i.name for i in nodemap[node]]) + '\n'
        outputHandle.write(out)
    """
    for hnode in nodemap:
        for gnode in nodemap[hnode]:
            out = hnode.name + '\t' + gnode.name + '\n'
            outputHandle.write(out)

    outputHandle.close()

def writeLeafMapping(nodemap, filename):
    """
    Writes out a mapping between host nodes and guest nodes for all leaf nodes in 
    the guest tree.
    """
    outputHandle = open(filename, 'w')

    for hnode in nodemap:
        for gnode in nodemap[hnode]:
            if gnode.children != []:
                continue
            out = hnode.name + '\t' + gnode.name + '\n'
            outputHandle.write(out)

    outputHandle.close()

#Duplication rate functions for use with guest tree generation
def s(x):
    """Sigmoid function designed to quickly reduce losses as domains are lost"""
    denom = 1 + exp(8-x)
    return .7 - .3 / denom

def s2(x):
    """Tighter sigmoid than s(x)"""
    denom = 1 + exp(10-x)
    return .7 - .6 / denom

def s3(x):
    """s2(x) but with strict cutoffs at 7 and 13 domains"""
    if x >= 13:
        return 0
    if x <= 7:
        return 1
    return s2(x)

def rec_cost(g, mapping, tds, dupcost=2, lcost=1, dtype='str'):

    duplications = set()
    for td in tds:
        if dtype == 'str':
            doms = set([g&i for i in td])
        else:
            doms = set(td)
        duplications = duplications.union(doms)

    def get_losses(a, b, ga, duplications):
        if a == b:
            return 0
        if a.name == 'g1_3':
            print (a in duplications)
        return int(a.get_distance(b, topology_only=True)) + 1 if ga in duplications else 0

    loss_cost = 0

    for node in g.traverse():
        if node.children == []:
            continue
        up, c1, c2 = mapping[node], mapping[node.children[0]], mapping[node.children[1]]
        loss_cost += lcost * (get_losses(up, c1, node, duplications) + get_losses(up, c2, node, duplications))

    dup_cost = dupcost * len(tds)

    return loss_cost + dup_cost

def reconcileDL(host, guest, leafmap):
    """
    Performs DL reconciliation on the guest tree and returns both the cost and the full
    mapping of guest -> host.

    Args:
    host (Tree): The host tree in ete3 format
    guest (Tree): The guest tree in ete3 format
    leafmapping (dict): guest -> host mapping containing all leaf nodes in the guest tree

    Output:
    cost (float): The cost of the min cost reconciliation with L = 1 and D = 2
    fullmap (dict): guest -> host mapping containing all nodes in the guest tree
    """
    fullmap = {}
    cost = 0
    losscost = 0
    seen = set()
    jobs = [leaf for leaf in guest]

    #Perform Reconciliation
    while jobs != []:
        job = jobs.pop(0)
        if job in leafmap:
            fullmap[job] = leafmap[job]
        else:
            a, b = job.children
            hostmap = fullmap[a].get_common_ancestor(fullmap[b])
            fullmap[job] = hostmap
        seen.add(job)
        guestParent = job.up
        if job.up == None:
            continue
        try:
            a, b = guestParent.children
        except:
            print (guest.get_ascii())
            raise Exception
        if a in seen and b in seen:
            jobs.append(guestParent)

    #Compute Cost
    for node in guest.traverse():
        #CASE 1: Leaf, no cost 
        if node.children == []:
            continue
        #CASE 2: Both children mapped to same node (DUPLICATION) easy duplication
        myhost = fullmap[node]
        lchild = fullmap[node.children[0]]
        rchild = fullmap[node.children[1]]

        if myhost == lchild and myhost == rchild:
            cost += 2
        #CASE 3: either lchild or rchild is mapped to the same but not the other (DUPLOSS)
        elif myhost == lchild:
            cost += 2 + myhost.get_distance(rchild, topology_only=True) + 1
            losscost += myhost.get_distance(rchild, topology_only=True) + 1
        elif myhost == rchild:
            cost += 2 + myhost.get_distance(lchild, topology_only=True) + 1
            losscost += myhost.get_distance(lchild, topology_only=True) + 1
        #CASE 4: Neither child is mapped to myhost (SPECIATION + LOSS)
        else:
            cost += myhost.get_distance(rchild, topology_only=True) + myhost.get_distance(lchild, topology_only=True)
            losscost += myhost.get_distance(rchild, topology_only=True) + myhost.get_distance(lchild, topology_only=True)

    return (cost - losscost) / 2, fullmap