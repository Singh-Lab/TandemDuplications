from ete3 import Tree
from ilp_generator import *
from gurobipy import *
from tqdm import tqdm
from time import time
import numpy as np
from TreeUtils import rec_cost, reconcileDL
from heuristic import pipeline, identified_tds, isolate_subtree, add_event_info
import argparse

def create_mapping(host, guest):
    lmap = {}
    for gnode in guest:
        hname = 'h' + gnode.name[1:].split("_")[0]
        hnode = host&hname
        lmap[gnode] = hnode
    return lmap

def extract_mapping(m, host, guest):
	mapping = {}
	for item in m.getVars():
		if item.varName[0] == "M" and item.x == 1:
			a,b = item.varName[2:].split("_")
			a,b = guest.search_nodes(label=str(a))[0].name, host.search_nodes(label=str(b))[0].name
			mapping[a] = b

	return mapping

def run_ilp(gene_tree, domain_tree):

    host = Tree(gene_tree, format=1)
    guest = Tree(domain_tree, format=1)
    lmap = create_mapping(host, guest)

    h = createTreeRepresentation(host)
    g = createTreeRepresentation(guest)
    d = createDistMatrix(host)
    e = createEMatrix(guest)
    mapping = createMapping(lmap)

    eqnames, rhs, coldict = createEqns(host, guest, e, h, g, mapping, d)
    write('treesolve.mps', eqnames, rhs, coldict)
    m = read('treesolve.mps')

    m.optimize()
    id_tds = identified_tds(m, guest)
    infer_map = extract_mapping(m, host, guest)

    f = open('tandem_duplications.out','w')
    for td in id_tds:
        f.write(str(list(td)) + '\n')
    
    f.close()

    f = open('mapfile.out','w')
    for key in infer_map:
        f.write(key + '\t' + infer_map[key] + '\n')

    f.close()

def run_heuristic(gene_tree, domain_tree):

    host, guest, intermediate_mapping, inferred_map, intermediate_tds, final_tds = pipeline(gene_tree, domain_tree)
    
    f = open('tandem_duplications.out','w')
    for td in final_tds:
        f.write(str([i.name for i in td]) + '\n')
    
    f.close()

    f = open('mapfile.out','w')
    for key in inferred_map:
        f.write(key.name + '\t' + inferred_map[key].name + '\n')

    f.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-type', 
                        type=str, 
                        help='The type of solver (exact/heuristic) to use', 
                        default='heuristic',
                        choices=['heuristic', 'exact'])

    parser.add_argument('gtree',
                        type=str,
                        help='Path to newick file containing the gene tree')

    parser.add_argument('dtree',
                        type=str,
                        help='Path to newick file containing the domain tree')

    args = parser.parse_args()

    solver = run_ilp if args.type == 'exact' else run_heuristic

    solver(args.gtree, args.dtree)