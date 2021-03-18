from ete3 import Tree
from ilp_generator import *
from gurobipy import *
from tqdm import tqdm
from time import time
import numpy as np
from TreeUtils import rec_cost, reconcileDL
from heuristic import pipeline, identified_tds, isolate_subtree, add_event_info
import argparse

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

def map_accuracy(real_map, infer_map):
	total = len(real_map)
	correct = 0.0
	if type(list(infer_map.keys())[0]) == str:
		for key in real_map:
			if real_map[key].name == infer_map[key.name]:
				correct += 1

	else:
		for key in real_map:
			if real_map[key].name == infer_map[key].name:
				correct += 1

	return correct / total

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

def extract_mapping(m, host, guest):
	mapping = {}
	for item in m.getVars():
		if item.varName[0] == "M" and item.x == 1:
			a,b = item.varName[2:].split("_")
			a,b = guest.search_nodes(label=str(a))[0].name, host.search_nodes(label=str(b))[0].name
			mapping[a] = b

	return mapping
			
def str_to_map(strmap, host, guest):
	mapping = {}
	for key in strmap:
		mapping[guest&key] = host&(strmap[key])
	return mapping

def gen_map(guest, gtree):
	g2g = {}
	for leaf in gtree:
		gname = 'g' + leaf.name[1:]
		g2g[leaf] = guest&gname
	return g2g

def jaccard(a,b):
    if len(a.union(b)) == 0:
        return 1.
    return 1. * len(a.intersection(b)) / len(a.union(b))

def fdr(a,b):
	if len(a) == 0:
		return 0.
	incorrect = 0.0
	for item in a:
		if item not in b:
			incorrect += 1

	return incorrect / len(a)

def tdr(a,b):
	if len(b) == 0:
		return 0.
	correct = 0.0
	for item in a:
		if item in b:
			correct += 1
	
	return correct / len(b)

def recall(a,b):
	if len(a) == 0:
		return 1.
	correct = 0.0
	for item in a:
		if item in b:
			correct += 1
	
	return correct / len(a)

def real_td_list(dupfile):
    real_dups = read_dupfile(dupfile)
    real_pairs = set()
    for dup in real_dups:
        for a in dup:
            for b in dup:
                if True or a != b:
                    real_pairs.add(tuple(sorted((a, b))))
    return real_pairs

def td_to_pair_list(tds):
	allpairs = set()
	for td in tds:
		for a in td:
			for b in td:
				if True or a != b:
					allpairs.add(tuple(sorted((a.name, b.name))))
	return allpairs

def test_ilp(edist):
    acc = []
    map_acc = []
    real_cost = []
    ilp_cost = []
    totalTime = 0

    #[3,5,12,13,17,39,43,44,47,49]
    for ex in tqdm(range(50)):
        path = 'data/eventdist_' + str(edist) + '/' + str(ex) + '/'
        #path = 'otherdata/dupval_0.5/' + str(ex) + '/'
        
        #Read files
        real_dups = read_dupfile(path + 'dupfile.txt')
        host = Tree(path + 'host.nwk', format=1)
        guest = Tree(path + 'guest.nwk', format=1)

        #add in positions
        positions = list(open(path + 'posfile.txt'))
        for line in positions:
            name, pos = line.split('\t')
            pos = int(pos.strip())
            (guest&name).add_feature('pos', pos)

        #Read mapfile
        fullmapping = read_mapping(host, guest, path + 'mapfile.txt')
        lmapping = list(open(path + 'lmapfile.txt'))
        lmap = {}
        for line in lmapping:
            hname, gname = line.strip().split('\t')
            try:
                hnode, gnode = host&hname, guest&gname
            except:
                continue
            lmap[gnode] = hnode

        real_cost.append(rec_cost(guest, fullmapping, real_dups))

        #Reconciliation Portion
        h = createTreeRepresentation(host)
        g = createTreeRepresentation(guest)
        d = createDistMatrix(host)
        e = createEMatrix(guest)
        mapping = createMapping(lmap)

        eqnames, rhs, coldict = createEqns(host, guest, e, h, g, mapping, d)
        write('treesolve.mps', eqnames, rhs, coldict)
        m = read('treesolve.mps') #pylint: disable=undefined-variable
        m.setParam('OutputFlag', False)

        start = time()
        m.optimize()
        runtime = time() - start
        id_tds = identified_tds(m, guest)
        infer_map = extract_mapping(m, host, guest)
        map_acc.append(map_accuracy(fullmapping, infer_map))
        acc.append(accuracy(real_dups, id_tds))

        totalTime += runtime
        ilp_cost.append(rec_cost(guest, str_to_map(infer_map, host, guest), id_tds))

    print ("Average Time Per Example:", totalTime / 50.0)
    print ("Mapping Accuracy:", np.mean(map_acc))
    print ("Tandem Duplication Accuracy", np.mean(acc))

def test_heuristic(edist):
    outputs = []
    corrects = []
    totals = []
    multrec_cost = []
    heuristic_cost = []
    jaccards = []

    start = time()

    for ex in tqdm(range(50)):
        path = 'data/eventdist_' + str(edist) + '/' + str(ex) + '/'
        #path = 'otherdata_025/dupval_0.75/' + str(ex) + '/'
        host = path + 'host.nwk'
        guest = path + 'guest.nwk'
        fullmap = path + 'mapfile.txt'

        host, guest, intermediate_mapping, inferred_map, intermediate_tds, final_tds = pipeline(host, guest)
        true_map = read_mapping(host, guest, fullmap)

        correct = sum([1 if true_map[key] == intermediate_mapping[key] else 0 for key in true_map])
        total = len(true_map)
        corrects.append(correct)
        totals.append(total)

        output = str(correct) + ' / ' + str(total)
        output += " ! " if correct != total else ""
        outputs.append(output)
        multrec_cost.append(rec_cost(guest, intermediate_mapping, intermediate_tds, dtype='Tree'))
        heuristic_cost.append(rec_cost(guest, inferred_map, final_tds, dtype='Tree'))

        a = td_to_pair_list(final_tds)
        b = real_td_list(path + 'dupfile.txt')
        jaccards.append(recall(a,b))

    print ('Mapping Correctness:\n')
    for output in outputs:
        print (output )

    print ('Time per iteration:', (time() - start) / 50)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-type', 
                        type=str, 
                        help='The type of solver (exact/heuristic) to use', 
                        default='heuristic',
                        choices=['heuristic', 'exact'])

    parser.add_argument('-edist', 
                        type=str, 
                        help='The simulated event distance to test on (0.1/0.075/0.05/0.025/0.01)', 
                        default='0.1',
                        choices=['0.1','0.075','0.05','0.025','0.01'])

    args = parser.parse_args()

    if args.type == 'exact':
        test_ilp(args.edist)
    else:
        test_heuristic(args.edist)
