#Provides several utility functions useful for reading/writing trees, etc
from ete3 import Tree
from math import exp
from Utils import sortBy, printProgressBar
import numpy as np
import os
from sys import platform
from ConfigParser import ConfigParser as CP

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

def writeFasta(names, sequences, filename, outgroup=False):
    """
    Writes out a fasta sequence to a file named <filename>. If outgroup is true,
    creates a fake outgroup (length 23 "AAAAAAAAAAAAAAAAAAAAAAA") and adds as the 
    last sequence in the fasta file
    """
    f = open(filename, 'w')

    #TODO: Make sure that this is actually a valid outgroup in all cases (should be)
    if outgroup:
        alphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', \
                        'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        fakeseq = ''.join(np.random.choice(alphabet, 23))
        f.write(">Outgroup\n" + fakeseq + '\n')

    for i in range(len(names)):
        f.write(">" + names[i] + '\n')
        f.write(sequences[i] + '\n')

    f.close()

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

#Utilities for finding domains

def findDomainsFile(infile, hmmfile, version=2):
    """
    Finds all domains in infile according to hmm in hmmfile using hmmsearch232.
    Assumes only one sequence is in the infile

    Args:
    infile  (str): The path to the fasta file to search for domains in
    hmmfile (str): Path to file containing hmm to search for

    Output:
    starts (list): List of starting positions of domains found
    ends   (list): List of ending positions of domains found
    seqs   (list): List of domain sequences found 
    """
    #domName = os.popen("grep 'NAME' " + hmmfile).read().split()[1]
    seqName = list(open(infile))[0][1:].strip()
    if version == 2:
        os.system("hmmsearch232 --domE 0.00001 " + hmmfile + " " + infile 
                    + " | grep '^ [ ]*" + seqName + "' > tmp/grepoutput.txt")
    else:
        os.system("hmmsearch --domE 0.00001 " + hmmfile + " " + infile 
                    + " | grep '^ [ ]*" + seqName + "' > tmp/grepoutput.txt")
    hits = list(open('tmp/grepoutput.txt'))
    starts, ends, seqs = [], [], []

    for line in hits:
        line = line.split()
        if not line[1].isdigit():
            continue
        starts.append(int(line[1]) - 1) #HMMER is 1-indexed :(
        ends.append(int(line[3]) - 1)
        seqs.append(line[2].upper())

    #HMMER doesn't sort by position in the sequence. This fixes that
    seqs = sortBy(seqs, starts)
    starts.sort()
    ends.sort()
        
    return starts, ends, seqs

def findMotifsFile(infile, mfile):
    """
    Finds all motif matches in infile according to motif in mfile using mast.
    Assumes only one sequence is in the infile

    Args:
    infile  (str): The path to the fasta file to search for domains in
    mfile (str): Path to file containing motif to search for

    Output:
    starts (list): List of starting positions of motifs found
    ends   (list): List of ending positions of motifs found
    seqs   (list): List of motif sequences found 
    """
    starts, ends, seqs = [], [], []
    sequence = list(open(infile))[1].strip()
    
    #Switch from MAST to FIMO
    os.system("fimo --text --thresh 1 " + mfile + " " + infile + " > tmp/mast_output.txt 2> crap.txt")
    f = list(open('tmp/mast_output.txt'))
    
    for line in f:
        if 'motif_id' not in line and line.strip() != '':
            temp = line.split('\t')
            starts.append(int(temp[3]) - 1) #MAST is 1-indexed :(
            ends.append(int(temp[4]) - 1)
            seqs.append(sequence[int(temp[3]) - 1 : int(temp[4])])


    return starts, ends, seqs

def findDomains(sequence, hmmfile, version=2):
    """
    Finds domains in the given sequence

    Args:
    sequence (str): sequence to search for domains in
    hmmfile (str): Path to file containing hmm to search for

    Output:
    starts (list): List of starting positions of domains found
    ends   (list): List of ending positions of domains found
    seqs   (list): List of domain sequences found 
    """
    g = open('tmp/tmp.fa','w')
    g.write('>seq\n' + sequence)
    g.close()
    
    return findDomainsFile('tmp/tmp.fa', hmmfile, version)

def findMotifs(sequence, mfile):
    """
    Finds motifs in the given sequence

    Args:
    sequence (str): sequence to search for motifs in
    hmmfile (str): Path to file containing hmm to search for

    Output:
    starts (list): List of starting positions of motifs found
    ends   (list): List of ending positions of motifs found
    seqs   (list): List of motif sequences found 
    """
    g = open('tmp/tmp.fa','w')
    g.write('>seq\n' + sequence)
    g.close()
    
    return findMotifsFile('tmp/tmp.fa', mfile)

def findSubsequences(sequence, mfile):
    if eval(CP.USE_HMMER): #pylint: disable=no-member
        return findDomains(sequence, mfile)
    return findMotifs(sequence, mfile)

def printDomSeq(sequence, hmmfile, minimal_rep = False):
    """
    prints the sequence with domains highlighted in red 
    (first character highlighted in green)

    Args:
    sequence (str): Protein sequence to print
    hmmfile  (str): Path to hmm file of domain to highlight
    mimimal_rep (bool): If true, prints a string of dashes and X's (nondomain and domain
                        sequences) rather than the full highlighted sequences
    """

    #Escape sequences used to colorize terminal output
    RED = '\033[91m'
    GREEN = '\033[92m'
    NORMAL = '\033[0m'

    #Find domains, check if sequence begins and/or ends with a domain
    domains = findSubsequences(sequence, hmmfile)[2]

    #split on all domains
    for domain in domains:
        sequence = sequence.replace(domain, "xxx")
    sequences = sequence.split("xxx")

    if minimal_rep:
        out = ''
        for i in range(len(domains)):
            out += '---' + RED + "XXX" + NORMAL
        out += '---'

        print out
        return
    
    #Reassemble full sequence post evolution
    out = ''
    for i in range(len(domains)):
        out += sequences[i] + GREEN + domains[i][0] + RED + domains[i][1:] + NORMAL
    out += sequences[-1] #if len(sequences) > len(domains) else RED + domains[-1] + NORMAL

    print out

def printDomSeq2(sequence, starts, ends):
    """
    prints the sequence with domains highlighted in red 
    (first character highlighted in green)
    Uses known domain start and end positions rather than finding them

    Args:
    sequence (str): Protein sequence to print
    hmmfile  (str): Path to hmm file of domain to highlight
    mimimal_rep (bool): If true, prints a string of dashes and X's (nondomain and domain
                        sequences) rather than the full highlighted sequences
    """

    #Escape sequences used to colorize terminal output
    RED = '\033[91m'
    GREEN = '\033[92m'
    NORMAL = '\033[0m'

    #Find domains, check if sequence begins and/or ends with a domain
    domains = [sequence[starts[i]:ends[i]+1] for i in range(len(starts))]

    #split on all domains
    for domain in domains:
        sequence = sequence.replace(domain, "xxx")
    sequences = sequence.split("xxx")
    
    #Reassemble full sequence post evolution
    out = ''
    for i in range(len(domains)):
        out += sequences[i] + '\n' + GREEN + domains[i][0] + RED + domains[i][1:] + NORMAL + '\n'
    out += sequences[-1] #if len(sequences) > len(domains) else RED + domains[-1] + NORMAL

    print out

def isValid(domain, hmmfile):
    """Checks if the input string is a valid zf-C2H2 domain"""

    #Can HMMER/MAST find it?
    if eval(CP.USE_HMMER): #pylint: disable=no-member
        seqs = findDomains(domain, hmmfile)[2]
        #Assumes that we are using zf if we are using an hmm instead of a motif
        valid = len(domain) == 23 and domain[2] == "C" and domain[5] == "C"
        valid &= domain[18] == "H" and domain[22] == "H"
    else:
        seqs = findMotifs(domain, hmmfile)[2]
        valid = True
    if len(seqs) == 0 or (len(seqs) > 0 and domain != seqs[0]):
        return False

    return valid

def raxml(infile, outext):
    """
    Runs raxml with set parameters on the input 

    Args:
    infile (str): fasta file to build tree from
    outext (str): output file extension to use. Output tree will be at
                  RAxML_bestTree.<outext>
    """

    #Remove previous run if it exists (RAxML will not clobber existing results)
    if "RAxML_bestTree.nwk" in os.listdir('.'):
        os.system('rm RAxML_*')

    command = 'raxml -s '
    command += infile + ' -n ' + outext + ' -m PROTGAMMAJTT -T 8 -p ' + str(np.random.randint(2000) + 1)
    command += ' > raxml_log.txt'
    os.system(command)

def raxml_score_from_file(benchfile, testfile, seqfile):
    """
    Uses RAxML to compute the SH score between the benchmark tree constructed by RAxML and
    a set of other trees. Outputs a binary vector where the ith entry is 1 if tree i is 
    significantly worse than the benchmark and 0 otherwise

    Args:
    benchfile (str): path to newick file containing the best tree found by RAxML
    testfile  (str): path to newick file containing all trees to test, one per line
    seqfile   (str): path to fasta file containing leaf sequences

    Output:
    scores (list): The list of ML scores given by raxml for each input tree
    worse  (list): A list of 0/1 entries specifying whether or not each tree is
                   significantly worse than the benchmark or not. 
    """
    #Run RAxML to find if tree in benchfile is significantly better than those in testfile
    #command = '/home/caluru/Downloads/standard-RAxML-master/raxmlHPC-PTHREADS-AVX2 '
    command = 'raxml '
    #Switch to -f h if this takes too long
    command += '-f H -t ' + benchfile + ' -z ' + testfile + ' -s ' + seqfile + ' -m PROTGAMMAJTT -T 8 -n sco'
    command += ' > raxml_log.txt'
    #TODO: Read results and select tree
    os.system(command)

    #Parse resulting logfile
    f = list(open('RAxML_info.sco'))
    start = 0
    for i in range(len(f)):
        if ' trees in File ' in f[i]:
            start = i

    f = f[start+3:]
    scores = []
    worse = []
    for line in f:
        if 'Tree: ' not in line:
            continue
        score = float(line.split()[3])
        scores.append(score)
        answer = line.split('Significantly Worse: ')[1].split()[0]
        worse.append(1 if answer == 'Yes' else 0)

    os.system('rm *.sco')
    return scores, worse

def raxml_score(benchTree, testTrees, seqfile):
    """
    Uses RAxML to compute the SH score between the benchmark tree constructed by RAxML and
    a set of other trees. Outputs a binary vector where the ith entry is 1 if tree i is 
    significantly worse than the benchmark and 0 otherwise. Takes in ete3 objects, writes
    them to the appropriate files and calls raxml_score_file

    Args:
    benchTree (Tree): The best tree found by raxml
    testTree  (list): List of trees to test against the best tree
    seqfile   (str): path to fasta file containing leaf sequences

    Output:
    scores (list): The list of ML scores given by raxml for each input tree
    worse  (list): A list of 0/1 entries specifying whether or not each tree is
                   significantly worse than the benchmark or not. 
    """
    g = open('bestTree.nwk', 'w')
    g.write(benchTree.write(format = 9) + '\n')
    g.close()

    g = open('otherTrees.nwk', 'w')
    for tree in testTrees:
        g.write(tree.write(format=9) + '\n')
    g.close()

    scores, worse = raxml_score_from_file('bestTree.nwk', 'otherTrees.nwk', seqfile)
    os.system('rm bestTree.nwk')
    os.system('rm otherTrees.nwk')

    return scores, worse

#TODO: only tested with Tree and dictionary inputs, not file inputs
def run_treefix(host, guest, lmap, sequences, short=True, suffix = ''):
    os.system('mkdir -p tfix' + suffix + '/config/; mkdir -p tfix' + suffix + '/data/0/')

    if type(host) == str:
        os.system('cp ' + host + ' tfix' + suffix + '/config/host.stree')
    else:
        writeTree(host, 'tfix' + suffix + '/config/host.stree')

    #guest is a path to a tree file
    if type(guest) == str:
        os.system('cp ' + guest + ' tfix' + suffix + '/data/0/0.nt.raxml.tree')
    #guest is a tree object
    else:
        writeTree(guest, 'tfix' + suffix + '/data/0/0.nt.raxml.tree')

    #lmap is a file
    if type(lmap) == str:
        os.system('cp ' + lmap + ' tfix' + suffix + '/config/guest.smap')
    #lmap is a dictionary of guest -> host nodes
    else:
        f = open('tfix' + suffix + '/config/guest.smap', 'w')
        for key in lmap:
            f.write(key.name + '\t' + lmap[key].name + '\n')
        f.close()

    #sequences is always a path to a fasta file
    os.system('cp ' + sequences + ' tfix' + suffix + '/data/0/0.nt.align')

    #Run TreeFix
    cmd = 'treefix -s tfix' + suffix + '/config/host.stree'
    cmd += ' -S tfix' + suffix + '/config/guest.smap' 
    cmd += ' -A nt.align '
    cmd += ' -o nt.raxml.tree' 
    cmd += ' -n nt.raxml.treefix.tree'
    cmd += ' -V 0'
    cmd += ' -l data/0/0.nt.raxml.treefix.log'
    cmd += ' -e " -m PROTGAMMAJTT"'
    cmd += ' tfix' + suffix + '/data/0/0.nt.raxml.tree'

    if not short:
        cmd += ' --niter=1000'
        cmd += ' --nquickiter=100'
        cmd += ' --freconroot=1'
    
    os.system(cmd)

    out = Tree('tfix' + suffix + '/data/0/0.nt.raxml.treefix.tree')
    os.system('rm -r tfix' + suffix)
    return out

def generateFakeSequence(domfile, l=100):
    """
    Generates a fake sequence with flanking regions of length l generated completely at random 
    """
    start = ''.join([np.random.choice(alphabet) for _ in range(l)])
    end = ''.join([np.random.choice(alphabet) for _ in range(l)])
    doms = list(open(domfile))[1::2]
    starts = [l]
    ends = [l + len(doms[0]) - 1]
    return starts, ends, start + np.random.choice(doms)[:-1] + end

###UTILS FOR RUNNING MULTREC

def map_isometric_nodes(template, iso_tree, lmap):
    """
    Given two isomorphic trees and a mapping of input -> template leaves, returns a full mapping
    of input -> template. WILL LIKELY BREAK IF TREES ARE NOT ISOMETRIC
    """
    fullmap = {}

    for leaf in iso_tree:
        inode = leaf
        tnode = lmap[leaf]
        while inode != None:
            fullmap[inode] = tnode
            inode = inode.up
            tnode = tnode.up
    
    return fullmap

def convert_trees(stree, gtree):
	s = Tree(stree, format=1)
	g = Tree(gtree, format=1)

	for leaf in g:
		leaf.name = "h" + leaf.name[1:]

	s = s.write(format=8)[:-1] + s.name + ';' #Stupid thing never writes the root's name for some reason
	g = g.write(format=9)

	return s, g

def run_multrec(stree, gtree):
    if platform == 'linux' or platform == 'linux2':
        cmd = "/home/caluru/Downloads/Multrec-master/Multrec/Multrec -d 2 -h 50 "
    else:
        cmd = "/Users/chaitanya/Downloads/Multrec-master/Multrec/Multrec -d 2 -h 50 "
    cmd += '-l 1 -spsep _ -o output.txt -g "' + gtree + '"' + ' -s "' + stree + '"'
    os.system(cmd)
    f = list(open('output.txt'))
    cost = int(f[1].strip())
    gindex = f.index("<GENETREES>\n") + 1
    gtree = f[gindex]
    return cost, gtree

def parse_multrec(multrec_string, host):
    mtree = Tree(multrec_string, format=1)
    levelorder = [node for node in mtree.traverse()][::-1]
    fullmap = {}

    for node in levelorder:
        hname = node.name.split("_")[0]
        if '-' in hname:
            hname = hname.split('-')[0]

        fullmap[node] = host&hname


    return mtree, fullmap

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
            print a in duplications
        return int(a.get_distance(b, topology_only=True)) + 1 if ga in duplications else 0

    loss_cost = 0

    for node in g.traverse():
        if node.children == []:
            continue
        up, c1, c2 = mapping[node], mapping[node.children[0]], mapping[node.children[1]]
        loss_cost += lcost * (get_losses(up, c1, node, duplications) + get_losses(up, c2, node, duplications))

    dup_cost = dupcost * len(tds)

    return loss_cost + dup_cost