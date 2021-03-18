#Given a starting tree, perform search around tree to find tree with high likelihood and low reconciliation score

from ete3 import Tree
import numpy as np
import os
import logging
from ilp_generator import *
from gurobipy import *
from TreeUtils import raxml, raxml_score, writeReconciliation

def spr(tree, subtree, new_sibling):
    """
    Performs a subtree prune and regraft operation moving a subtree to 
    a new location in the tree

    Arguments:
    tree (Tree): The full tree in which to perform an spr move
    subtree (Tree): The subtree that will be pruned and grafted
    new_sibling (Tree): The new sibling of the subtree being moved

    >>> t = Tree('((A,D),(B,C));')
    >>> subtree = t&"A"
    >>> new_sibling = t&"B"
    >>> t = spr(t, subtree, new_sibling)
    >>> newtree = Tree('(D,(C,(B,A)));')
    >>> rf = newtree.robinson_foulds(t)[0]
    >>> assert(rf == 0)
    """
    
    if tree == subtree:
        raise ValueError
    if subtree == new_sibling or subtree.up == new_sibling:
        raise ValueError
    if subtree.up == new_sibling.up:
        raise ValueError
    if tree.get_common_ancestor(subtree, new_sibling) == subtree:
        raise ValueError

    #CASE 1 (int -> int):
    if subtree.up != tree and new_sibling != tree:
        #Add node between new_sibling and its parent
        temp = Tree()
        temp.up = new_sibling.up
        temp.children = [subtree, new_sibling]
        new_sibling.up = temp
        if temp.up.children[0] == new_sibling:
            temp.up.children[0] = temp
        else:
            temp.up.children[1] = temp

        #Remove subtree from its current location
        old_parent = subtree.up
        subtree.up = temp
        temp.name = old_parent.name

        #Remove old parent
        ancestor = old_parent.up
        if old_parent.children[0] == subtree:
            other_child = old_parent.children[1]
        else:
            other_child = old_parent.children[0]

        other_child.up = ancestor
        if ancestor.children[0] == old_parent:
            ancestor.children[0] = other_child
        else:
            ancestor.children[1] = other_child

    #CASE 2 (cor -> int)
    elif subtree.up == tree:

        old_root = tree
        if tree.children[0] == subtree:
            tree = tree.children[1]
        else:
            tree = tree.children[0]
        tree.up = None

        old_root.up = new_sibling.up
        old_root.children = [subtree, new_sibling]
        new_sibling.up = old_root
        if old_root.up.children[0] == new_sibling:
            old_root.up.children[0] = old_root
        else:
            old_root.up.children[1] = old_root

    #CASE 3 (int -> root)
    else:

        temp = Tree()
        temp.up = None
        temp.children = [tree, subtree]
        tree.up = temp
        tree = temp
        temp.name = subtree.up.name

        old_parent = subtree.up
        subtree.up = tree

        #Remove old parent
        ancestor = old_parent.up
        if old_parent.children[0] == subtree:
            other_child = old_parent.children[1]
        else:
            other_child = old_parent.children[0]

        other_child.up = ancestor
        if ancestor.children[0] == old_parent:
            ancestor.children[0] = other_child
        else:
            ancestor.children[1] = other_child

    return tree

def pick_spr(tree):
    """
    ASSUMES INPUT TREE IS ROOTED
    Picks (and performs) a random spr move by picking a subtree to move 
    and a location to move it to
-
    Arguments:
    tree: The tree to perform a random spr for
    """
    nodes = [i for i in tree.traverse()][1:]

    remaining = set()
    while len(remaining) == 0:
        subtree = np.random.choice(nodes)
        invalid = set([i for i in subtree.traverse()])
        invalid.add(subtree.up)
        invalid.update(subtree.up.children)
        remaining = set([i for i in tree.traverse()]) - invalid

    ns = np.random.choice(list(remaining))
    tree = spr(tree, subtree, ns)
    return tree, (subtree.name, ns.name)

def pick_sprs(tree, n):
    """
    Tries to pick n unique sprs of a tree and return a list of those sprs.
    If that many are difficult to find, may give up early and return fewer

    Arguments:
    tree: The tree to find random sprs for
    n: The number of unique trees to attempt to find

    Output:
    the list of sprs
    """
    i = 0
    for node in tree.traverse():
        node.add_feature('label', str(i))
        i += 1

    sprs = []
    seen = set()
    i = 0
    failcount = 0

    while i < n:

        newTree = tree.copy()
        newTree, used = pick_spr(newTree)

        if used not in seen:

            seen.add(used)
            sprs.append(newTree)
            i += 1
            failcount = 0

        else:
            failcount += 1
            if failcount >= 100:
                return sprs
                
    logging.debug('Found ' + str(len(sprs)) + ' sprs out of ' + str(n))
    return sprs

def generate_rootings(tree):
    """
    Takes a rooted tree as input and generates all possible rerootings 
    of that tree. Does not change the input tree

    Arguments:
    tree (Tree): The rooted tree to produce alternate rootings of
    """
    trees = []
    copy = tree.copy()
    copy.unroot()

    #Every node needs a unique name for this to work
    i = np.random.randint(32768)
    for node in copy.traverse():
        if 'g' not in node.name:
            node.name = str(i)
            i += 1

    for node in copy.traverse():
        if node == copy:
            continue
        temp = copy.copy()
        temp.set_outgroup(temp&(node.name))
        trees.append(temp)

    return trees

#TODO: Test this function
def reconcile(host, guest, leafmap):
    """
    Performs TDL reconciliation on the guest tree and returns both the cost and the full
    mapping of guest -> host.

    Args:
    host (Tree): The host tree in ete3 format
    guest (Tree): The guest tree in ete3 format
    leafmapping (dict): guest -> host mapping containing all leaf nodes in the guest tree

    Output:
    cost (float): The cost of the min cost reconciliation with L = 1 and 
                  D(n) = 2 + .5(n-1)
    fullmap (dict): guest -> mapping containing all nodes in the guest tree
    """
    h = createTreeRepresentation(host)
    g = createTreeRepresentation(guest)
    d = createDistMatrix(host)
    e = createEMatrix(guest)
    mapping = createMapping(leafmap)

    eqnames, rhs, coldict = createEqns(host, guest, e, h, g, mapping, d)
    write('treesolve.mps', eqnames, rhs, coldict)

    m = read('treesolve.mps') # pylint: disable=undefined-variable
    m.Params.outputflag = 0
    m.optimize()
    cost = m.getObjective().getValue()

    fullmap = extractMapping(m, host, guest)

    os.system('rm treesolve.mps')
    return cost, fullmap

#TODO: Test this function
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

def reroot(host, guest, leafmap, recModule=reconcileDL):
    """
    Roots the input guest tree by checking all possible rootings for the lowest 
    reconciliation score. If the midpoint rooting has the lowest possible score,
    it will be used. #TODO: make the second sentence true.

    Args:
    host (Tree): The host tree to reconcile against
    guest (Tree): The guest tree to reroot by reconciliation
    leafmap (dict): guest -> host mapping including all guest leaves
    recModule (func): Function to reconcile guest with host. Must take as input
                      the host, guest, and leafmap arguments passed here and 
                      output (cost, fullmap) where cost is the float cost of the 
                      reconciliation and fullmap is a guest -> host mapping 
                      including all guest nodes.

    Output:
    The rooting of the guest tree with the lowest reconciliation cost.
    """
    trees = generate_rootings(guest)
    costs = []
    
    for tree in trees:
        #Generate new mapping
        newmap = {}
        for node in leafmap:
            newguest = tree&(node.name)
            newmap[newguest] = leafmap[node]
        costs.append(recModule(host, tree, newmap)[0])

    index = np.argmin(costs)
    return trees[index]

def perform_strict_search(sequences, host, guest, leafmap, num_iter=100):
    """
    Performs a search in tree space surrounding the highest scoring guest 
    tree according to raxml. Tries to find a guest with a similar likelihood 
    value and a better reconciliation score.

    Steps:
    1) Determine reconciliation score of guest w.r.t. host
    2) Generate 100 SPR moves, test all for reconciliation and raxml scores
    3) Pick a better tree, move to it
    4) Repeat 2-3 100 times
    """

    bestTree = guest
    bestScore = reconcileDL(host, guest, leafmap)[0]

    #Base nodemap on names, not actual nodes
    nodemap = {}
    for node in leafmap:
        nodemap[node.name] = leafmap[node].name

    for iteration in range(num_iter):
        logging.info('Iteration number ' + str(iteration))
        sprs = pick_sprs(guest, 200)
        scores = raxml_score(bestTree, sprs, sequences)[1]
        good_trees = [sprs[i] for i in range(len(sprs)) if scores[i] == 0]
        logging.debug('Found ' + str(len(good_trees)) + ' close trees')

        rec_scores = []
        for i in range(len(good_trees)):
            logging.debug('evaluated tree ' + str(i))
            tree = good_trees[i]
            lmap = {}
            for node in tree:
                lmap[node] = host&(nodemap[node.name])
            good_trees[i] = reroot(host, tree, lmap, reconcileDL) #Maybe use reconcile() instead? (Test performance)

            tree = good_trees[i]
            lmap = {}
            for node in tree:
                lmap[node] = host&(nodemap[node.name])
            rec_scores.append(reconcileDL(host, tree, lmap)[0]) #Replace this with reconcile() after testing

        """
        index = np.argmin(rec_scores)
        newScore = rec_scores[index]
        """

        #Sort into better, worse
        betterTrees, worseTrees = [], []
        betterScores, worseScores = [], []

        for i in range(len(rec_scores)):
            if rec_scores[i] >= bestScore:
                worseTrees.append(good_trees[i])
                worseScores.append(rec_scores[i])
            else:
                betterTrees.append(good_trees[i])
                betterScores.append(rec_scores[i])

        dice_roll = np.random.random()
        if dice_roll < 0.1 or len(betterTrees) == 0:
            index = np.random.randint(len(worseTrees))
            logging.info('Picking a worse tree, new: ' + str(worseScores[index]) + ' old: ' + str(bestScore))
            bestTree = worseTrees[index]
            bestScore = worseScores[index]
        elif dice_roll < 0.2:
            index = np.random.randint(len(betterTrees))
            logging.info('Picking a better tree, new: ' + str(betterScores[index]) + ' old: ' + str(bestScore))
            bestTree = betterTrees[index]
            bestScore = betterScores[index] 
        else:
            index = np.argmin(betterScores)
            logging.info('Picking best tree, new: ' + str(betterScores[index]) + ' old: ' + str(bestScore))
            bestTree = betterTrees[index]
            bestScore = betterScores[index]

        """
        if newScore > bestScore:
            logging.debug('Did not find a better tree')
        else:
            logging.info('Found better tree, new: ' + str(newScore) + ' old: ' + str(bestScore))
            bestTree = good_trees[index]
            bestScore = newScore
        """

    return bestTree

def perform_search(sequences, host, guest, leafmap, num_iter=100, len_search_path=100):

    bestTree = guest
    bestScore = reconcileDL(host, guest, leafmap)[0]
    seen = set() #Set of trees already seen by algo

    #Base nodemap on names, not actual nodes
    nodemap = {}
    for node in leafmap:
        nodemap[node.name] = leafmap[node].name

    for iteration in range(num_iter):

        #May change in inner loop, don't change bestScore or bestTree
        curBestScore = bestScore
        curBestTree = bestTree

        pool = [curBestTree]
        costs = [curBestScore]
        
        for i in range(len_search_path):

            #Propose new tree
            newTree, treeHash = pick_spr(curBestTree.copy())

            #Don't waste time on a tree that has already been seen
            if treeHash in seen:
                continue

            seen.add(treeHash)

            #Regenerate mapping
            lmap = {}
            for node in newTree:
                lmap[node] = host&(nodemap[node.name])

            #Reroot with low probability
            if np.random.random() < 0.05:
                newTree = reroot(host, newTree, lmap, reconcileDL)
                #Need to redo this because reroot returns a new tree
                lmap = {}
                for node in newTree:
                    lmap[node] = host&(nodemap[node.name])

            rec_score = reconcileDL(host, newTree, lmap)[0]

            pool.append(newTree)
            costs.append(rec_score)

            if rec_score < curBestScore \
                or (rec_score == curBestScore and np.random.random() < .5) \
                or (rec_score > curBestScore and np.random.random() < .1):
                curBestScore = rec_score
                curBestTree = newTree

        #Remove any tree that doesn't pass the RAxML threshold
        raxml_costs = raxml_score(guest, pool, sequences)[1]
        for i in xrange(len(raxml_costs)-1, -1, -1):
            if raxml_costs[i] != 0:
                pool.pop(i)
                costs.pop(i)

        #pick the best tree for the next iteration
        if len(costs) != 0:
            index = np.argmin(costs)
            bestScore = costs[index]
            bestTree = pool[index]

        logstring = 'Iteration Number ' + str(iteration)  + ": "
        if index != 0:
            logstring += "Better tree found, " + str(bestScore)
        else:
            logstring += "Better tree not found"
        logging.info(logstring)

    return bestScore, bestTree

def best_of_n(sequences, host, guest, leafmap, num_iter=100, len_search_path=100, n=5):
    bestTree = None
    bestScore = float('inf')

    for _ in range(n):
        h = host.copy()
        g = guest.copy()
        bs, bt = perform_search(sequences, h, g, genMap(h, g), num_iter, len_search_path)
        if bs < bestScore:
            bestScore = bs
            bestTree = bt

    return bestScore, bestTree

def genMap(host, guest):
    #{guest -> host}
    nodemap = {}
    for leaf in guest:
        hname = 'h' + leaf.name.split("_")[0][1:]
        nodemap[leaf] = host&hname
    return nodemap

def name(t):
    i=0
    for node in t.traverse():
        if 'g' not in node.name:
            node.name = str(i)
            i += 1

