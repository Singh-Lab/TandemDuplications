# -*- coding: utf-8 -*-
#Useful functions for internal testing (unrelated to tree building)

import sys
import os

#From: https://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console 
# Print iterations progress
def printProgressBar(iteration, total, prefix='Progress', suffix='', decimals=1, bar_length=60):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        bar_length  - Optional  : character length of bar (Int)
    """
    str_format = "{0:." + str(decimals) + "f}"
    percents = str_format.format(100 * (iteration / float(total)))
    filled_length = int(round(bar_length * iteration / float(total)))
    bar = '█' * filled_length + '-' * (bar_length - filled_length)

    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),

    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()

#From: https://stackoverflow.com/questions/6618515/sorting-list-based-on-values-from-another-list
def sortBy(a,b):
    """Sorts list a by the values in list b. Not in place."""
    return [x for (_,x) in sorted(zip(b,a))]

def msa(infile, outfile):
    """
    os.system('cp orthogroup.fa genes.fa')
    """
    if len(list(open(infile, 'r'))) == 2:
                if infile != outfile:
                        os.system('cp ' + infile + ' ' + outfile)
    else:
        os.system('clustalo --threads=8 -i ' + infile + ' -o ' + outfile + ' -t Protein --force')