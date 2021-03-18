# Tandem Duplications

This tool finds tandem duplications of domains using the method described in INSERT_REFERENCE_HERE. We take as input a gene tree and a species tree, as well as the relative ordering of domains on extant sequences and output a list of tandem duplications. Below are instructions on how to install dependencies, test the code, and run your own examples.

## Setup

This tool runs using Python 2.7+, and requires installations of the ete3 Tree package and the gurobi ILP solver. The ete3 Tree package can be downloaded using pip or easyinstall:

`pip install --upgrade ete3`  
or  
`easy_install -U ete3`

For full instructions, see [the download page](http://etetoolkit.org/download/)

We use Gurobi's python api as our ILP solver. In the future, we hope to add support for FOSS ILP solvers as well. However, Gurobi is free to use for academic purposes. See [this page](https://www.gurobi.com/academia/academic-program-and-licenses/) for full instructions. Briefly, the setup requires 3 steps:

1. Download Gurobi from the [downloads](https://www.gurobi.com/downloads/) page
2. Get a Gurobi academic license [here](https://www.gurobi.com/downloads/end-user-license-agreement-academic/)
3. Install your license using `grbgetkey <YOUR_KEY_HERE>`

## Testing

You can test the code using our test script on the simulated examples in the data folder. Our test script takes two optional arguments, type and edist. If run with no arguments, the default options are used.

* `type`: The type of solver to use. Choose between an `exact` ILP solution to the TDL Reconciliation Problem, or a fast `heuristic`. (Default: heuristic)
* `edist`: The average event distance between events in the simulated examples. The smaller the event distance, the larger the number of domain level events per domain tree. Examples are given for event distances of `0.1`, `0.075`, `0.05`, `0.025`, `0.01`. (Default: 0.1)

To run the tests, run 

`python test.py -type <exact/heuristic> -edist <0.1/0.075/0.05/0.025/0.01>`

When running with the exact solver, note that event distances less than 0.05 may take a very long time to run.

## Input Format

To infer tandem duplications and single losses in a domain tree, our tool requires as input a gene tree, a domain tree, a mapping file and a position file. Note that the tool can instead be used to infer events at the gene level by using a species tree and domain tree respectively. The input formats are as follows:

* **Gene Tree**: The input gene tree in newick format. Only names for the leaf genes are required.<br>  Example: `'((A,B),C);'`
* **Domain Tree**: The input domain tree in newick format. Only names for the leaf genes are required.<br>  Example: `'((A,B),C);'` 
