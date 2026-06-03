

 GOParsimony 
=============

GOParsimony is a fast, robust bioinformatics package designed for parsimonious Gene Ontology (GO) 
functional annotation. 

It calculates hypergeometric enrichment statistics and maps the results to the hierarchical GO Directed Acyclic Graph, 
filtering them with a "parsimonious saltation" procedure across the DAG to assign cascades of "Precedence" scores. 

"Precedence" is simply a way to rank the terms according to their functional importance. This ranking enables the 
drawing of informative graphs. First-level precedence (precedence =1) is the parsimoniously least redundant list, 
but further precedence levels need to be added in to generate clusters indicating possible functional modules. 

The package uses a comprehensive Python visualization script to generate publication-ready figures, including coloured 
Precedence DAGs, Jaccard-distance clustered 2D heatmaps of shared GO terms that are useful for picking out possible functional 
modules, and comparative Venn diagrams and bar charts.

The package consists of three main components:

GOParsimony.c:         The core C program for parsing the OBO ontology and calculating statistics.
precedence-graphs.py:  The Python visualization script.
run-GOParsimony.sh:    A convenient shell script to execute the entire pipeline automatically.

It is licensed with a 3-clause BSD licence. 



INSTALL / COMPILE
================= 

To use GOParsimony, you must compile the C program and ensure your Python environment has the required scientific libraries installed.

Step 1: Compile the C Program
----------------------------- 
GOParsimony is written in standard C and requires the math library. 
Open your terminal and run the following command to compile it into an executable named "GOParsimony":

gcc -o GOParsimony GOParsimony.c -lm


Step 2: Install Python Dependencies
----------------------------------- 
The visualization script requires Python 3 and several scientific libraries. You can install all required packages 
using pip3 (or pip):

pip3 install pandas numpy scipy matplotlib seaborn graphviz matplotlib-venn


Step 3: Install the Graphviz System Binary
------------------------------------------ 
The Python graphviz library requires the underlying Graphviz engine to be installed on your operating system to render the DAG PDFs.

- On macOS: 
   brew install graphviz

   If you do not have Homebrew installed, first run the following command in your terminal:
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    and follow the instructions to add Homebrew to your PATH. 


- On Ubuntu/Debian linux: 
   sudo apt-get install graphviz

   If you do not have Homebrew installed, run: 
   sudo apt-get update && sudo apt-get install graphviz 



Step 4: Make the Shell Script Executable (Optional but recommended)
----------------------------------------
If you plan to use the automated bash script to run the pipeline, give it execution permissions:

chmod +x run-GOParsimony.sh
(Note: Ensure the variable C_EXEC="./GOParsimony" is set correctly inside the shell script).



USAGE
===== 

Inputs: 
-------

GOParsimony requires three primary inputs to run an enrichment analysis:

 (1) A Gene Ontology definition file (e.g., go.obo or go-basic.obo). 

     A Gene Ontology definition (OBO) file is provided as default (go-basic.obo), that was downloaded 
     from https://geneontology.org/docs/download-ontology/ . 

     It needs to be uncompressed as it is compressed to go-basic.obo.zip . 

     The user can specify one with the -o option. 

 (2) Background annotation file (a 2-column TSV/TXT mapping ID/accession to GO term, or a standard GAF file).

     Annotations in GAF format can be downloaded from https://geneontology.org/docs/download-go-annotations/ 
     and the links listed on that page. 

 (3) A Sample file (a TXT file listing IDs/accessions for proteins/genes (consistent with IDs/accessions in the 
     background annotation file).


Optional inputs include:

A GO Slims file (-m) to tag specific broadly-categorized terms. 

A GOSlims file 'goslim_generic.obo' is bundled with the package. GOSlims OBO files can be downloaded from 
https://geneontology.org/docs/download-ontology/ . 

A highlight list (-l) containing specific GO term IDs (one per line) to visually emphasize with green borders in the graphical output.



Flags used:
----------- 
-o: Path to the Gene Ontology definition (OBO format) file
-b: Path to the background annotations (in GAF or two-column {ID + GOterm} format) 
-s: Path to the sample list of IDs or accessions (consistent with the background annotations file) 
-p: Prefix for output files (e.g., "my_experiment")
-x: Maximum precedence threshold for the concise graphs and heatmaps (default is 5)
-m: (Optional) Path to GO Slims file
-l: (Optional) Path to Highlighted terms list
-e: (Optional) To keep the intermediate "Expanded" annotation text files on disk



Outputs:
-------- 

The pipeline produces three core data tables and a suite of PDF visualizations:

([prefix] refers to the option prefix that can be specified, see below)


[prefix]_Enrichment_Results.tsv:              The statistical results for all tested terms.

[prefix]_Significantly_Annotated_List.tsv:    A filtered cross-reference of significant proteins/genes and their terms.

[prefix]_Precedence_Edges.tsv:                The structural relationships between significant terms.

[prefix].Precedence-Graphs.pdf:               The complete, color-coded Precedence DAGs.

[prefix].Concise_Precedence_Graphs.pdf:       A filtered DAG showing only terms up to your specified max precedence, 
                                               using dotted lines to indicate skipped ancestors.

[prefix].GOTerm-table.pdf:                    A color-matched, sorted data table.

[prefix].2D-Heatmap.pdf:                      A heatmap dual-clustered using Jaccard distance to show shared GO terms 
                                               across significant proteins/genes, and to help infer possible functional 
                                               modules.

[prefix].Comparison-Graphs.pdf:               Analytical bar charts and Venn diagrams comparing added terms and sets 
                                               (generated if max precedence > 1 or if GOSlims or highlight files are provided).


Examples:
--------- 

The easiest way to run the analysis is by running the provided shell script, which handles passing files between 
the C program and the Python script automatically.

The user can try out these example runs using the example inputs in the example/ subdirectory. 

Please make sure to uncompress go-basic.obo.zip and example/yeast.background.gaf.zip . 

Example 1: 

./run-GOParsimony.sh -o go-basic.obo -b example/yeast.background.two-column.txt -s example/yeast-AQrich-cluster.sample.txt -p Example-1

This is a basic run using two-column input files for the background file of annotations 
and sample list of IDs/accessions (consistent with the background file). 


Example 2: 

./run-GOParsimony.sh -b example/yeast.background.gaf -s example/yeast-AQrich-cluster.sample.txt  -m goslim_generic.obo -p Example-2 -x 3 

(To run this example you need to uncompress example/yeast.background.gaf first.)
This example uses a GAF format file for background annotation input, the GOSlims file goslims.obo, 
and specifies a maximum precedence of 3 for graph drawing. 


Example 3: 

Running the Components Manually

If you prefer to run the steps individually, you can execute the C program and Python script sequentially, as in the 
following: 

Step A: Run the C parser
./GOParsimony  -b example/yeast.background.gaf -s example/yeast-AQrich-cluster.sample.txt -p Example-3 

(To run this example you need to uncompress example/yeast.background.gaf first.)

Step B: Run the Python visualization script

python3 precedence-graphs.py --results Example-3_Enrichment_Results.tsv  --edges Example-3_Precedence_Edges.tsv  --sig_prots Example-3_Significantly_Annotated_List.tsv  --out Example-3  --max_prec 3 





