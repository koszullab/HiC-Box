HiC-Box
=======

**HiC-Box** is a HiC data processing pipeline and visualizer, written mostly in Python (as of 2018, it is only compatible to Python 2; with deprecation coming up in 2020, work is undergoing to convert it to Python 3). It uses bowtie2 as a backend to align input paired-end reads onto an input genome, and derives a contact map from the alignment data and the position of restriction sites along the genome (the restriction enzyme is also given as input).

The restriction fragments are then binned (and the corresponding matrix sum-pooled) at different scales, hence building a *pyramid*. Each *level* of the pyramid is a contact map at a different resolution. Each such level can be browsed with the box.

HiC-Box generates datasets that are compatible with [GRAAL](http://github.com/koszullab/GRAAL) for reassembly. Both softwares operate on the same data template and their codebase is redundant but since GRAAL has specific requirements that can be difficult to deploy, they are kept on separate repos. Nevertheless, if you wish to use GRAAL on your own genome and paired-end reads, you will need to use the box to convert them into GRAAL-digestible input data. [Here](https://github.com/koszullab/GRAAL#datasets) are examples of what such prepared datasets should look like.

Dependencies
============

Python packages
---------------

* numpy
* scipy
* matplotlib
* h5py
* Biopython
* pymongo
* [mirnylib](https://bitbucket.org/mirnylab/mirnylib)
* cython (mirnylib dependency)
* wxpython
* networkx

Other software
--------------

* bowtie2

How to use
==========

* Make sure you have read/write/execute permission on the main folder and all its subfolders.
* Run main.py
* Choose output folder
* Load a fasta genome
* Make sure a bowtie index is built and placed in a directory named "index" within the main folder; this can be done automatically by clicking on *index* (bowtie path must be specified in the advanced options or a bowtie folder must be in the main folder)
* Load fastq paired-end read files
* Load restriction enzyme (or manually type in the site sequence)
* Check advanced parameters and tweak if needed
* Align and wait for *ready for computation* terminal output (this may take a while)
* Click on **Pyramid** to proceed to the visualizer (visualizer may also be directly summoned by running main_window.py)
* Load data and browse to the inside of a folder named *analysis* within the output folder
* Build pyramid (this may take a while) and load it
* Visualize using options as needed

Advanced
========

Visualizer options
------------------

* When building a pyramid, the *sparsity filter* will remove bins whose total contact count is too low. 
* Saturation threshold can be either absolute (fixed number), relative (percentile, append % at the end) or auto-adjusted depending on the needed rendering.
* Matrices can be manually loaded "raw" (from a text file) and binned, but no genome information will be displayed.
* "Kb binning" will regroup bins such that the total length of each bin is closest to 10kb.
* The norm called *sparsity* is an alias for fragment-wise order 1 normalization.

Output templating
=================

The box is designed to go through all steps of a typical Hi-C pipeline, from alignment to visualization or a possible reassembly by [GRAAL](http://github.com/koszullab/GRAAL). However, if you wish to customize the alignment step, you may directly provide a sam file named **0.sam** in the output folder. HiC-Box will skip the alignment part and proceed to the contact matrix generation.

Troubleshooting
===============

* Currently the PCR duplicate/tag detection step is very unoptimized and may crash the box if there are too many reads. For this reason (and because filtering/trimming reads on the fly every time can be clunky anyway) it is preferable that you trim and filter your reads prior to launching the pipeline. An example of script that does this can be found [here](https://github.com/koszullab/DADE/blob/master/pcr_duplicate_Hiseq.pl). Once your reads are prepared, be sure to set *tag length* to 0 in the options.




