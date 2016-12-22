Other Plotting
**************

plot_kmer
---------

Plot signal distribution acorss kmers. The command allows the selection of the number of bases up and downstream around the central genomic base.

.. figure::  _images/kmer_dist.jpg
   :align:   center
   :scale: 50%
   
   K-mer signal distribution plotting example.

.. figure::  _images/kmer_dist_read.jpg
   :align:   center
   :scale: 50%
   
   K-mer signal distribution connecting each read example.


cluster_most_significant
------------------------

Cluster raw signal traces at bases with significant differences. Clustering is done on the difference between two runs (generally a native and amplified sample) centered on the most significantly different bases. Note that this command can be quite slow for many points, so it is advised to only use few points for plotting (<1000).

.. figure::  _images/cluster.jpg
   :align:   center
   :scale: 50%
   
   Signal clustering plotting example.
