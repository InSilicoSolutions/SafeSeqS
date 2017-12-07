SafeSeqS – Safe Sequencing System
---------------------------------------------

SafeSeqS Python pipeline for analysis of SafeSeqS DNA sequencing data.   https://www.ncbi.nlm.nih.gov/pubmed/21586637


The package supports two input types:    


SafeSeqS produces the following outputs:

The identification of mutations that are present in a small fraction of DNA templates.

Documentation
-------------

Please see the `documentation <https://github.com/xxxx/SafeSeqS/wiki>`_ on github for complete usage details.

Installation
------------

Use pip to install the SafeSeqS package:

.. code-block:: bash

    $ pip install SafeSeqS

SafeSeqS depends on several other packages and these will be installed automatically by pip if they
are not already installed.  

**Required packages:**

* scipy
* pypiwin32 (for Windows machines)
* resource (for Linux machines)

Note the pip scipy installation can sometimes have issues.  If so, alternate methods of installing scipy may be preferable.  See scipy `install <http://www.scipy.org/install.html>`_ for details.  Then rerun pip install SafeSeqS.
