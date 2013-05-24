This is a collection of scripts to map Pfam domains to targets from the ChEMBL database and to heuristically label all attainable targets with the protein domain that mediates small molecule binding. Mappings that are up-to-date with the current version of ChEMBL are provided on a dedicated [website](http://www.ebi.ac.uk/~fkrueger/mapChEMBLPfam/).

The script master.py specifies the workflow and with a local MySQL instance of ChEMBL, you can run through the entire set of scripts with the command:

    python master.py release user password host port

(Substituting release with the ChEMBL release number you are working with - the remaining variables should be replaced by your MySQL user credentials.)

The mapping of small molecule binding to structural domains is detailed in a [conference supplement](http://www.biomedcentral.com/bmcbioinformatics/supplements) in BMC Bioinformatics. If you are interested in the exact commit used to obtain the results - it is on the master branch and will be named "publication commit".

This commit reproduces the data and figures presented in the paper.
