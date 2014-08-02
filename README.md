This is a collection of scripts to map Pfam domains to targets from the ChEMBL database and to heuristically label all attainable targets with the protein domain that mediates small molecule binding. Mappings that are up-to-date with the current version of ChEMBL are provided on a dedicated [website](http://www.ebi.ac.uk/~fkrueger/mapChEMBLPfam/).

The script master.py specifies the workflow and with a local MySQL instance of ChEMBL, you can run through the entire set of scripts with the command:

    python master.py 

The parameters to this script are provided in a config file called mpf.yaml. You will have to generate this file, using example.yaml as a template.

The mapping of small molecule binding to structural domains is detailed in a [conference supplement](http://www.biomedcentral.com/bmcbioinformatics/supplements) in BMC Bioinformatics. If you are interested in the exact commit used to obtain the results - it is on the master branch and will be named "publication commit".
