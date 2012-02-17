"""
  Function:  getLigandsForTarget
  --------------------

  momo.sander@ebi.ac.uk
"""  
def getLigandsForActivity(activity, release, user, pword, host, port): 


  import queryDevice

  
  ligands = queryDevice.queryDevice("SELECT DISTINCT act.molregno, standard_value,\
	 standard_type, standard_units, canonical_smiles, relation, act.activity_id \
										\
                                        FROM activities act \
                                        JOIN compound_structures cs\
                                            ON act.molregno = cs.molregno \
     				 WHERE activity_id = %s" \
                                        %activity, release, user, pword, host, port)
  return ligands

                    





def getLigandsForTarget(target, release, user, pword, host, port): 


  import queryDevice

  
  ligands = queryDevice.queryDevice("SELECT DISTINCT act.molregno, standard_value,\
	 standard_type, standard_units, canonical_smiles, act.relation, act.activity_id \
										\
                                        FROM activities act \
                                        JOIN assay2target a2t \
                                            ON act.assay_id = a2t.assay_id\
                                        JOIN target_dictionary td \
                                            ON a2t.tid = td.tid \
                                        JOIN assays ass \
                                            ON ass.assay_id = act.assay_id \
					JOIN compound_structures cs \
					    ON cs.molregno=act.molregno \
									\
     				 WHERE td.protein_accession = '%s' \
        			 AND ass.assay_type='B'  \
        			 AND act.relation ='='    \
        			 AND a2t.multi=0  \
        			 AND a2t.complex=0 \
        			 AND a2t.relationship_type = 'D'\
         			 AND act.standard_type IN('Ki','Kd','IC50', \
				'EC50','-Log Ki','pKd' , 'pA2', 'pI', 'pKa')" \
                                        %target, release, user, pword, host, port)
  return ligands

