import gffutils

def parseCG(protein):
    """ Pull out the CG number which is burried inside of the protein
        attributes. """
    attrib = protein.attributes['Dbxref']

    # If there are multiple different attributes this will be a list.
    if type(attrib) == list:
        # CG number is part of a string that starts with 'FlyBase_Annotation_IDs' 
        index = [i for i, s in enumerate(attrib) if 'FlyBase_Annotation_IDs' in s]
        for i in index:
            cgNumberIso = attrib[i].split(':')[1]

    # If there is only a single attribute then this will be a string.
    elif type(attrib) == str:
        cgNumberIso = attrib.split(':')[1]

    # return CG number
    return(cgNumberIso)

def getFbgn(Fbtr):
    """ Identify the FBgn number by using the FBtr number """
    parent = next(db.parents(Fbtr,featuretype='gene'))
    return(parent.id,parent.attributes['Name'])




if __name__ == '__main__':
    # Connect to GFF database
    db = gffutils.FeatureDB('/home/jfear/tmp/dsrp/dmel-all-r5.48.gff.db')

    # Since the feature types I am trying to get to are protein isoform names
    # (CG numbers) I will start by pulling out all of the proteins.
    proteins = db.features_of_type('protein')

    with open('/home/jfear/tmp/dsrp/symbol2proteinCG_test.csv','w') as OUT:
        OUT.write('symbol,FBgn,FBtr,protein_CG_number\n')
        
        for protein in proteins:
            # Skip mitochondrial proteins
            if protein.chrom != 'dmel_mitochondrion_genome':
                try:
                    proteinCgNumber = parseCG(protein)
                    Fbtr = protein.attributes['Derives_from']
                    Fbgn, symbol = getFbgn(Fbtr)
                    myout = ','.join([symbol, Fbgn, Fbtr, proteinCgNumber]) + "\n"
                    OUT.write(myout)
                except:
                    print 'exception: Something funky is going on! ', attrib
