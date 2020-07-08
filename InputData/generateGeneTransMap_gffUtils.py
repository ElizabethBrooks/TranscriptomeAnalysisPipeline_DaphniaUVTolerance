library(reticulate)

# check which python I am using
py_discover_config()
import sys
print(sys.version)
import gffutils
import itertools
import os
os.listdir()
db = gffutils.create_db("PA42.3.0.annotation.18440.gff", ":memory:", force = True,merge_strategy="merge", id_spec={'gene': 'Dbxref'})
list(db.featuretypes())
tx_and_gene=[]
## instead, loop over all features in the database
with open("tx2gene_NCBI.tsv", "w") as f:
        for feature in db.all_features():
                transcript = feature.attributes.get('transcript_id', [None])[0]
                gene = feature.attributes.get('gene', [None])[0]
                if gene and transcript and ([transcript, gene] not in tx_and_gene):
                        tx_and_gene.append([transcript, gene])
                        f.write(transcript + "\t" + gene + "\n")