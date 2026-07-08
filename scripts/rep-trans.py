import sys
from multitax import GtdbTx, NcbiTx, CustomTx
from ganon.report import parse_rep

rep_file = sys.argv[1]
gtdb_tax_file = sys.argv[2]
ncbi_tax_file = sys.argv[3]
out = sys.argv[4]

gtdb_tax = GtdbTx.from_customtx(CustomTx(version="232", files=gtdb_tax_file, cols=["node", "parent", "rank", "name"]))
ncbi_tax = NcbiTx.from_customtx(CustomTx(version="x", files=ncbi_tax_file, cols=["node", "parent", "rank", "name"]))
#ncbi_tax = NcbiTx(files=ncbi_tax_file)
#ncbi_tax.filter(["2", "2157"], desc=True)
#ncbi_tax.write("ncbi_arc_bac.tax", cols=["node", "parent", "rank", "name"])
ncbi_tax.build_lca()
gtdb_tax.build_translation(ncbi_tax, representatives=False, file="/home/pirov/code/multitax/data/gtdb/232_acc_rep_lin_ncbi.tsv.gz")
reports, counts = parse_rep(rep_file, normalize=False)

new_rep = {}

for taxid, val in reports["H1"].items():
    translated_tax = []
    
    for t in gtdb_tax.lineage(taxid)[::-1]:
        translated_tax = gtdb_tax.translate(t)
        if translated_tax:
            break
    
    if len(translated_tax)==1:
        new_taxid = translated_tax.pop()
    elif len(translated_tax)>1:
        new_taxid = ncbi_tax.lca(translated_tax)
    else:
        print(taxid, "no translation find, assigning to root")
        new_taxid = ncbi_tax.root_node
    
    if new_taxid not in new_rep:
        new_rep[new_taxid] = {"direct_matches": 0, "unique_reads": 0 , "lca_reads": 0, "rank": ncbi_tax.rank(new_taxid), "name": ncbi_tax.name(new_taxid)}
    
    new_rep[new_taxid]["direct_matches"] += val["direct_matches"]
    new_rep[new_taxid]["unique_reads"] += val["unique_reads"]
    new_rep[new_taxid]["lca_reads"] += val["lca_reads"]

with open(out, "w") as file:
    for taxid, val in new_rep.items():
        print("H1", taxid, val["direct_matches"], val["unique_reads"], val["lca_reads"], val["rank"], val["name"], sep="\t", file=file)
    print("#total_classified", counts["total"]["reads"], sep="\t", file=file)
    print("#total_unclassified", counts["total"]["unclassified"], sep="\t", file=file)
