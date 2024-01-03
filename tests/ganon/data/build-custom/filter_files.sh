# copy desired files to files/ and files/more/

# Get full files
mkdir base_files/
wget --quiet --output-document - "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt" | tail -n+3 > base_files/assembly_summary.txt
wget --quiet --output-document - "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq_historical.txt" | tail -n+3 >> base_files/assembly_summary.txt
wget --quiet --output-document - "ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt" | tail -n+3 >> base_files/assembly_summary.txt
wget --quiet --output-document - "ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank_historical.txt" | tail -n+3 >> base_files/assembly_summary.txt
wget "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz" --output-document base_files/new_taxdump.tar.gz
wget "https://data.gtdb.ecogenomic.org/releases/latest/bac120_taxonomy.tsv.gz" --output-document base_files/bac120_taxonomy.tsv.gz
wget "https://data.gtdb.ecogenomic.org/releases/latest/ar53_taxonomy.tsv.gz" --output-document base_files/ar53_taxonomy.tsv.gz
wget "https://data.gtdb.ecogenomic.org/releases/latest/bac120_metadata.tsv.gz" --output-document base_files/bac120_metadata.tsv.gz
wget "https://data.gtdb.ecogenomic.org/releases/latest/ar53_metadata.tsv.gz" --output-document base_files/ar53_metadata.tsv.gz
tar xf base_files/bac120_metadata.tar.gz bac120_metadata_r207.tsv
tar xf base_files/ar53_metadata.tar.gz ar53_metadata_r207.tsv
tar xf base_files/new_taxdump.tar.gz taxidlineage.dmp nodes.dmp names.dmp
mv bac120_metadata_r207.tsv ar53_metadata_r207.tsv taxidlineage.dmp nodes.dmp names.dmp base_files/

# extract accessions
zgrep -h -o '^>[^ ]*' files/*.fna.gz | sed 's/>//' > seq_acc.txt
zgrep -h -o '^>[^ ]*' files/more/*.fna.gz | sed 's/>//' >> seq_acc.txt
ls -1 files/*.fna.gz | grep -o "GC[FA]_[0-9.]*" > ass_acc.txt
ls -1 files/more/*.fna.gz | grep -o "GC[FA]_[0-9.]*" >> ass_acc.txt

# filter assembly_summary.txt 
join <(sort ass_acc.txt) <(sort -t$'\t' -k 1,1 base_files/assembly_summary.txt) -t$'\t' | sort | uniq > assembly_summary.txt

# filter taxdump.tar.gz
# get all taxids and lineage
cut -f 6 assembly_summary.txt | sort | uniq > taxids.txt
cut -f 6 assembly_summary.txt | xargs -I {} grep "^{}[^0-9]" base_files/taxidlineage.dmp | cut -f 3 | tr ' ' '\n' | sort | uniq | sed -r '/^\s*$/d' >> taxids.txt 
# add root and nodes
head -n 1 base_files/nodes.dmp > nodes.dmp
cat taxids.txt | xargs -I {} grep "^{}[^0-9]" base_files/nodes.dmp >> nodes.dmp
head -n 2 base_files/names.dmp > names.dmp
cat taxids.txt | xargs -I {} grep "^{}[^0-9]" base_files/names.dmp >>  names.dmp
# make new and filtered taxdump.tar.gz
touch merged.dmp
tar -czf taxdump.tar.gz nodes.dmp names.dmp merged.dmp
rm nodes.dmp names.dmp merged.dmp taxids.txt
# create simulated species_genome_size.txt.gz
echo "#species_taxid	min_ungapped_length	max_ungapped_length	expected_ungapped_length	number_of_genomes	method_determined" > species_genome_size.txt
cut -f 7 assembly_summary.txt | awk -v base=985251 'BEGIN{OFS="\t"}{print $1, base*NR, base*NR, base*NR, 1, "automatic"}' >> species_genome_size.txt
gzip species_genome_size.txt


# filter gtdb tax
cat ass_acc.txt | xargs -I {} zgrep "^[GR][BS]_{}" base_files/bac120_taxonomy.tsv.gz | gzip > bac120_taxonomy.tsv.gz
cat ass_acc.txt | xargs -I {} zgrep "^[GR][BS]_{}" base_files/ar53_taxonomy.tsv.gz | gzip > ar53_taxonomy.tsv.gz
# make metadata
head -n 1 base_files/bac120_metadata_r207.tsv > bac120_metadata_r207.tsv
head -n 1 base_files/ar53_metadata_r207.tsv > ar53_metadata_r207.tsv
join <(cut -f 1 assembly_summary.txt | sort) <(awk 'BEGIN{FS=OFS="\t"}{print substr($1,4,length($1)), $0}' base_files/bac120_metadata_r207.tsv | sort -t$'\t' -k 1,1) -t$'\t' | cut -f 2- >> bac120_metadata_r207.tsv
join <(cut -f 1 assembly_summary.txt | sort) <(awk 'BEGIN{FS=OFS="\t"}{print substr($1,4,length($1)), $0}' base_files/ar53_metadata_r207.tsv | sort -t$'\t' -k 1,1) -t$'\t' | cut -f 2- >> ar53_metadata_r207.tsv
gzip -c bac120_metadata_r207.tsv > bac120_metadata.tsv.gz
gzip -c ar53_metadata_r207.tsv > ar53_metadata.tsv.gz
rm bac120_metadata_r207.tsv ar53_metadata_r207.tsv

# make nucl_gb.accession2taxid.gz
../../../../scripts/ganon-get-seq-info.sh -i seq_acc.txt -e | awk 'BEGIN{FS=OFS="\t"; print "accession", "accession.version", "taxid", "gi"}{split($1,acc,"."); print acc[1], $1, $3, "0" }' | gzip > nucl_gb.accession2taxid.gz

rm seq_acc.txt ass_acc.txt
#rm -rf base_files
