# Table

`ganon table` filters and summarizes several reports obtained with `ganon report` into a table. Filters for each sample or for averages among all samples can also be applied.

## Examples

Given several `.tre` from `ganon report`:

### Counts of species

```bash
ganon table --input *.tre --output-file table.tsv --rank species
```

### Abundance of species

```bash
ganon table --input *.tre --output-file table.tsv --output-value percentage --rank species
```

### Top 10 species (among all samples)

```bash
ganon table --input *.tre --output-file table.tsv --output-value percentage --rank species --top-all 10
```

### Top 10 species (from each samples)

```bash
ganon table --input *.tre --output-file table.tsv --output-value percentage --rank species --top-sample 10
```

### Filtering results 

```bash
ganon table --input *.tre --output-file table.tsv --output-value percentage --rank species --min-count 0.0005
```

This will keep only results with a min. abundance of `0.05%`.
