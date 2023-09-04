# New Bacillus genus species

## Assembly
The new species were assembled using the Perfect Genome Assembly tutorial.

## Phylogeny
We used the `gtdb-tk` we build the phylogenetic tree of the assemblies.
More specificaly the `de_novo_wf` command.

The taxonomy was transformed with the following oneliner:

```
gawk -F"\t" '{taxon=$1; split($2,taxonomy, "; ") ; for (i in taxonomy){ split(taxonomy[i], classification,"__"); print taxon "\t" taxonomy[i] "\t" classification[1] "\t" classification[2]}}' gtdbtk.bac120.decorated.tree-taxonomy > gtdbtk.taxonomy_long.tsv
```

Then the results are loaded in R.

## Coding environment

All computations were performed in the Zorbas HPC facility of IMBBC-HCMR.


