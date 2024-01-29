# Choosing the right presentation format

By default Deemian uses `detailed_conf_first` for the readable output format.
There are six available format that can be used:

| Priority | Detailed | Clustered | Summarized |
|---|---|---|---|
| conformation first     | `detailed_conf_first` | `clustered_conf_first` | `summarized_conf_first` |
| interaction type first | `detailed_type_first` | `clustered_type_first` | `summarized_type_first` |

By now you should be familiar with the `detailed_conf_first` format as it is the one used in getting started and tutorials, now is time to find out the clustered and summarized one.
Clustered means that instead of listing the interaction one by one, the interaction will be clustered based on residue id and chain.
While the summarized one means that rather than showing the detail of interaction it will display the number of interaction per residue pair.

To quickly compare the three of them lets use the script from [getting started](../gettingstarted) and add a few lines to get three output format at once:


:::{card}
:class-header: sd-px-1 sd-py-1

n1-oseltamivir-5nzn-three-format.txt
^^^
```{code-block} text
:emphasize-lines: 15, 16, 17

molecule 5nzn.pdb [
    select protein_A = chain A and protein
    select oseltamivir = chain A and resname G39
    assign bond oseltamivir template CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)N)C(=O)O
]

measure protein_ligand [
    ionizable positive true
    ionizable negative true
    interactions all
    between oseltamivir and protein_A
]

present protein_ligand [
    interactions detailed_conf_first n1_g39_5nzn_detailed.txt
    interactions clustered_conf_first n1_g39_5nzn_clustered.txt
    interactions summarized_conf_first n1_g39_5nzn_summarized.txt
    deemiandata n1_g39_5nzn_ionizable_corrected.dd
]
```
:::


Running the script above will generate three output in readable formats:


::::::{tab-set}

:::::{tab-item} n1_g39_5nzn_detailed.txt
:class-label: sd-shadow-sm

```text
Deemian version: 0.0.0.post71+bad1636

interaction for "oseltamivir:protein_A":
conf          1
    ELECTROSTATIC as_cation:

                         oseltamivir protein_A

    id                   19          281     282     529
    atom_name            N4          OE1     OE2     OD1
    res_name             G39         GLU     GLU     ASP
    res_num.chain        503.A       119.A   119.A   151.A
    distance                         4.206   2.744   3.157

    ELECTROSTATIC as_anion:

                         oseltamivir protein_A

    id                   1           271     272     273     2202    2203    2204
    atom_name            O1A         CZ      NH1     NH2     CZ      NH1     NH2
    res_name             G39         ARG     ARG     ARG     ARG     ARG     ARG
    res_num.chain        503.A       118.A   118.A   118.A   368.A   368.A   368.A
    distance                         3.646   2.800   3.584   3.625   2.832   3.523

    id                   2           1624    1625    1626    2202    2203    2204
    atom_name            O1B         CZ      NH1     NH2     CZ      NH1     NH2
    res_name             G39         ARG     ARG     ARG     ARG     ARG     ARG
    res_num.chain        503.A       293.A   293.A   293.A   368.A   368.A   368.A
    distance                         3.547   3.181   3.015   3.702   3.677   2.833
```

:::::

:::::{tab-item} n1_g39_5nzn_clustered.txt
:class-label: sd-shadow-sm

```text
Deemian version: 0.0.0.post71+bad1636

interaction for "oseltamivir:protein_A":
conf          1
    ELECTROSTATIC as_cation:

                         oseltamivir protein_A

                         G39         GLU             ASP
                         503.A       119.A           151.A

    id                   19          281     282     529
    atom_name            N4          OE1     OE2     OD1
    distance                         4.206   2.744   3.157

    ELECTROSTATIC as_anion:

                         oseltamivir protein_A

                         G39         ARG                     ARG
                         503.A       118.A                   368.A

    id                   1           271     272     273     2202    2203    2204
    atom_name            O1A         CZ      NH1     NH2     CZ      NH1     NH2
    distance                         3.646   2.800   3.584   3.625   2.832   3.523

                         G39         ARG                     ARG
                         503.A       293.A                   368.A

    id                   2           1624    1625    1626    2202    2203    2204
    atom_name            O1B         CZ      NH1     NH2     CZ      NH1     NH2
    distance                         3.547   3.181   3.015   3.702   3.677   2.833
```

:::::

:::::{tab-item} n1_g39_5nzn_summarized.txt
:class-label: sd-shadow-sm

```text

Deemian version: 0.0.0.post71+bad1636
interaction for "oseltamivir:protein_A":
conf          1

    ELECTROSTATIC as_cation  has interactions with   interaction numbers
    oseltamivir              protein_A

    G39 503.A                GLU 503.A                       2
                             ASP 503.A                       1

    ELECTROSTATIC as_anion   has interactions with   interaction numbers
    oseltamivir              protein_A

    G39 503.A                ARG 503.A                       3
                             ARG 503.A                       3
                             ARG 503.A                       6
```

:::::

::::::


As shown above, we can see that in detailed mode and clustered mode the results are pretty similar.
The only difference is how the interactions are getting clustered based on the residue id in clustered mode.
The summarized mode give much simpler output as it leave the atom id and atom name details and only present the number of interactions.

The other format choice is *conformation first* versus *interaction type first*.
This option is only useful when there are multiple conformation being analyzed.
Sometime we would like to compare the interaction between conformations for specific interaction type (such as electrostatic), in this case using `detailed_type_first` or `clustered_type_first` would be more effective.
