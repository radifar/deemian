# Quick Tour

For example using the following script as the input:

```text
molecule 5nzn.pdb [
	select protein_A = chain A and protein
	select oseltamivir = chain A and resname G39
	assign bond oseltamivir template CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)N)C(=O)O
]

measure protein_ligand [
	interactions all
	ionizable positive true
	ionizable negative true
	between oseltamivir and protein_A
]

present protein_ligand [
	interactions protein_ligand.txt
	deemiandata protein_ligand.dd
]
```

***

Will generate two output file, the Deemian data file (`protein_ligand.dd`) which can be opened with Deemian Viewer:

```{figure} ../_static/Deemian_viewer_n1_example.png
:align: center

**Figure 1.** Deemian Viewer shows the interactive 3D molecule visualization on the left and interaction data on the right.
```

***

And the following readable output (`protein_ligand.txt`):

```text
Deemian version: 0.1.0

interaction for "oseltamivir:protein_A":
conf          1
    ELECTROSTATIC as_cation:

                         oseltamivir protein_A

    id                   4           562     563     1053
    atom_name            N4          OE1     OE2     OD1
    res_name             G39         GLU     GLU     ASP
    res_num.chain        503.A       119.A   119.A   151.A
    distance                         4.206   2.744   3.157

    ELECTROSTATIC as_anion:

                         oseltamivir protein_A

    id                   0           539     540     541     4295    4296    4297
    atom_name            O1A         CZ      NH1     NH2     CZ      NH1     NH2
    res_name             G39         ARG     ARG     ARG     ARG     ARG     ARG
    res_num.chain        503.A       118.A   118.A   118.A   368.A   368.A   368.A
    distance                         3.646   2.800   3.584   3.625   2.832   3.523

    id                   1           3182    3183    3184    4295    4296    4297
    atom_name            O1B         CZ      NH1     NH2     CZ      NH1     NH2
    res_name             G39         ARG     ARG     ARG     ARG     ARG     ARG
    res_num.chain        503.A       293.A   293.A   293.A   368.A   368.A   368.A
    distance                         3.547   3.181   3.015   3.702   3.677   2.833
```
