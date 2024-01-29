# Analyze the human ACE2 structural recognition by SARS-Cov-2 Omicron variant

The COVID-19 pandemic was caused by a virus known as SARS-CoV-2.
It emerged in late December 2019 and, in just two years, over 260 million confirmed cases were reported and more than five million people lost their lives due to the virus.
During this time, many variants of the virus appeared, including the Omicron variant, which has around 50 mutations relative to the original Wuhan-Hu-1 variant.
More than 30 of these mutations were found in the spike protein (see Figure 1), which is a key factor in the virus's transmission and immune evasion.

:::{card}
:img-top: ../../_static/diagnostics-13-00559-g001.png
:class-card: sd-px-5 sd-mx-5
:class-img-top: sd-width-100

**Figure 1**. Overview of the mutations of SARS-CoV-2 Alpha (B.1.1.7) and Omicron (B.1.1.529) variants.

+++

Source: SARS-CoV-2 Omicron (B.1.1.529) Variant: A Challenge with COVID-19 (Afshar et.al., 2023) [link](https://doi.org/10.3390/diagnostics13030559)
:::


In this tutorial, we will focus on identifying some of the significant electrostatic interactions that occur between the human ACE-2 protein and the spike protein of the SARS-CoV-2 Omicron variant.
We will demonstrate how easy it is to create multiple selections and analyze the interactions between them.
Such segmentation of the interaction helps researchers communicate their findings more effectively.

We will also explore different output format.
By default the readable output is in `detailed format`, and in this tutorial we will try different output formats such as `clustered format` and `summarized format`.

## Constructing the script

First, lets get the SARS-CoV-2 Omicron receptor binding domain complexed with human ACE2 deposited on RCSB PDB under the ID [7u0n](https://www.rcsb.org/structure/7u0n).
As usual, we will write the script in three step.
In the first step, we will open `7u0n.pdb` and create three selection, two for each of the ACE2 binding hotspot and one for receptor binding motif (RBM) of Omicron spike protein.


:::{card}
:class-header: sd-px-1 sd-py-1

omicron-rbd-ace2-7u0n-multiselect.txt
^^^
```text
molecule 7u0n.pdb [
    select ace2_hotspot31 = chain A and protein and resid 1 to 55
    select ace2_hotspot353 = chain A and protein and resid 341 to 364
    select spike_rbm = chain E and protein and resid 470 to 510
]
```
:::


Here, `chain A` refers to the ACE2 protein, while `chain E` refers to the RBD of spike protein.
These informations are available both on RCSB PDB page and `7u0n.pdb` file.
To be sure, you can also do a quick check via molecule visualization tool of your choice (PyMOL, VMD, or else).
The `resid` keyword here means the residue id, it followed by `start_id to end_id` values, where `start_id` is the selection starting index and `end_id` is the selection end index.
The choice of residue id is selected based on the smallest unit of secondary structure that involved in the receptor binding process.

Now that we have create three selections, the next step is to setup the interaction analysis in `measure` step.
Here, we interested in evaluating four interaction pairs:

1. Interaction between spike protein and Hotspot 31 of ACE2 protein.
2. Interaction between spike protein and Hotspot 353 of ACE2 protein.
3. Internal interaction of ACE2 protein, specifically between Hotspot 31 and 353.
4. Internal interaction of RBM in spike protein.

Which translate to these instructions:


:::{card}
:class-header: sd-px-1 sd-py-1

omicron-rbd-ace2-7u0n-multiselect.txt
^^^
```{code-block} text
:emphasize-lines: 11-14

molecule 7u0n.pdb [
    select ace2_hotspot31 = chain A and protein and resid 1 to 55
    select ace2_hotspot353 = chain A and protein and resid 341 to 364
    select spike_rbm = chain E and protein and resid 470 to 510
]

measure ace2_spike_rbd [
    interactions all
    ionizable positive true
    ionizable negative true
    between spike_rbm and ace2_hotspot31 as rbm_h31
    between spike_rbm and ace2_hotspot353 as rbm_h353
    between ace2_hotspot31 and ace2_hotspot353 as internal_ace2
    between spike_rbm and spike_rbm as internal_rbm
]
```
:::


The final step is easy, we just have to declare the output file names.


:::{card}
:class-header: sd-px-1 sd-py-1

omicron-rbd-ace2-7u0n-multiselect.txt
^^^
```{code-block} text
:emphasize-lines: 18, 19

molecule 7u0n.pdb [
    select ace2_hotspot31 = chain A and protein and resid 1 to 55
    select ace2_hotspot353 = chain A and protein and resid 341 to 364
    select spike_rbm = chain E and protein and resid 470 to 510
]

measure ace2_spike_rbd [
    interactions all
    ionizable positive true
    ionizable negative true
    between spike_rbm and ace2_hotspot31 as rbm_h31
    between spike_rbm and ace2_hotspot353 as rbm_h353
    between ace2_hotspot31 and ace2_hotspot353 as internal_ace2
    between spike_rbm and spike_rbm as internal_rbm
]

present ace2_spike_rbd [
    interactions ace2_spike_rbd_detailed.txt
    deemiandata ace2_spike_rbd_detailed.dd
]
```
:::


## Running the script and analyzing the results


Given the script above, if we run it with Deemian we will get the following results:


::::::{tab-set}

:::::{tab-item} ace2_spike_rbd_detailed.txt
:class-label: sd-shadow-sm

```text
Deemian version: 0.0.0.post71+bad1636

interaction for "rbm_h31":
conf          1
    ELECTROSTATIC as_cation:

                         spike_rbm   ace2_hotspot31

    id                   10999       141
    atom_name            NH1         OE1
    res_name             ARG         GLU
    res_num.chain        493.E       35.A
    distance                         4.363

    id                   11043       163
    atom_name            NE          OD1
    res_name             ARG         ASP
    res_num.chain        498.E       38.A
    distance                         4.213

    id                   11046       163
    atom_name            NH2         OD1
    res_name             ARG         ASP
    res_num.chain        498.E       38.A
    distance                         4.027

    id                   11096       156
    atom_name            CE1         OE2
    res_name             HIS         GLU
    res_num.chain        505.E       37.A
    distance                         4.483




interaction for "rbm_h353":
    No interaction detected

interaction for "internal_ace2":
conf          1
    ELECTROSTATIC as_anion:

                         ace2_hotspot31 ace2_hotspot353

    id                   155         2724
    atom_name            OE1         NZ
    res_name             GLU         LYS
    res_num.chain        37.A        353.A
    distance                         3.550




interaction for "internal_rbm":
    No interaction detected
```

:::::

:::::{tab-item} ace2_spike_rbd_detailed.dd
:class-label: sd-shadow-sm

:::{card}
:img-top: ../../_static/omicron_ace2_visualized.png

**Figure 2.** The `ace2_spike_rbd_detailed.dd` visualized with Deemian Viewer.
:::

:::::

::::::

From the result above, we can see that there are no electrostatic interaction detected for RBM-Hotspot 353 pair and internal RBM.
On the other hand there are four electrostatic interaction identified in RBM and Hotspot 31 pair and one electrostatic interaction in the internal ACE2.

It is important to note that both ARG 493 and ARG 498 are originated from Q493R and Q498R mutation respectively.
These mutations are known to affect the binding affinity and contribute to Omicron variant immune evassiveness.
As for the internal ACE2 interaction, we can see that the LYS 353, which usually interact with the RBM now interact with the GLU 37 of ACE2.

This demonstration shows that using Deemian we can quickly identify both interprotein and intraprotein interaction.
Which are important in viral protein, protein engineering, and vaccine research.
We also see that in Deemian script we could chain selection using `and` keyword as in `protein and chain A and resid 5 to 30`.
