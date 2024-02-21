# Roadmap

Deemian is still in its early phase and thus only has minimum features.
However, there are many features planned, including:

- **Basic interactions**: Hydrophobic, hydrogen bond, electrostatic interaction, aromatic interaction.
- **Other interactions**: Polar interaction (loose H-bond), pi-cation interaction, halogen bond, carbonyl interaction, metal complex/coordination, etc.
- **Bridged interactions**: Water-bridged interaction and metal-bridged interaction.
- **Anti-interactions**: Steric clash interaction, same-charge interaction, polar-non polar interaction.
- **Element and atom type interactions**: Calculate the interaction between elements or Sybil atom type [(Ballester, Schreyer, and Blundel, 2014)](https://dx.doi.org/10.1021/ci500091r).
- **Bitstring AKA Bitvector**: The ability to present interaction as bitstring which allows for interaction comparison. Also add the ability to calculate the interaction similarity.
- **Better compatibility with Docking results**: Create commands to make virtual screening result analysis easier.
- **Better compatibility with MD results**: Create commands to make MD result analysis easier.
- **Interaction Pseudo-Atom (IPA)**: Present interaction as [IPA (Desaphy et.al., 2013)](https://dx.doi.org/10.1021/ci300566n) with three kind of placement, at first subject, at the middle of interaction, and at second subject.
- **Hashing capability**: Integrating hashing and scoring to enable [TIFP (Desaphy et. al., 2013)](https://dx.doi.org/10.1021/ci300566n) and [PLEC (Wójcikowski et. al., 2018)](https://dx.doi.org/10.1093/bioinformatics/bty757).
- **Structure alignment**: Add the capability to align 3D structure. Useful for accumulating IPA from different protein-ligand complexes, protein family (e.g. Kinase family) comparison, virus mutation study.
- **Group/substructure interaction**: Allows interaction analysis between chemical groups or large substructure.
- **Interacting surface area**: Calculate the interacting surface area for each interaction type.
- **Integration with Jupyter**: Allows the Deemian data to be processed and analyzed within Jupyter Notebook/Lab. Also allows visualization with NGLView.
- **Extension with Python**: Allows user to create custom instruction via Python scripting.
- **Support various file format**: Theoretically any file format that is supported by RDKit and NGL can be loaded.
- **Support CGMD results**: especially Martini 3, which gain a lot of traction lately and useful in large scale simulation.
