# Roadmap

Deemian is still in its early phase and thus only has minimum features.
However, there are many features planned, including:

- **Basic interactions**: hydrophobic, hydrogen bond, electrostatic interaction, aromatic interaction.
- **Other interactions**: polar interaction (weak H-bond), pi-cation interaction, halogen bond, metal complex/coordination, etc.
- **Bridged interactions**: water-bridged interaction and metal-bridged interaction.
- **Bitstring AKA Bitvector**: The ability to present interaction as bitstring which allows for interaction comparison. Also add the ability to calculate the interaction similarity.
- **Interaction Pseudo-Atom (IPA)**: Present interaction as IPA with three kind of placement, at first subject, at the middle of interaction, and at second subject.
- **Structure alignment**: Add the capability to align 3D structure. Useful for accumulating IPA from different protein-ligand complexes, protein family (e.g. Kinase family) comparison, virus mutation study.
- **Integration with Jupyter**: Allows the Deemian data to be processed and analyzed within Jupyter Notebook/Lab. Also allows visualization with NGLView too.
- **Extension with Python**: Allows creating custom instruction via Python scripting.
- **Support various file format**: Theoretically any file format that is supported by RDKit and NGL can be loaded.
- **Better compatibility with Docking results**: Create special commands to make virtual screening result analysis easier.
- **Better compatibility with MD results**: Create special commands to make MD result analysis easier.
- **Support CGMD results**: especially Martini 3, which gain a lot of traction lately and useful in large scale simulation.
