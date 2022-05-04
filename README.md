# RIP-MD3
RIP-MD to be used under python3.

To perform their function, proteins adopt a three dimensional structure that is mainly determined by non-covalent interactions. The structure of a protein is not static, residues may undergo conformational rearrangements, and in doing so, create, stabilize or break non-covalent interactions. Molecular Dynamics (MD) is a computational technique widely used to simulate these movements at atomic resolution. However, it is very difficult to gather relevant structural and functional information from MD simulations given the data-intensive nature of the technique. Consequently, several tools are often required to perform these highly complex analyses. Among these are Residue Interaction Networks (RINs), which have been used to facilitate the study of static protein structures. In a RIN, nodes represent Amino Acids (AAs) and the connections between nodes depict the non-covalent interactions that occur in the three dimensional structure of the protein. For this reason, we created RIP-MD, a tool that create RINs for a protein structure.

RIP-MD generates RINs for static protein structures or from MD trajectory files. The non-covalent interactions defined by RIP-MD include Hydrogen bonds, Salt bridges, Van der Waals contacts, cation-π, π-π, Arginine-Arginine and Coulomb interactions. In addition, RIP-MD also computesinteractions based on distances between C α s and disulphide bridges. Results of the analysis can be displayed in an user friendly interface. Moreover, the user can take advantage of the VMD visualization capacities, whereby through some effortless steps, it is possible to select and visualize interactions described for a single, several or all residues in a MD trajectory. Network and descriptive table files are also generated, allowing their further study in other specialized platforms. Furthermore RIP-MD generates correlation plots, where relationships between the dynamic behavior of different parts of the protein can be determined and quantified.

For more information, please read the user manual.

