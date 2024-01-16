# Why Deemian DSL?

Molecular interaction analysis is essential in many research areas, including drug discovery, vaccine design, molecular genetic study, antibacterial/antiviral agent resistance mechanism, enzyme engineering, viral organ/host tropism, and more. With myriad use cases in molecular interaction analysis, currently, no tool can deal with all of those cases. Frequently, A molecular interaction analysis tool is designated for a specific research area, such as PDB structure analysis (e.g., Arpeggio, Protein Structure Atlas, PLIP, and FingeRNAt), structure-based virtual screening (e.g., IChem and PyPLIF HIPPOS), and dynamic binding study using MD (e.g., MD-IFP).

One possible approach to cover a wide range of use cases is by giving the user the power to specify the interaction analysis via a programming language. Such an approach is doable using a library/framework for molecular interaction analysis such as ProLIF. However, there is a caveat as it requires the user to have a certain level of mastery in Python language. And this is where Deemian enters the stage.

Deemian is a Domain Specific Language (DSL) for molecular interaction analysis. Therefore, it allows the user to specify the instructions for molecule file reading, molecule selection, preparation, molecular interaction configuration, and result presentation. Deemian could empower computational chemists or bioinformaticians to do in-depth molecular interaction analysis by providing users with simple and intuitive commands.

The other benefits of using a DSL for molecular interaction analysis are better readability, reproducibility, and reusability. As it is easy to read, understand, replicate, and modify Deemian script from another study.
