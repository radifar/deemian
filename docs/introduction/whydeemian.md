# Why Deemian DSL?

Molecular interaction analysis is a critical aspect of various research fields such as drug discovery, vaccine design, enzyme engineering, and more.
However, there is no one-size-fits-all tool that can handle all types of analyses.
These tools are usually designed for specific research areas such as PDB structure analysis, structure-based virtual screening, or dynamic binding studies using MD.

One way to address this issue is to allow users to specify the interaction analysis through programming.
A library/framework for molecular interaction analysis, such as ProLIF, can achieve this.
However, using such a library requires a certain level of mastery in the Python programming language, which can be challenging for some users.

This is where Deemian comes in.
Deemian is a Domain-Specific Language (DSL) for molecular interaction analysis.
In general, using DSLs offers several advantages, such as:

- **Expressiveness**: DSLs are tailored to the specific needs and vocabulary of a particular domain, making them more expressive and concise for tasks within that domain. This can lead to clearer and more readable code.
- **Abstraction**: DSLs allow you to abstract away unnecessary details, focusing on the key aspects of your problem domain. This can simplify development and improve code maintainability.
- **Productivity**: DSLs can enhance productivity by providing specialized constructs and syntax that align closely with the problem you’re solving, reducing the amount of boilerplate code needed.
- **Collaboration**: DSLs can bridge the communication gap between domain experts and developers. Non-programmers can often understand and contribute to DSL code, fostering better collaboration between technical and non-technical team members.

In simple terms, Deemian is a DSL for molecular interaction analysis, which allows users to specify instructions for molecule file reading, molecule selection, preparation, molecular interaction configuration, and result presentation.
It is designed to make in-depth molecular interaction analysis easier for computational chemists or bioinformaticians by providing simple and intuitive commands.
Additionally, Deemian scripts are easy to read, replicate, and modify, which makes them more reusable and reproducible.
