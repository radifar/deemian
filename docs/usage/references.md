# References

## Step, Instruction, and Keyword

In Deemian DSL, the scripting is segmented into three **step**s, where each step acts as a container to execute a series of instructions.
If you are familiar with programming, step is similar to class instantiation in object-oriented programming.
A step is usually structured like this:

```text
step_name step_identifier [
    instruction_1
    instruction_2
    ...
]
```

`step_name` is the name of the step (`molecule`, `measure`, or `present`).
While `step_identifier` is the identifier of the step, it is similar to `variable_name` in programming.
Defining the right `step_identifier` is important, and it will be explained in detail in each step below.

The `step` contains a series of instructions enclosed in a square bracket.
The spaces before each instruction are optional but recommended to make the script more readable.
As a rule of thumb, use four spaces before each instruction.

Instruction consists of instruction keyword, identifier, and variable.
The structure/syntax of each instruction is different and will be explained in detail below.

### 1. molecule step

`molecule` step act as a container for molecule processing, which cover molecule file reading, molecule selection, and molecule bond assignment.
The `step_identifier` for `molecule` step correspond to the file name that would be opened.
Therefore it is possible to open more than one file using multiple `molecule` step, which is useful for virtual screening where multiple ligands and poses were paired with single target.

#### select `selection_id` = `selection_string`

`select` keyword is used to create a new selection using `selection_string` and then assigned the selection to `selection_id`, which could be used for further processing (bond assignment) or directly used as interacting subject in `measure` step later.

The `selection_string` is either a single selection or a chain of selections.
A chain of selections consist of multiple selection connected with boolean keyword such as **and**, **or**, and **not**.
At the moment, Deemian only support **and** and **not** keyword.

A selection could be a keyword-value pair or a macro.
Keyword value pair is consisted of a keyword, such as chain, resname, and resid, which then followed by its acceptable values.
The acceptable values for each keyword is different:

- **chain**: upper case single letter such as 'A', 'B', or 'C', and it has to match the chain you are about to select from the molecule file being read. Example `chain A`
- **resname**: upper case three letter such as 'ALA', 'COX', 'G39'. Should match the resname of the residue name you are about to select from the file being read. Example `resname G39`
- **resid**: positive integer, which designate the residue id/number that being selected. At the moment it only accept the resid in the form of resid range, for example: `resid 5 to 30`.

Macro is shorthand for a molecule category, such as `protein`, `dna`, `hydrophobic`, `sugar`, `ion`, and more. But at the moment Deemian only accept `protein` macro.

By combining the the keyword-value pair and macro we could construct a `selection string` such as:
- `not protein`
- `protein and chain A`
- `chain A and resname COX`
- `chain A and resid 30 to 50`

#### assign bond `selection_id` template `smile_string`

`assign bond` keyword is used for bond assignment task in case the bond order of a molecule is not properly recognized.
This is usually happen when the molecule file being read do not carry hydrogen coordinates or connectivity information.
`assign bond` keyword then can be used to modify a selection (usually a small molecule) then use a template in the form of SMILES to correct the bond orders.

An example of `assign bond` usage: `assign bond oseltamivir template CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)N)C(=O)O`

### 2. measure step

`measure` step act as a container for molecule interaction analysis configuration, which cover interaction type selection, interaction type customization, selection pairing, and conformer selection.
The `step_identifier` for `measure` is arbitrary, later the present step will look for the same `step_identifier` to retrieve the interaction analysis result and present it in the form of readable output and deemian data output.

#### interactions `options`

`interactions` keyword used to select the interaction type included in the analysis.
The possible options for `interactions` keyword are `all`, `electrostatic`, `nonpolar`, `hydorgen_bond`, and `pi`.
At the moment Deemian only capable of doing electrostatic analysis, therefore the only acceptable options are `all` and `electrostatic`.

An example `interactions` usage: `interactions all` and `interactions electrostatic`

#### ionizable `charge` `boolean`

`ionizable` keyword used to customize the charge identification which would later be used for electrostatic interaction analysis.
The default value for charge identification is using apparent charge, which means only deprotonated acid, protonated base, and ions are regarded as charged atoms or groups.
It is possible to treat potentially ionizable groups as charged groups by using the `ionizable` instruction.
For example, using `ionizable positive true` will treat neutral amine, histidine, guanidine moieties as positively charged moieties.
While using `ionizable negative true` will treat neutral carboxylic, phosphoric, and sulphuric moieties as negatively charged moieties.

#### between `subject_1` and `subject_2` [as `identifier`]

`between ___ and ___` keywords used to create the selection pair that will be analyzed.
Here the first value in the blank will be referred as subject_1 and the second value as subject_2.
The sequence is important as it will be reflected in the resulting presentation.
The rule of thumb is to put the subject of interest as first value.
For example if we are analyzing the interaction of oseltamivir and neuraminidase, it is better to put oseltamivir as first value like so `between oseltamivir and neuraminidase`.

The resulting interaction will be refered using the combination of subject_1 and subject_2 identifier in the resulting presentation.
For example, using the `between oseltamivir and neuraminidase` instruction will generate an interaction referred as `oseltamivir:neuraminidase`.
To use custom interaction identifier, use the optional `as identifier` keywords, for example: `between oseltamivir and neuraminidase as oseltamivir_N1`.

#### conformation `start` to `end`

Sometimes we analyze multiconformer/multimodel structures, in this case we can select which conformer/model we would like to analyze using the `conformation start to end` instruction, where `start` is the starting conformer number, and `end` is the last conformer number.
For example to analyze the first to twentieth conformer/model we can use the following instruction: `conformation 1 to 20`.

### 3. present step

`present` step acts as a container for result presentation, which include the presentation format for readable output, and the output file name.
The `step_identifier` should match the `measure` `step_identifier` it consumes.
When the `present` step find the `step_identifier` it matches, it will retrieve the analysis data from the corresponding `measure` step and execute the instructions contained within the step to extract and present the data as instructed.

#### interactions [`output_format`] `file_name`

`interactions` keyword in `present` step will generate the formatted readable output of interaction analysis.
This keyword followed by optional `output_format` and mandatory `file_name` value.
By default the `output_format` is `detailed_conf_first`, currently the available formats are:
- `detailed_conf_first`
- `detailed_type_first`
- `clustered_conf_first`
- `clustered_type_first`
- `summarized_conf_first`
- `summarized_type_first`

#### deemiandata `file_name`

`deemiandata` keyword is used to generate the Deemian data which can be visualized with Deemian Viewer.
