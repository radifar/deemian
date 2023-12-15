from lark import Tree

from deemian.engine.builder import DeemianDataBuilder


def director(steps: Tree, data_builder: DeemianDataBuilder):
    for step in steps.children:
        if step.data == "molecule_op":
            molecule_op_instructions = step.children
            mol_filename = molecule_op_instructions[0].value
            data_builder.read_molecule(mol_filename)

            for instruction in molecule_op_instructions[1:]:
                if instruction.type == "selection":
                    name = instruction.name
                    selection = instruction.selection
                    data_builder.assign_selection(name, selection, mol_filename)

                elif instruction.type == "bond_correction":
                    name = instruction.name
                    template = instruction.template
                    data_builder.correct_bond(name, template)

        elif step.data == "measurement":
            measurement_instructions = step.children
            measure_identifier = measurement_instructions[0].value

            for instruction in measurement_instructions[1:]:
                if instruction.type == "interaction_list":
                    data_builder.set_interactions(instruction.interactions)

                elif instruction.type == "include_ionizable":
                    charge = instruction.charge
                    boolean = instruction.boolean
                    data_builder.set_ionizable(charge, boolean)

                elif instruction.type == "interacting_subject":
                    subject_1 = instruction.subject_1
                    subject_2 = instruction.subject_2
                    name = instruction.name
                    data_builder.set_interacting_subjects(subject_1, subject_2, name)

                elif instruction.type == "conformation":
                    data_builder.set_conformation(instruction.number)

                elif instruction.type == "conformation_range":
                    start = instruction.start
                    end = instruction.end
                    data_builder.set_conformation_range(start, end)

            else:
                data_builder.calculate_interactions(measure_identifier)
        elif step.data == "presentation":
            presentation_instructions = step.children
            presentation_identifier = presentation_instructions[0].value

            for instruction in presentation_instructions[1:]:
                if instruction.type == "readable_output":
                    data_builder.write_readable_output(instruction.out_file, presentation_identifier)

                elif instruction.type == "deemian_data":
                    data_builder.write_deemian_data(instruction.out_file, presentation_identifier)
