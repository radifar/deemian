from lark import Tree

from deemian.engine.builder import DeemianData


def director(steps: Tree, data: DeemianData):
    for step in steps.children:
        if step.data == "molecule_op":
            instructions = step.children
            mol_filename = instructions.pop(0).value
            data.add_molecule(mol_filename)

            for inst in instructions:
                if inst.type == "selection":
                    data.add_selection(inst.name, inst.selection, mol_filename)

                elif inst.type == "bond_correction":
                    data.correct_bond(inst.name, inst.template)

        elif step.data == "measurement":
            instructions = step.children
            measurement_id = instructions.pop(0).value
            measurement = data.add_measurement(measurement_id)

            for inst in instructions:
                if inst.type == "interaction_list":
                    measurement.interactions.extend(inst.interactions)

                elif inst.type == "include_ionizable":
                    measurement.set_ionizable(inst.charge, inst.boolean)

                elif inst.type == "interacting_subject":
                    measurement.interacting_subjects[inst.name] = (inst.subject_1, inst.subject_2)

                elif inst.type == "conformation":
                    measurement.conformation.extend(inst.number)

                elif inst.type == "conformation_range":
                    measurement.set_conformation_range(inst.start, inst.end)

            data.calculate_interactions(measurement_id)

        elif step.data == "presentation":
            instructions = step.children
            presentation_id = instructions.pop(0).value

            for inst in instructions:
                if inst.type == "interaction_output":
                    data.write_readable_output(presentation_id, inst.out_file, inst.format)

                elif inst.type == "deemian_data":
                    data.write_deemian_data(presentation_id, inst.out_file)
