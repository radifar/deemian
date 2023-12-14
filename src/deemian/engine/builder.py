from dataclasses import dataclass, field


@dataclass
class DeemianData:
    selections: dict = field(default_factory=lambda: {})
    interactions: list = field(default_factory=lambda: [])
    ionizable: dict = field(default_factory=lambda: {"positive": False, "negative": False})
    interacting_subjects: dict = field(default_factory=lambda: {})
    conformation: list = field(default_factory=lambda: [])
    interaction_details: dict = field(default_factory=lambda: {})
    readable_output: dict = field(default_factory=lambda: {})


class DeemianDataBuilder:
    def __init__(self, deemian_data: DeemianData) -> None:
        self.deemian_data = deemian_data

    def assign_selection(self, name: str, selection: list[tuple], molecule: str):
        selection.insert(0, molecule)
        self.deemian_data.selections[name] = selection

    def correct_bond(self, name: str, template: str):
        return (name, template)

    def set_interactions(self, interactions: list[str]):
        self.deemian_data.interactions = interactions

    def set_ionizable(self, charge: str, boolean: str):
        if boolean == "true":
            boolean = True
        elif boolean == "false":
            boolean = False
        self.deemian_data.ionizable[charge] = boolean

    def set_interacting_subjects(self, subject_1: str, subject_2: str, name: str):
        self.deemian_data.interacting_subjects[name] = (subject_1, subject_2)

    def set_conformation(self, number: str):
        self.deemian_data.conformation = [int(number)]

    def set_conformation_range(self, start: str, end: str):
        self.deemian_data.conformation = list(range(int(start), int(end) + 1))

    def calculate_interactions(self):
        return 1

    def generate_readable_output(self):
        return 1

    def generate_deemian_data(self):
        return self.deemian_data
