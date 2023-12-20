from collections import namedtuple

from lark import Lark, Tree, Token, Transformer


def parser(text: str, grammar: str = "grammar.lark") -> Tree:
    """Generate parser using Lark and grammar then parse the input text"""
    lark_parser = Lark.open(grammar, rel_to=__file__, parser="lalr")

    return lark_parser.parse(text)


class InstructionTransformer(Transformer):
    Selection = namedtuple("Selection", "name selection type", defaults=["selection"])
    BondCorrection = namedtuple("BondCorrection", "name template type", defaults=["bond_correction"])
    ReadableOutput = namedtuple("ReadableOutput", "results out_file type", defaults=["readable_output"])
    DeemianData = namedtuple("DeemianData", "out_file type", defaults=["deemian_data"])
    InteractionList = namedtuple("InteractionList", "interactions type", defaults=["interaction_list"])
    IncludeIonizable = namedtuple("IncludeIonizable", "charge boolean type", defaults=["include_ionizable"])
    Conformation = namedtuple("Conformation", "number type", defaults=["conformation"])
    ConformationRange = namedtuple("ConformationRange", "start end type", defaults=["conformation_range"])
    InteractingSubject = namedtuple(
        "InteractingSubject", "subject_1 subject_2 name type", defaults=["interacting_subject"]
    )

    def selection(self, args):
        selection_phrase = []
        for arg in args:
            if isinstance(arg, Tree):
                selection_phrase.append(arg.data)
                selection_phrase.extend([child.value for child in arg.children])
            elif isinstance(arg, Token):
                selection_phrase.append(arg.value)

        return tuple(selection_phrase)

    def combine_selection(self, args):
        selection_phrase = []
        for arg in args:
            if isinstance(arg, Tree):
                selection_phrase.append(arg.data)
                selection_phrase.extend([child.value for child in arg.children])
            elif isinstance(arg, Token):
                selection_phrase.append(arg.value)

        return tuple(selection_phrase)

    def assign_selection(self, args):
        name = args[0].value
        selection_chain = []
        for arg in args[1:]:
            selection_chain.append(arg)

        return self.Selection(name, selection_chain)

    def bond_correction(self, args):
        name = args[0].value
        template = args[1].value

        return self.BondCorrection(name, template)

    def interaction_list(self, args):
        return self.InteractionList([arg.value for arg in args])

    def include_ionizable(self, args):
        return self.IncludeIonizable(args[0].value, args[1].value)

    def interacting_subject(self, args):
        subject_1 = args[1].value
        subject_2 = args[3].value
        name = subject_1 + ":" + subject_2

        return self.InteractingSubject(subject_1, subject_2, name)

    def interacting_subject_with_alias(self, args):
        subject_1 = args[1].value
        subject_2 = args[3].value
        name = args[4].value

        return self.InteractingSubject(subject_1, subject_2, name)

    def conformation(self, args):
        conformation_list = [int(arg.value) for arg in args]
        return self.Conformation(conformation_list)

    def conformation_range(self, args):
        return self.ConformationRange(args[0].value, args[1].value)

    def readable_output(self, args):
        results = [arg.value for arg in args[:-1]]
        out_file = args[-1].value

        return self.ReadableOutput(results, out_file)

    def deemian_data(self, args):
        return self.DeemianData(args[0].value)
