from jinja2 import Environment, PackageLoader


def generate_report(measurement, template_form):
    template_dict = dict(
        detailed_conf_first="01_detailed_conf_first.txt",
        detailed_type_first="02_detailed_type_first.txt",
        clustered_conf_first="03_clustered_conf_first.txt",
        clustered_type_first="04_clustered_type_first.txt",
        summarized_conf_first="05_summarized_conf_first.txt",
        summarized_type_first="06_summarized_type_first.txt",
    )

    template_file = template_dict[template_form]
    environment = Environment(loader=PackageLoader("deemian.writer", "templates"))
    template = environment.get_template(template_file)

    results = measurement.calculation_results
    interacting_subjects = measurement.interacting_subjects
    report = template.render(
        results=results, interacting_subjects=interacting_subjects, zip=zip, len=len, sorted=sorted
    )

    return report


def write_readable(report, out_file):
    with open(out_file, mode="w", encoding="utf-8") as result:
        result.write(report)
        print(f"... wrote {out_file}")
