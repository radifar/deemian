import typer
from rich.console import Console

from deemian.engine.builder import DeemianData
from deemian.engine.processor import parser, InstructionTransformer
from deemian.engine.director import director


console = Console()
app = typer.Typer()


@app.command(short_help="run deemian script")
def run(script_name: str):
    console.print(f"[bold magenta]Running {script_name}[/bold magenta]")
    with open(script_name) as f:
        text = f.read()

    command_tree = parser(text)
    transformer = InstructionTransformer()
    command_tree = transformer.transform(command_tree)

    deemian_data = DeemianData()
    director(command_tree, deemian_data)


# one command one callback, just a temporary helper command
# delete this when there are multiple commands
# https://typer.tiangolo.com/tutorial/commands/one-or-multiple/
@app.callback()
def callback():
    pass
