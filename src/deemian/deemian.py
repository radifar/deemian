import time

import typer
from rich.console import Console

from deemian import __version__ as deemian_version
from deemian.engine.builder import DeemianData
from deemian.engine.processor import parser, InstructionTransformer
from deemian.engine.director import director


console = Console()
app = typer.Typer()


@app.command(short_help="run deemian script")
def run(script_name: str):
    now = time.time()
    console.print(f"[bold deep_pink3]   Running {script_name}[/bold deep_pink3]")
    print(f"   on Deemian version: {deemian_version}")
    with open(script_name) as f:
        text = f.read()

    command_tree = parser(text)
    transformer = InstructionTransformer()
    command_tree = transformer.transform(command_tree)

    deemian_data = DeemianData()
    deemian_data.metadata.start_time = now
    director(command_tree, deemian_data)
    console.print(f"\nrunning time: {deemian_data.metadata.running_time}\n")


# one command one callback, just a temporary helper command
# delete this when there are multiple commands
# https://typer.tiangolo.com/tutorial/commands/one-or-multiple/
@app.callback()
def callback():
    pass
