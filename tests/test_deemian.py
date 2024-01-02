from unittest.mock import patch

from typer.testing import CliRunner

from deemian.deemian import app


runner = CliRunner()


def test_app():
    with patch("deemian.deemian.director"):
        result = runner.invoke(app, ["run", "tests/data/n1-oseltamivir-5nzn.txt"])
        assert result.exit_code == 0
