# Installation

Before installing Deemian, please ensure that your machine has Python 3.10 or 3.11 installed.
If you don't have it, you can install Conda or Mamba, which are the easiest solutions.

After that, you need to create a virtual environment for Deemian installation so that your base Python environment remains clean.
If you have Python 3.10 or 3.11, you can use the venv module. Alternatively, you can use Conda or Mamba.

::::{tab-set}

:::{tab-item} venv
:sync: venv
:class-label: sd-shadow-sm

```bash
python -m venv venv
```
:::

:::{tab-item} Conda
:sync: conda
:class-label: sd-shadow-sm

```bash
conda create -n learn_deemian python=3.10
```
:::

:::{tab-item} Mamba
:sync: mamba
:class-label: sd-shadow-sm

```bash
mamba create -n learn_deemian python=3.10
```
:::

::::


Then activate the virtual environment:

::::::{tab-set}

:::::{tab-item} venv
:sync: venv
:class-label: sd-shadow-sm

::::{tab-set}

:::{tab-item} Linux/MacOS

```bash
source venv/bin/activate
```

:::

:::{tab-item} Windows cmd

```batch
venv\Scripts\activate.bat
```

:::

:::{tab-item} Windows PowerShell

```powershell
venv\Scripts\Activate.ps1
```

:::

::::

:::::

:::::{tab-item} Conda
:sync: conda
:class-label: sd-shadow-sm

```bash
conda activate learn_deemian
```
:::::

:::::{tab-item} Mamba
:sync: mamba
:class-label: sd-shadow-sm

```bash
mamba activate learn_deemian
```
:::::

::::::

To install Deemian, download the latest version from the [Deemian release page](https://github.com/radifar/deemian/releases).
Extract the zip / tar.gz file and run the `pip install` command, using the directory name as the last argument:

```bash
pip install deemian-0.1.0
```

Now, Deemian is ready to use in your virtual environment.
You can check out the [Getting Started](gettingstarted) or [Tutorial](tutorial) and follow the instructions to ensure that it is installed correctly.
