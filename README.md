## DiffEx_SCA1_SCA2

This repository includes the jupyter notebook [`notebooks/DiffEx_SCA1_SCA2.ipynb`](https://github.com/chopravi/DiffEx_SCA1_SCA2/blob/master/notebooks/DiffEx_SCA1_SCA2.ipynb) which contains code for analyzing ion channel gene dysregulation overlap between mouse models of spinocerebellar ataxia and for grouping channel genes according to their regulation by an Atxn1 Cic binding mutant.

This repository also includes the `diffex` package, which support analysis in the jupyter notebook.

In order to run the jupyter notebook, you'll need to have [Python 3](https://www.python.org/downloads/) installed, then follow the installation and usage instructions below.

### Installation

These instructions are meant for Windows powershell, but if you're a bash user it should be straightforward to translate.

Open a powershell prompt, paste in the instructions below and press enter to run.

```sh
# Copy this repository into a folder called
# 'Diffex_SCA1_SCA2' in your current directory
git clone https://github.com/chopravi/DiffEx_SCA1_SCA2.git
# Make sure virtual env is installed
pip install virtualenv
# Move to the Diffex_SCA1_SCA2 folder and create the .diffex virtual environment
cd .\DiffEx_SCA1_SCA2
python -m venv .diffex
# Activate the .diffex environment
.\.diffex\Scripts\activate
# Install requirements to .diffex environment
pip install -r requirements.txt
# Add .diffex to the list of envs that jupyter notebook run kernels in
python -m ipykernel install --user --name=.diffex
# Install this package to the .diffex environment
python setup.py install
```

### Usage

Once you've completed installation, you'll be able to run the analysis [notebook](https://github.com/chopravi/DiffEx_SCA1_SCA2/blob/master/notebooks/DiffEx_SCA1_SCA2.ipynb) on your machine.

To open the notebook, make sure you are in the `DiffEx_SCA1_SCA2` source directory (which you still will be if you haven't exited powershell after the installation step above). Also, make sure the `.diffex` virtual environment is activated as written above.

```sh
.\.diffex\Scripts\activate
```

Then, navigate to the notebooks directory and run the notebook

```sh
cd .\notebooks
jupyter notebook .\DiffEx_SCA1_SCA2.ipynb
```

This will open the notebook in your default web browser. Make sure that the notebook is running in the `.diffex` kernel. You can check by looking at the upper righthand corner of the window below 'Logout'. If it doesn't say '.diffex', you can switch to the `.diffex` kernel by selecting Kernel > Change Kernel > .diffex. If '.diffex' does not appear in that menu, make sure you've run the installation step above. 

You can then run the all code in the notebook by selecting Cell > Run All. This will output several plots and .csv files that are already included in the repository.

For additional questions regarding this code, please contact authors below:

Ravi Chopra (chopra.r@wustl.edu)
John Cooper (jpcoope@utexas.edu)
