## DiffEx_SCA1_SCA2

Python 3 code for analyzing gene dysregulation overlap between disease models

### Installation

These instructions are specifc to Windows, so run them in powershell.

```sh
git clone https://github.com/chopravi/DiffEx_SCA1_SCA2.git
# Make sure virtual env is installed
pip install virtualenv
# Move the the project folder and create the .diffex environment
cd .\DiffEx_SCA1_SCA2
python -m venv .diffex
# Activate the .diffex environment
.\.diffex\Scripts\activate
# Install requirements to .diffex environment
pip install -r requirements.txt
# Add .diffex to the list of envs that jupyter notebook can access
python -m ipykernel install --user --name=.diffex
# Install this package to the .diffex environment
python setup.py install
```

### Usage

All functionality is found in `notebooks/DiffEx_SCA1_SCA2.ipynb`. Make sure to run the line `python -m ipykernel install --user --name=.diffex` and make sure your jupyter notebook kernel is running in `.diffex`

For additional questions regarding this code, please contact authors below:

Ravi Chopra (chopra.r@wustl.edu)
John Cooper (jpcoope@utexas.edu)
