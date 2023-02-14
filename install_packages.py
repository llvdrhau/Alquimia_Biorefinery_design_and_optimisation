import os
import sys

# Get the path to the directory where the script is located
script_dir = os.path.dirname(os.path.abspath(__file__))

# Set the name of the requirements.txt file
requirements_file = 'requirements.txt'

# Get the path to the Python interpreter that is currently running the script
python_path = sys.executable
# Get the path to the Python interpreter that is currently running the script
#python_path = sys.executable
#print(python_path)
#python_path = '/Users/lucasvanderhauwaert/Documents/Programing/Alquimia/venv'

# Get the path to the virtual environment directory
venv_dir = os.path.abspath(os.path.join(os.path.dirname(python_path), '..', ''))
venv_activate = os.path.join(venv_dir, 'bin', 'activate')

# Activate the virtual environment
os.system(f'source {venv_activate}')

# Generate a list of all installed packages and their versions
os.system(f'pip freeze > {script_dir}/{requirements_file}')

# Install all the required packages in the Python interpreter
os.system(f'pip install -r {script_dir}/{requirements_file}')


