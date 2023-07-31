"""
Run through the full case study (V2) calling all scripts that need to be called

Lucas Van der Hauwaert
"""
import subprocess
import sys

# List of script filenames to run
script_filenames = [
# start by fixing the models
r'C:\Users\lucas\PycharmProjects\Alquimia\SBML screening\Fixing_models\fix_other_propionibacteria_models.py',

# make the surrogate models for the pure culture reactors
r'C:\Users\lucas\PycharmProjects\Alquimia\case study propionate\surrogate_models_SBML_v2.py',

# make the surrogate models for the mixed culture reactor
r'C:\Users\lucas\PycharmProjects\Alquimia\case study propionate\surrogate_model_open_fermentation.py',

# make and solve the superstructure
r'C:\Users\lucas\PycharmProjects\Alquimia\case study propionate\superstructure_v2.py'

]

# Run each script sequentially
for script_filename in script_filenames:
    subprocess.run([sys.executable, script_filename], check=True, cwd=sys.path[0])

