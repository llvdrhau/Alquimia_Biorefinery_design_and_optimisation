"""
Run through the full case study (V2) calling all scripts that need to be called

Lucas Van der Hauwaert
"""
import subprocess
import sys


# get the working directory
basePath = sys.path[0]
# find the word 'Alquimia' in the path and split the path at this point
# to get the path to the Alquimia folder
basePath = basePath.split('Alquimia')[0] + 'Alquimia\\'

pathFixingModels = basePath + 'SBML screening\\Fixing_models\\fix_all_propionibacteria_models.py'
pathSurrogateModles  = basePath + 'case study propionate\\surrogate_models_SBML_v2.py'
pathSurrogateModelOpenFermentation = basePath + 'case study propionate\\surrogate_model_open_fermentation.py'
pathSuperstructure = basePath + 'case study propionate\\superstructure_v2.py'


# List of script filenames to run
script_filenames = [
pathFixingModels,     # start by fixing the models
pathSurrogateModles,   # make the surrogate models for the pure culture reactors
pathSurrogateModelOpenFermentation, # make the surrogate models for the mixed culture reactor
pathSuperstructure # make and solve the superstructure
]

# Run each script sequentially
for script_filename in script_filenames:
    subprocess.run([sys.executable, script_filename], check=True, cwd=sys.path[0])

