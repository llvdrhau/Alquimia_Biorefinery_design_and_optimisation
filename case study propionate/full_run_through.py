"""
Run through the full case study calling all scripts that need to be called

Lucas Van der Hauwaert
"""
import subprocess

# start by fixing the models
subprocess.run(['python', r'C:\Users\lucas\PycharmProjects\Alquimia\SBML screening\Fixing_models\fix_other_propionibacteria_models.py'])

# make the surogate models
subprocess.run(['python', r'C:\Users\lucas\PycharmProjects\Alquimia\case study propionate\surrogate_models_SBML_v2.py'])

# make and solve the superstructure
subprocess.run(['python', r'C:\Users\lucas\PycharmProjects\Alquimia\case study propionate\superstructure_v2.py'])

