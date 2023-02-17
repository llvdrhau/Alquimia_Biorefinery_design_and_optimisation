from pyomo.environ import *

from pyomo.environ import *

# Create an instance of the abstract model
model = AbstractModel()

# Define sets
model.I = Set()
model.J = Set()

# Define parameters
model.a = Param(model.I)
model.b = Param(model.J)

# Define variables
model.x = Var(model.I, domain=PositiveReals)

# Define constraints
def constraint_rule(model, j):
    return sum(model.a[i] * model.x[i] for i in model.I) >= 1
model.constraint = Constraint(model.J, rule=constraint_rule)

# Define objective function
def objective_rule(model):
    return sum(model.b[j] * model.x[j] for j in model.J)
model.objective = Objective(rule=objective_rule, sense=minimize)

# Create an instance of the model and set values for its components
instance = model.create_instance()
instance.I = [1, 2, 3]
instance.J = [1, 2]
instance.a.store_values = {1: 1.0, 2: 2.0, 3: 3.0}
instance.b.store_values = {1: 4.0, 2: 5.0}

instance.pprint()