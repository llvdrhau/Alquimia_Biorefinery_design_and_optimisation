import pandas as pd
import keras
import numpy as np

# df_features = pd.read_excel('Gelatinedata.xlsx', sheet_name = 'features')
# df_outputs = pd.read_excel('Gelatinedata.xlsx', sheet_name = 'outputs')

model = keras.models.load_model('NN_gelatine')
a = model.summary()
print(a)
inputs = np.array([10,4,6])
inputs = inputs.reshape([1,-1])
outputs = model.predict(inputs)
print(outputs)


#make a test to see if the plots of the test are the same as the one where the model was trained but on the right scale.

