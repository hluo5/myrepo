import pandas as pd
import numpy as np
from lmfit import Model,Parameters, minimize, report_fit, fit_report
from scipy import stats
import matplotlib.pyplot as plt

#define the set of equations without H (1-O, 2-Mg, 3-Si)
def equation1(x, a, b, c):
    pres, temp = x
    logD = a + b/temp + c*(pres/temp)
    #logD = a + b*np.exp(pres/c)
    return logD

#load data
data_1 = pd.read_excel('D_data.xlsx')
column_names = data_1.columns
x_data_1 = np.array([
    data_1['pres'], 
    data_1['temp'] 
])
y_data_1 = np.array([
    data_1['logD_ThO2'] 
])

def model_1(x, a, b, c):
    return equation1(x, a, b, c)

model_1 = Model(model_1)

# build parameters
params = Parameters()
params.add('a', value=0.0)
params.add('b', value=0.0, vary=False)
#params.add('b', value=0.0)
#params.add('c', value=1.0, vary=False)
params.add('c', value=0.0)

result = model_1.fit(y_data_1, params, x=x_data_1)

print(report_fit(result))
report_text = fit_report(result)
with open("fit_report.txt", "w") as f:
    f.write(report_text)

#check P values
parameters = result.params
p_values = []
for param_name, param_value in parameters.items():
    if param_value.value != 0.0:
        stderr = param_value.stderr
        test_statistic = param_value.value / stderr
        p_value = 2 * (1 - stats.norm.cdf(abs(test_statistic)))
        p_values.append((param_name, p_value))
    else:
        p_values.append((param_name, None))
for param_name, p_value in p_values:
    print(f'{param_name}: p-value = {p_value}')
with open("fit_report.txt", "a") as f:
    for value in p_values:
        f.write(str(value) + "\n")

#plot out the predicted data vs. input data
y_pred_1 = model_1.eval(params=result.params, x=x_data_1)

plt.figure(figsize=(10, 5))
plt.scatter(y_data_1, y_pred_1)
plt.plot([np.min(y_data_1), np.max(y_data_1)], [np.min(y_data_1), np.max(y_data_1)], 'k--', lw=2)
plt.xlabel('Actual y_data_1')
plt.ylabel('Predicted y_data_1')
plt.title('Model 1')

plt.tight_layout()
plt.savefig('my_plot.pdf', dpi=600)
plt.show()

# save data
y_data_1 = y_data_1.reshape(-1, 1)
y_pred_1 = y_pred_1.reshape(-1, 1)
data1_to_save = np.column_stack((y_data_1, y_pred_1))
np.savetxt('y_data_1_and_pred_1.txt', data1_to_save)
