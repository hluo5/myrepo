import pandas as pd
import numpy as np
from lmfit import Model,Parameters, minimize, report_fit, fit_report
from scipy import stats
import matplotlib.pyplot as plt

#define the set of equations without H (1-O, 2-Mg, 3-Si)
#def equation1(x, a21, b21, c21, oo1, omg1, osi1):
#    pres, temp, xo, xmg, xsi = x
#    lnGamma_O = -oo1*1873/temp*np.log(1-xo) - (omg1*1873/temp*xmg*(1+np.log(1-xmg)/xmg-1/(1-xo)) + osi1*1873/temp*xsi*(1+np.log(1-xsi)/xsi-1/(1-xo))) + (omg1*1873/temp*xmg**2*xo*(1/(1-xo)+1/(1-xmg)+xo/(2*(1-xo)**2)-1) + osi1*1873/temp*xsi**2*xo*(1/(1-xo)+1/(1-xsi)+xo/(2*(1-xo)**2)-1))
#    lnKo = a21 + b21/temp + c21*pres/temp - lnGamma_O
#    return lnKo

#def equation2(x, a31, b31, c31, mgmg1, omg1, mgsi1, oo1, osi1):
#    pres, temp, xo, xmg, xsi = x
#    lnGamma_O = -oo1*1873/temp*np.log(1-xo) - (omg1*1873/temp*xmg*(1+np.log(1-xmg)/xmg-1/(1-xo)) + osi1*1873/temp*xsi*(1+np.log(1-xsi)/xsi-1/(1-xo))) + (omg1*1873/temp*xmg**2*xo*(1/(1-xo)+1/(1-xmg)+xo/(2*(1-xo)**2)-1) + osi1*1873/temp*xsi**2*xo*(1/(1-xo)+1/(1-xsi)+xo/(2*(1-xo)**2)-1))
#    lnKmg = a31 + b31/temp + c31*pres/temp - lnGamma_O - (-mgmg1*1873/temp*np.log(1-xmg) - (omg1*1873/temp*xo*(1+np.log(1-xo)/xo-1/(1-xmg)) + mgsi1*1873/temp*xsi*(1+np.log(1-xsi)/xsi-1/(1-xmg))) + (omg1*1873/temp*xo**2*xmg*(1/(1-xmg)+1/(1-xo)+xmg/(2*(1-xmg)**2)-1) + mgsi1*1873/temp*xsi**2*xmg*(1/(1-xmg)+1/(1-xsi)+xmg/(2*(1-xmg)**2)-1)))
#    return lnKmg

#def equation3(x, a41, b41, c41, sisi1, osi1, mgsi1, oo1, omg1):
#    pres, temp, xo, xmg, xsi = x
#    lnKsi = a41 + b41/temp + c41*pres/temp - (-sisi1*1873/temp*np.log(1-xsi) - (osi1*1873/temp*xo*(1+np.log(1-xo)/xo-1/(1-xsi)) + mgsi1*1873/temp*xmg*(1+np.log(1-xmg)/xmg-1/(1-xsi))) + (osi1*1873/temp*xo**2*xsi*(1/(1-xsi)+1/(1-xo)+xsi/(2*(1-xsi)**2)-1) + mgsi1*1873/temp*xmg**2*xsi*(1/(1-xsi)+1/(1-xmg)+xsi/(2*(1-xsi)**2)-1)))
#    return lnKsi


#define the set of equations with H (4-H, 5-O, 6-Mg, 7-Si)
def equation4(x, a1, b1, c1, hh, ho, hmg, hsi, oo, omg, osi):
    pres, temp, xh, xo, xmg, xsi = x
    lnKh = a1 + b1/temp + c1*pres/temp - 2*(-hh*1873/temp*np.log(1-xh) - (ho*1873/temp*xo*(1+np.log(1-xo)/xo-1/(1-xh)) + hmg*1873/temp*xmg*(1+np.log(1-xmg)/xmg-1/(1-xh)) + hsi*1873/temp*xsi*(1+np.log(1-xsi)/xsi-1/(1-xsi))) + (ho*1873/temp*xo**2*xh*(1/(1-xh)+1/(1-xo)+xh/(2*(1-xh)**2)-1) + hmg*1873/temp*xmg**2*xh*(1/(1-xh)+1/(1-xmg)+xh/(2*(1-xh)**2)-1) + hsi*1873/temp*xsi**2*xh*(1/(1-xh)+1/(1-xsi)+xh/(2*(1-xh)**2)-1)))
    return lnKh

def equation5(x, a2, b2, c2, oo, ho, omg, osi):
    pres, temp, xh, xo, xmg, xsi = x
    lnGamma_O = -oo*1873/temp*np.log(1-xo) - (ho*1873/temp*xh*(1+np.log(1-xh)/xh-1/(1-xo)) + omg*1873/temp*xmg*(1+np.log(1-xmg)/xmg-1/(1-xo)) + osi*1873/temp*xsi*(1+np.log(1-xsi)/xsi-1/(1-xo))) + (ho*1873/temp*xh**2*xo*(1/(1-xo)+1/(1-xh)+xo/(2*(1-xo)**2)-1) + omg*1873/temp*xmg**2*xo*(1/(1-xo)+1/(1-xmg)+xo/(2*(1-xo)**2)-1) + osi*1873/temp*xsi**2*xo*(1/(1-xo)+1/(1-xsi)+xo/(2*(1-xo)**2)-1))
    lnKo = a2 + b2/temp + c2*pres/temp - lnGamma_O
    return lnKo

def equation6(x, a3, b3, c3, mgmg, hmg, omg, mgsi, oo, ho, osi):
    pres, temp, xh, xo, xmg, xsi = x
    lnGamma_O = -oo*1873/temp*np.log(1-xo) - (ho*1873/temp*xh*(1+np.log(1-xh)/xh-1/(1-xo)) + omg*1873/temp*xmg*(1+np.log(1-xmg)/xmg-1/(1-xo)) + osi*1873/temp*xsi*(1+np.log(1-xsi)/xsi-1/(1-xo))) + (ho*1873/temp*xh**2*xo*(1/(1-xo)+1/(1-xh)+xo/(2*(1-xo)**2)-1) + omg*1873/temp*xmg**2*xo*(1/(1-xo)+1/(1-xmg)+xo/(2*(1-xo)**2)-1) + osi*1873/temp*xsi**2*xo*(1/(1-xo)+1/(1-xsi)+xo/(2*(1-xo)**2)-1))
    lnKmg = a3 + b3/temp + c3*pres/temp - lnGamma_O - (-mgmg*1873/temp*np.log(1-xmg) - (hmg*1873/temp*xh*(1+np.log(1-xh)/xh-1/(1-xmg)) + omg*1873/temp*xo*(1+np.log(1-xo)/xo-1/(1-xmg)) + mgsi*1873/temp*xsi*(1+np.log(1-xsi)/xsi-1/(1-xmg))) + (hmg*1873/temp*xh**2*xmg*(1/(1-xmg)+1/(1-xh)+xmg/(2*(1-xmg)**2)-1) + omg*1873/temp*xo**2*xmg*(1/(1-xmg)+1/(1-xo)+xmg/(2*(1-xmg)**2)-1) + mgsi*1873/temp*xsi**2*xmg*(1/(1-xmg)+1/(1-xsi)+xmg/(2*(1-xmg)**2)-1)))    
    return lnKmg

def equation7(x, a4, b4, c4, sisi, hsi, osi, mgsi, oo, ho, omg):
    pres, temp, xh, xo, xmg, xsi = x
    lnKsi = a4 + b4/temp + c4*pres/temp - (-sisi*1873/temp*np.log(1-xsi) - (hsi*1873/temp*xh*(1+np.log(1-xh)/xh-1/(1-xsi)) + osi*1873/temp*xo*(1+np.log(1-xo)/xo-1/(1-xsi)) + mgsi*1873/temp*xmg*(1+np.log(1-xmg)/xmg-1/(1-xsi))) + (hsi*1873/temp*xh**2*xsi*(1/(1-xsi)+1/(1-xh)+xsi/(2*(1-xsi)**2)-1) + osi*1873/temp*xo**2*xsi*(1/(1-xsi)+1/(1-xo)+xsi/(2*(1-xsi)**2)-1) + mgsi*1873/temp*xmg**2*xsi*(1/(1-xsi)+1/(1-xmg)+xsi/(2*(1-xsi)**2)-1)))
    return lnKsi


# load data
#data_1 = pd.read_excel('OMgSiFe_exc&dis_1.xlsx')
#column_names = data_1.columns
#x_data_1 = np.array([
#    data_1['pres'], 
#    data_1['temp'], 
#    data_1['xo'], 
#    data_1['xmg'], 
#    data_1['xsi']
#])
#y_data_1 = np.array([
#    data_1['lnKo'], 
#    data_1['lnKmg'], 
#    data_1['lnKsi']
#])

data_2 = pd.read_excel('HOMgSiFe_exc&dis_1.xlsx')
column_names = data_2.columns
x_data_2 = np.array([
    data_2['pres'],
    data_2['temp'],
    data_2['xh'],
    data_2['xo'],
    data_2['xmg'],
    data_2['xsi']
])
y_data_2 = np.array([
    data_2['lnKh'],
    data_2['lnKo'],
    data_2['lnKmg'],
    data_2['lnKsi']
])

# build combined models
#def combined_model_1(x, a21, b21, c21, a31, b31, c31, a41, b41, c41, oo1, omg1, osi1, mgmg1, mgsi1, sisi1):
#    return [equation1(x, a21, b21, c21, oo1, omg1, osi1),
#            equation2(x, a31, b31, c31, mgmg1, omg1, mgsi1, oo1, osi1),
#            equation3(x, a41, b41, c41, sisi1, osi1, mgsi1, oo1, omg1)]

def combined_model_2(x, a1, b1, c1, a2, b2, c2, a3, b3, c3, a4, b4, c4, hh, ho, hmg, hsi, oo, omg, osi, mgmg, mgsi, sisi):
    return [equation4(x, a1, b1, c1, hh, ho, hmg, hsi, oo, omg, osi),
            equation5(x, a2, b2, c2, oo, ho, omg, osi),
            equation6(x, a3, b3, c3, mgmg, hmg, omg, mgsi, oo, ho, osi),
            equation7(x, a4, b4, c4, sisi, hsi, osi, mgsi, oo, ho, omg)]

#model_1 = Model(combined_model_1)
model_2 = Model(combined_model_2)

# build parameters
params = Parameters()
#params.add('a21', value=-2.0)
#params.add('b21', value=0.0)
#params.add('c21', value=66.0)
#params.add('a31', value=-6.0)
#params.add('b31', value=-7.0)
#params.add('c31', value=0.0)
#params.add('a41', value=0.0)
#params.add('b41', value=20.0)
#params.add('c41', value=-108.0)
#params.add('oo1', value=-5.0)
#params.add('omg1', value=-16.0)
#params.add('osi1', value=-4.0)
#params.add('mgmg1', value=0.0)
#params.add('mgsi1', value=0.0)
#params.add('sisi1', value=0.0)

#params.add('a1', value=0.0, vary=False)
params.add('a1', value=0.0)
params.add('b1', value=0.0, vary=False)
params.add('c1', value=0.0)
#params.add('a2', value=-2, vary=False)
params.add('a2', value=0.0, vary=False)
params.add('b2', value=0.0, vary=False)
#params.add('c2', value=66, vary=False)
params.add('c2', value=0.0, vary=False)
#params.add('a3', value=-6.0, vary=False)
params.add('a3', value=0.0)
#params.add('b3', value=-7000)
params.add('b3', value=0.0, vary=False)
params.add('c3', value=0.0)
params.add('a4', value=0.0)
#params.add('b4', value=-21000, vary=False)
params.add('b4', value=0.0, vary=False)
#params.add('c4', value=-108, vary=False)
params.add('c4', value=0.0)
params.add('hh', value=0.0)
params.add('ho', value=0.0)
params.add('hmg', value=0.0, vary=False)
params.add('hsi', value=0.0)
params.add('oo', value=0.0, vary=False)
#params.add('omg', value=-16.0, vary=False)
params.add('omg', value=0.0, vary=False)
params.add('osi', value=-4.0)
params.add('mgmg', value=0.0)
params.add('mgsi', value=0.0)
params.add('sisi', value=0.0)

#build shared parameters
#params['a21'].set(expr='a2')
#params['b21'].set(expr='b2')
#params['c21'].set(expr='c2')
#params['a31'].set(expr='a3')
#params['b31'].set(expr='b3')
#params['c31'].set(expr='c3')
#params['a41'].set(expr='a4')
#params['b41'].set(expr='b4')
#params['c41'].set(expr='c4')
#params['oo1'].set(expr='oo')
#params['omg1'].set(expr='omg')
#params['osi1'].set(expr='osi')
#params['mgmg1'].set(expr='mgmg')
#params['mgsi1'].set(expr='mgsi')
#params['sisi1'].set(expr='sisi')

#build joint residual to optimize all parameters in the two combined models simultaneously
#def joint_residual(params, x_data_1, y_data_1, x_data_2, y_data_2):
#    res_1 = model_1.eval(params=params, x=x_data_1) - y_data_1
#    res_2 = model_2.eval(params=params, x=x_data_2) - y_data_2
#    return np.concatenate([res_1, res_2], axis=None)

#result = minimize(joint_residual, params, args=(x_data_1, y_data_1, x_data_2, y_data_2))
result = model_2.fit(y_data_2, params, x=x_data_2)

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

# plot out the predicted data vs. input data
#y_pred_1 = model_1.eval(params=result.params, x=x_data_1)
y_pred_2 = model_2.eval(params=result.params, x=x_data_2)

plt.figure(figsize=(10, 5))
#plt.subplot(1, 2, 1)
#plt.scatter(y_data_1, y_pred_1)
#plt.plot([np.min(y_data_1), np.max(y_data_1)], [np.min(y_data_1), np.max(y_data_1)], 'k--', lw=2)
#plt.xlabel('Actual y_data_1')
#plt.ylabel('Predicted y_data_1')
#plt.title('Model 1')

#plt.subplot(1, 2, 2)
plt.scatter(y_data_2, y_pred_2)
plt.plot([np.min(y_data_2), np.max(y_data_2)], [np.min(y_data_2), np.max(y_data_2)], 'k--', lw=2)
plt.xlabel('Actual y_data_2')
plt.ylabel('Predicted y_data_2')
plt.title('Model 2')

plt.tight_layout()
plt.savefig('my_plot.png', dpi=600)
plt.show()

# save data
#y_data_1 = y_data_1.reshape(-1, 1)
#y_pred_1 = y_pred_1.reshape(-1, 1)
y_data_2 = y_data_2.reshape(-1, 1)
y_pred_2 = y_pred_2.reshape(-1, 1)
#data1_to_save = np.column_stack((y_data_1, y_pred_1))
data2_to_save = np.column_stack((y_data_2, y_pred_2))
#np.savetxt('y_data_1_and_pred_1.txt', data1_to_save)
np.savetxt('y_data_2_and_pred_2.txt', data2_to_save)
