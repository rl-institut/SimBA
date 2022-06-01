""" Sensitivity analysis based on monte carlo simulation
"""
import scipy
from scipy.stats import gamma
from ebus_toolbox.util import read_arguments
import pandas as pd

# define sturgeon values

# 1. buffer times

# reading data
delay_data = pd.read_csv('data/monte_carlo/Fitted_Delay_at_terminus_over_weekdays.csv', sep=';', decimal='.')

# generate buffer times from gamma distribution
list_r = []

for i in range(delay_data.shape[0]):
    gamma_a = delay_data.loc[i, 'gamma_a']
    if delay_data.loc[i, 'gamma_scale'] < 0:
        gamma_scale = delay_data.loc[i, 'gamma_scale']*-1
    else:
        gamma_scale = delay_data.loc[i, 'gamma_scale']
    r = scipy.stats.gamma.rvs(gamma_a, loc=0, scale=1/gamma_scale, size=1)
    list_r.append(r[0])

delay_data['r'] = list_r

# averaging of values over hour of the day
means = delay_data.groupby('hour').mean()

# create dictionary for cfg file 
delay = {'0': means.r[0],
         '1': means.r[1],
         '2': means.r[2],
         '3': means.r[3],
         '4': means.r[4],
         '5': means.r[5],
         '6': means.r[6],
         '7': means.r[7],
         '8': means.r[8],
         '9': means.r[9],
         '10': means.r[10],
         '11': means.r[11],
         '12': means.r[12],
         '13': means.r[13],
         '14': means.r[14],
         '15': means.r[15],
         '16': means.r[16],
         '17': means.r[17],
         '18': means.r[18],
         '19': means.r[19],
         '20': means.r[20],
         '21': means.r[21],
         '22': means.r[22],
         '23': means.r[23],
         'else': 2
         }

# configs erstellen


def read_cfg():
    args = read_arguments()
    return "converged"


if __name__ == '__main__':
    read_cfg()

# config manipulieren/ergänzen

args.buffer_time = delay

# gc_power_opps, gc_power_deps (Netzauslastung)
# cs_power_opps, cs_power_deps_depb, cs_power_deps_oppb (Technik)

# ergänzung: temp, verspätung, aufall hpc, mehrkilometer, (batteriealterung), zeitl. abläufe depot

# neue config speichern
