""" Sensitivity analysis based on monte carlo simulation
"""
import scipy
from scipy.stats import gamma
from ebus_toolbox.util import read_arguments
import pandas as pd
import numpy as np

# define sturgeon values

# 1. buffer times

# reading data
delay_data = pd.read_csv('data/monte_carlo/Fitted_Delay_at_terminus_over_weekdays.csv', sep=';', decimal='.')

# generate buffer times from gamma distribution
list_r = []

for i in range(delay_data.shape[0]):
    gamma_a = delay_data.loc[i, 'gamma_a']
    if delay_data.loc[i, 'gamma_scale'] < 0:
        gamma_scale = delay_data.loc[i, 'gamma_scale']*0
    else:
        gamma_scale = delay_data.loc[i, 'gamma_scale']
    r = scipy.stats.gamma.rvs(gamma_a, loc=0, scale=gamma_scale, size=1)
    list_r.append(r[0])

delay_data['r'] = list_r

# averaging of values over hour of the day
delay = delay_data.groupby('hour')['r'].mean()

# create dictionary for cfg file 

delay = delay.to_dict()

# 2. temperature

# reading data
weather_data = pd.read_csv('data/monte_carlo/produkt_tu_stunde_19510101_20211231_00433.csv', sep=';', decimal='.')

weather_data['MESS_DATUM'] = weather_data['MESS_DATUM'].astype('str')

for i in range(weather_data.shape[0]):
    weather_data.loc[i, 'year'] = weather_data.loc[i, 'MESS_DATUM'][:4]
    weather_data.loc[i, 'month'] = weather_data.loc[i, 'MESS_DATUM'][4:6]
    weather_data.loc[i, 'day'] = weather_data.loc[i, 'MESS_DATUM'][6:8]
    weather_data.loc[i, 'hour'] = weather_data.loc[i, 'MESS_DATUM'][8:10]

weather_data['Date'] = pd.to_datetime(weather_data[['year', 'month', 'day', 'hour']])
weather_data = weather_data[weather_data['year'] == '2021']
weather = pd.DataFrame(weather_data.groupby('hour')['TT_TU'].mean())
weather = weather.to_dict()

# generate random values from normal distribution
mean_2 = weather_data.mean()

mu = mean_2.TT_TU  # mean
sigma = weather_data['TT_TU'].std()  # standard deviation
s = np.random.normal(mu, sigma, 1)

# configs erstellen

args = read_arguments()

# config manipulieren/ergänzen

args.buffer_time = delay

# gc_power_opps, gc_power_deps (Netzauslastung)
# cs_power_opps, cs_power_deps_depb, cs_power_deps_oppb (Technik)

# ergänzung: temp, verspätung, aufall hpc, mehrkilometer, (batteriealterung), zeitl. abläufe depot

# neue config speichern
