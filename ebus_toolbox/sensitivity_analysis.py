""" Sensitivity analysis based on monte carlo simulation
"""
import scipy
from scipy.stats import gamma
import pandas as pd
import numpy as np
from random import *
import random

# define sturgeon values

# 1. buffer times (delay on track)


def get_buffer_times():
    """Gets buffer times based on gamma function

    :return: delay
    """
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
    delay = round(delay_data.groupby('hour')['r'].mean())

    # preparation for dictionary

    delay = delay.reset_index()
    delay['hour'] = delay['hour'].astype(int)
    delay['hour'] = delay['hour'].astype(str)
    delay['r'] = delay['r'].astype(int)

    # create dictionary for cfg file

    delay = dict(zip(delay.hour, delay.r))
    delay['else'] = 0

    return delay

# 2. temperature


def get_temperature():
    """Gets temperature profile of a random day

    :return: day
    """
    # reading data
    weather_data = pd.read_csv('data/monte_carlo/produkt_tu_stunde_19510101_20211231_00433.csv', sep=';', decimal='.')

    weather_data['MESS_DATUM'] = weather_data['MESS_DATUM'].astype('str')

    # get year, month, day and hour from data

    for i in range(weather_data.shape[0]):
        weather_data.loc[i, 'year'] = weather_data.loc[i, 'MESS_DATUM'][:4]
        weather_data.loc[i, 'month'] = weather_data.loc[i, 'MESS_DATUM'][4:6]
        weather_data.loc[i, 'day'] = weather_data.loc[i, 'MESS_DATUM'][6:8]
        weather_data.loc[i, 'hour'] = weather_data.loc[i, 'MESS_DATUM'][8:10]

    weather_data['Date'] = pd.to_datetime(weather_data[['year', 'month', 'day', 'hour']])
    weather_data = weather_data.astype({'year': 'int', 'month': 'int', 'day': 'int'})

    # select random day by year, month and day

    year = randint(2010, 2021)
    month = randint(1, 12)
    if month == 2:
        day = randint(1, 28)
    elif month == 4 or 6 or 9 or 11:
        day = randint(1, 30)
    else:
        day = randint(1, 31)

    day = weather_data[(weather_data.year == year) & (weather_data.month == month) & (weather_data.day == day)]
    day = day[['Date', 'TT_TU']]
    day = day.reset_index(drop=True)
    day = day.drop(['Date'], axis=1)

    day = day.to_dict()

    return day

# gc_power_opps, gc_power_deps (Netzauslastung)
# cs_power_opps, cs_power_deps_depb, cs_power_deps_oppb (Technik)

# cs_power_opps

def get_reduced_power():
    "creates list of different values which are smaller then the default power "
    power_opps = np.arange(start=50, stop=350, step=50)
    reduced_power_opps = random.choice(power_opps)


    # cs_power_deps_depb
    power_deps_depb = np.arange(start=30, stop=90, step=30)
    reduced_power_depb = random.choice(power_deps_depb)


    # cs_power_deps_depb
    power_deps_oppb = np.arange(start=30, stop=120, step=30)
    reduced_power_oppb = random.choice(power_deps_depb)

    return reduced_power_opps, reduced_power_depb, reduced_power_oppb

# 4. network utilization
# gc_power_opps, gc_power_deps (Netzauslastung)

def get_grid_utilization():
    """
    Gets grid utilization
    :return:
    """
    return

# 5. extra mileage


def get_extra_mileage():
    """
    Gets extra mileage
    :return:
    """
    return

# 6. default hpc


def get_default_hpc():
    """
    Gets default hpc
    :return:
    """
    return

# 7. battery aging


def get_battery_aging():
    """
    Gets battery aging
    :return:
    """
    return

# 8. depot delay


def get_depot_delay():
    """
    Gets depot delay
    :return:
    """
    return
