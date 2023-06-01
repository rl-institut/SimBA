""" Sensitivity analysis based on monte carlo simulation
"""
import scipy
from scipy.stats import gamma
import pandas as pd
import numpy as np
from random import *
import random
import secrets

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

    range_hour = ['0-1', '1-2', '2-3', '3-4', '4-5', '5-6', '6-7', '7-8', '8-9', '9-10', '10-11', '11-12',
                  '12-13', '13-14', '14-15', '15-16', '16-17', '17-18', '18-19', '19-20', '20-21', '21-22',
                  '22-23', '23-0']

    delay['hour'] = range_hour

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
    day['Date'] = day['Date'].dt.hour
    day.rename(columns={'Date': 'hour', 'TT_TU': 'temperature'}, inplace=True)
    day = day.reset_index(drop=True)

    day['hour'] = pd.date_range(start='1/1/2023 00:00:00', end='1/1/2023 23:00:00', periods=24)
    day['hour'] = pd.to_datetime(day['hour'], format).apply(lambda x: x.time())

    return day

# 3. reduced power

def get_reduced_power():
    """Gets reduced charging power for different bus types
    :return: reduced_power_opps, reduced_power_depb, reduced_power_oppb
    """
    power_opps = np.arange(start=50, stop=400, step=50)
    reduced_power_opps = secrets.choice(power_opps)

    power_deps_depb = np.arange(start=30, stop=90, step=30)
    reduced_power_depb = secrets.choice(power_deps_depb)

    power_deps_oppb = np.arange(start=30, stop=120, step=30)
    reduced_power_oppb = secrets.choice(power_deps_oppb)

    return reduced_power_opps, reduced_power_depb, reduced_power_oppb


# 4. default hpc


def get_default_hpc():
    """
    Gets default hpc by choosing from a preset of electrified_stations.json
    :return: default_hpc
    """
    list_json = ["data/bvg/optimized_stations.json",
                 "data/bvg/optimized_stations_1.json",
                 "data/bvg/optimized_stations_2.json",
                 "data/bvg/optimized_stations_3.json",
                 "data/bvg/optimized_stations_4.json"]

    default_hpc = secrets.choice(list_json)

    return default_hpc

# 5. battery aging


def get_battery_aging():
    """
    Gets battery aging by choosing from a preset of vehicle_types.json
    :return: battery_aging
    """
    list_json = ["data/bvg/vehicle_types_95.json",
                 "data/bvg/vehicle_types_90.json",
                 "data/bvg/vehicle_types_85.json",
                 "data/bvg/vehicle_types_80.json",
                 "data/bvg/vehicle_types.json"]

    battery_aging = secrets.choice(list_json)

    return battery_aging

# 6. depot delay


def get_depot_delay():
    """
    Gets depot delay
    :return:
    """

    depot_delay = pd.read_csv('data/monte_carlo/depot_delay_data.csv', sep=',', decimal='.')

    data = pd.DataFrame()
    data['hour'] = depot_delay['hour'].unique()
    data = data.sort_values(by='hour')
    data = data.reset_index(drop=True)
    means = depot_delay.groupby('hour')['standing_time'].mean()
    means = means.reset_index(drop=True)
    stds = depot_delay.groupby('hour')['standing_time'].std()
    stds = stds.reset_index(drop=True)
    data['mean'] = means
    data['std'] = stds

    list_s = []

    for i in range(data.shape[0]):
        mu = data.loc[i, 'mean']
        sigma = data.loc[i, 'std']
        s = np.random.normal(mu, sigma, 1)
        list_s.append(s[0])

    data['s'] = list_s
    data.s = data.s.round()
    data.drop(['mean', 'std'], axis=1, inplace=True)

    delay_dep = dict(zip(data.hour, data.s))
    delay_dep['else'] = 0

    return delay_dep
