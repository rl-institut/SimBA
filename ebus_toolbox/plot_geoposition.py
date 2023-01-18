""" Module to implement plotting functionality of busses with georeferences"""
import datetime

import ebus_toolbox.util as util
import pickle
import csv
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import ebus_toolbox.trip
    import ebus_toolbox.schedule
    import spiceev.scenario

import requests
import json
from time import sleep
import io

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import math

# Mathematical function we need to plot
from matplotlib.colors import LinearSegmentedColormap
import matplotlib
import matplotlib.animation as animation
matplotlib.use("TkAgg")
import os
shared_axis=None
FIGSIZE=(14,5)
LINEWIDTH=2
DELIMITER = ';'
FRAMEALPHA=0.8
ALPHA_POINTS=0.1


###
rli_dblue =(0, 46 / 255, 80 / 255)
rli_green=(68 / 255, 175 / 255, 105 / 255)
rli_orange=(254 / 255, 127 / 255, 45 / 255)
rli_yellow=(241 / 255, 196 / 255, 15 / 255)
rli_lblue=(34/ 255, 116 / 255, 165 / 255)
rli_dgrey=(51 / 255, 88 / 255, 115 / 255)
### Not from official coperate design
rli_red =(192/255,0,0)
rli_green2 =(0, 176/255,80/255)
rli_brown = (132 / 255, 60 / 255, 12 / 255)
rli_black=(0,0,0)
rli_colors = [rli_dblue, rli_green, rli_orange, rli_yellow,rli_lblue, rli_dgrey,rli_black,rli_red, rli_brown,rli_green2]
###

cdict = {'red':   [[0.0,  0, 0],
                   [0.2,  34/255, 34/255],
                   [0.4,  68/255, 68/255],
                   [0.6,  254/255, 254/255],
                   [0.8,  241/255, 241/255],
                   [1.0,  255/255, 255/255]],
         'green': [[0.0,  46/255, 46/255],
                   [0.2,  116/255, 116/255],
                   [0.4,  175/255, 175/255],
                   [0.6,  127/255, 127/255],
                   [0.8,  196/255, 196/255],
                   [1.0,  255/255, 255/255]],
         'blue':  [[0.0,  80/255, 80/255],
                   [0.2,  165/255, 165/255],
                   [0.4, 105/255, 105/255],
                   [0.6, 45/255, 45/255],
                   [0.8, 15/255, 15/255],
                   [1.0, 255/255, 255/255]]}
rli_rainbow_cmp = LinearSegmentedColormap('testCmap', segmentdata=cdict, N=256)

cdict = {'red':   [[0.0,  0, 0],
                   [0.33,  0, 0],
                   [0.66,  34/255, 34/255],
                   [1,  255/255, 255/255]],
         'green': [[0.0,  0, 0],
                   [0.33,  46/255, 46/255],
                   [0.66,  116/255, 116/255],
                   [1,  255/255, 255/255]],
         'blue':  [[0.0,  0, 0],
                   [0.33,  80/255, 80/255],
                   [0.66,  165/255, 165/255],
                   [1, 255/255, 45/255]]}
rli_blue_cmp = LinearSegmentedColormap('testCmap', segmentdata=cdict, N=256)

with open("schedule_opt.pickle", "rb") as file:
    schedule = pickle.load(file)

with open("scenario_opt.pickle", "rb") as file:
    scenario = pickle.load(file)

with open("args_buffered_all_depb.pickle", "rb") as file:
    args = pickle.load(file)

args.station_data_path= "C:/Users/paul.scheer/Python/bus_toolbox/eBus-Toolbox/data/buffered_all_stations.csv"

class station:
    def __init__(self, name, lat, lon, elevation):
        self.name = name
        self.lat = lat
        self.lon = lon
        self.elevation = elevation

    def get_lat_lon(self, next_station, rel_pos):
        lat_pos= (next_station.lat-self.lat)*rel_pos+self.lat
        lon_pos = (next_station.lon - self.lon) * rel_pos + self.lon
        return lat_pos, lon_pos


with open(str(args.station_data_path), "r", encoding='utf-8') as f:
    delim = util.get_csv_delim(args.station_data_path)
    reader = csv.DictReader(f, delimiter=delim)
    stations = dict()
    for row in reader:
        elevation = float(row['elevation'])
        name = str(row['Endhaltestelle'])
        lon = float(row['lon'])
        lat = float(row['lat'])
        stations[name]=station(name, lat, lon, elevation)

def get_sorted_rotations(v_id, schedule):
    rots=[]
    for rot in schedule.rotations.values():
        if rot.vehicle_id==v_id:
            rots.append(rot)
    return sorted(rots, key= lambda x: x.departure_time)


def plot_merge_animate_battery(data: pd.DataFrame, z_ax=None,
                               save=False, repeat=True, vehicle_black=False, track_black=True):
    v_max=1
    v_min=0
    fig = plt.figure(figsize=(14, 8))

    fig.suptitle('Vehicle Tracking', fontsize=14)

    sub2 = fig.add_subplot(111)
    ax = plt.gca()
    data.xs("lat", level="Data", axis=1)
    lat_boundary = (data.xs("lat", level="Data", axis=1).min().min(),data.xs("lat", level="Data", axis=1).max().max())
    lon_boundary = (data.xs("lon", level="Data", axis=1).min().min(),data.xs("lon", level="Data", axis=1).max().max())

    # Get all the vehicle ids from data
    vehicles =    list(data.columns.levels[0])
    data_rounded =    data.copy()

    round_decimals = 3    #4 decimals approx 110m in germany
    all_lats=[]
    all_lons=[]
    # for v_id in vehicles:
    #     data_rounded.loc[:, (v_id, "lat")] = data_rounded[v_id]["lat"].round(round_decimals)
    #     data_rounded.loc[:, (v_id, "lon")]= data_rounded[v_id]["lon"].round(round_decimals)
    #     all_lats.extend(data_rounded[v_id]["lat"])
    #     all_lons.extend(data_rounded[v_id]["lon"])
    #
    # all_geos=np.array((all_lats, all_lons))
    # all_geos_unique=np.unique(all_geos, axis=1)
    # # Plot track on the canvas
    # if track_black:
    #     plot_dict = dict(color=(0, 0, 0), alpha=0.1)
    # else:
    #     plot_dict = dict(c=data[v_id]['soc'] * 100, alpha=0.01)
    #
    # lns2_1 = sub2.scatter(all_geos_unique[0,:],all_geos_unique[1,:],
    #                       **plot_dict, linestyle='-',
    #                       linewidth=0, vmin=v_min, vmax=v_max)
    # lns2_1.set_clim(vmin=v_min, vmax=v_max)


    for v_id in vehicles:
        if track_black:
            plot_dict = dict(color=(0.9, 0.9, 0.9), alpha=1)
            # lns2_1 = sub2.plot(data[v_id]['lat'], data[v_id]['lon'],
            #                    **plot_dict, linestyle='-',
            #                    linewidth=2)
            lns2_1 = sub2.scatter(data[v_id]['lat'], data[v_id]['lon'],
                               **plot_dict, linestyle='-',
                               linewidth=2, s=0, vmin=v_min, vmax=v_max)
        else:
            plot_dict = dict(c=data[v_id]['soc'] * 100, alpha=0.01)
            lns2_1 = sub2.scatter(data[v_id]['lat'], data[v_id]['lon'],
                               **plot_dict, linestyle='-',
                               linewidth=0, vmin=v_min, vmax=v_max)
            lns2_1.set_clim(vmin=v_min, vmax=v_max)

    # Helper plot so we can plot a colorbar aftwards, not possible with only line plot
    # lns2_1 = sub2.scatter([], [],
    #                       **plot_dict, linestyle='-',
    #                       linewidth=0, vmin=v_min, vmax=v_max)


    cbar = fig.colorbar(lns2_1, ax=ax, pad=0.0)
    cbar.set_alpha(1)
    cbar.draw_all()
    cbar.set_label("SOC")

    sub2.set_ylabel('Position')
    sub2.set_xlabel('Position')

    artist_objects=[]
    for _ in vehicles:
        lns2_1 = sub2.scatter([0, 0], [0, 0], [0, 100])
        artist_objects.append(lns2_1)

    d_bound=0.001
    sub2.set_ylim((lon_boundary[0]-d_bound, lon_boundary[1]+d_bound))
    sub2.set_xlim((lat_boundary[0]-d_bound, lat_boundary[1]+d_bound))

    def animate(i, data , counter=[]):
        roll_v = 3
        len_v = 5
        vehicles = list(data.columns.levels[0])
        if not counter:
            counter.append(i)
        else:
            counter[0] +=1
        i =counter[0]
        start_index=i*roll_v
        end_index=start_index+len_v

        for k, v_id in enumerate(vehicles):
            artist_objects[k].remove()
            if vehicle_black:
                plot_dict = dict(color=(0, 0, 0))
            else:
                mean_soc = data[v_id]["soc"][start_index:end_index].mean()
                l = len(data[v_id]["lat"][start_index:end_index])
                plot_dict = dict(c=mean_soc * np.ones(l))
            artist_objects[k] = sub2.scatter(data[v_id]["lat"][start_index:end_index],
                                             data[v_id]["lon"][start_index:end_index],
                                             **plot_dict,
                                             vmin=0, vmax=1,
                                             linestyle='-',
                                             linewidth=0)
        if end_index >= len(data[v_id]["soc"]):
            counter[0]=0

        return artist_objects

    ax.axis('off')
    # create animation using the animate() function
    myAnimation = animation.FuncAnimation(fig,
                                          lambda i: animate(i, data),
                                          frames=10,
                                          interval=50, blit=False, repeat=repeat)
    #
    sub2.set_ylabel('Position')
    sub2.set_xlabel('Position')

    # y_d = (y_max - y_min) * 0.02
    # x_d = (x_max - x_min) * 0.02
    # sub2.set_ylim((y_min - y_d, y_max + y_d))
    # sub2.set_xlim((x_min - x_d, x_max + x_d))

    # Toggle Comment of for saving
    if save:
        myAnimation.save('SOC_Animation_all_OLIA_w_locations.gif', writer='imagemagick')

    fig.tight_layout()

    plt.show()

def find_current_trip(trips, current_time):
    for i, trip in enumerate(trips):
        try:
            next_trip = trips[i + 1]
            if next_trip.departure_time>current_time:
                return trip, i
        except IndexError:
            return trip, i
    else:
        raise Exception("no trip found for time: " + str(current_time))

def get_rel_time_of_trip(current_time: datetime.datetime, current_trip: 'ebus_toolbox.trip.Trip'):
    rel_time=(current_time-current_trip.departure_time)/\
             (current_trip.arrival_time-current_trip.departure_time)
    return min(1,max(rel_time,0))

start_time=scenario.start_time
time_step = datetime.timedelta(hours=1/scenario.stepsPerHour)

vehicle_data_frames=pd.DataFrame()
c = 0
for v_id, soc_data in scenario.vehicle_socs.items():
    rotations=get_sorted_rotations(v_id, schedule)
    trips = [trip for rot in rotations for trip in rot.trips ]
    lats = []
    lons = []
    trip_nr=0
    for counter, soc in enumerate(soc_data):
        try:
            current_time=start_time+counter*time_step
            current_trip, found_trip_index = find_current_trip(trips[trip_nr:],current_time)
            trip_nr += found_trip_index
            rel_time_of_trip=get_rel_time_of_trip(current_time, current_trip)
            current_station = stations[current_trip.departure_name]
            next_station = stations[current_trip.arrival_name]
            lat, lon = current_station.get_lat_lon(next_station,rel_time_of_trip)
            lats.append(lat)
            lons.append(lon)
        except IndexError:
            pass
    columns=[(v_id, "soc"), (v_id, "lat") , (v_id, "lon")]
    data = np.array([soc_data, lats, lons]).transpose()
    vehicle_data_frames=pd.concat((vehicle_data_frames, pd.DataFrame(data ,columns=columns)),axis=1)
    c +=1
    if c>4:
        break


vehicle_data_frames.columns = pd.MultiIndex.from_tuples(vehicle_data_frames.columns, names=['Vehicle_id','Data'])
plot_merge_animate_battery(vehicle_data_frames)
