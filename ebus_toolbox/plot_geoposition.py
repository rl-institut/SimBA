""" Module to implement plotting functionality of busses with georeferences"""
import datetime
import typing
import pathlib
from matplotlib.offsetbox import TextArea, AnnotationBbox

import ebus_toolbox.util as util
import pickle
import csv
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import ebus_toolbox.trip
    import ebus_toolbox.schedule
    import spiceev.scenario

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Mathematical function we need to plot
from matplotlib.colors import LinearSegmentedColormap
import matplotlib
import matplotlib.animation as animation

matplotlib.use("TkAgg")
import os

shared_axis = None
FIGSIZE = (14, 5)
LINEWIDTH = 2
DELIMITER = ';'
FRAMEALPHA = 0.8
ALPHA_POINTS = 0.1

###
rli_dblue = (0, 46 / 255, 80 / 255)
rli_green = (68 / 255, 175 / 255, 105 / 255)
rli_orange = (254 / 255, 127 / 255, 45 / 255)
rli_yellow = (241 / 255, 196 / 255, 15 / 255)
rli_lblue = (34 / 255, 116 / 255, 165 / 255)
rli_dgrey = (51 / 255, 88 / 255, 115 / 255)
### Not from official coperate design
rli_red = (192 / 255, 0, 0)
rli_green2 = (0, 176 / 255, 80 / 255)
rli_brown = (132 / 255, 60 / 255, 12 / 255)
rli_black = (0, 0, 0)
rli_colors = [rli_dblue, rli_green, rli_orange, rli_yellow, rli_lblue, rli_dgrey, rli_black,
              rli_red, rli_brown, rli_green2]
###

cdict = {'red': [[0.0, 0, 0],
                 [0.2, 34 / 255, 34 / 255],
                 [0.4, 68 / 255, 68 / 255],
                 [0.6, 254 / 255, 254 / 255],
                 [0.8, 241 / 255, 241 / 255],
                 [1.0, 255 / 255, 255 / 255]],
         'green': [[0.0, 46 / 255, 46 / 255],
                   [0.2, 116 / 255, 116 / 255],
                   [0.4, 175 / 255, 175 / 255],
                   [0.6, 127 / 255, 127 / 255],
                   [0.8, 196 / 255, 196 / 255],
                   [1.0, 255 / 255, 255 / 255]],
         'blue': [[0.0, 80 / 255, 80 / 255],
                  [0.2, 165 / 255, 165 / 255],
                  [0.4, 105 / 255, 105 / 255],
                  [0.6, 45 / 255, 45 / 255],
                  [0.8, 15 / 255, 15 / 255],
                  [1.0, 255 / 255, 255 / 255]]}
rli_rainbow_cmp = LinearSegmentedColormap('testCmap', segmentdata=cdict, N=256)

cdict = {'red': [[0.0, 0, 0],
                 [0.33, 0, 0],
                 [0.66, 34 / 255, 34 / 255],
                 [1, 255 / 255, 255 / 255]],
         'green': [[0.0, 0, 0],
                   [0.33, 46 / 255, 46 / 255],
                   [0.66, 116 / 255, 116 / 255],
                   [1, 255 / 255, 255 / 255]],
         'blue': [[0.0, 0, 0],
                  [0.33, 80 / 255, 80 / 255],
                  [0.66, 165 / 255, 165 / 255],
                  [1, 255 / 255, 45 / 255]]}
rli_blue_cmp = LinearSegmentedColormap('testCmap', segmentdata=cdict, N=256)

f_path=pathlib.Path(__file__).parent.parent
with open(f_path / "schedule_rebased_BVG_BFI.pickle", "rb") as file:
    schedule = pickle.load(file)

with open(f_path / "scenario_rebased_BVG_BFI.pickle", "rb") as file:
    scenario = pickle.load(file)

with open(f_path / "args_rebased_BVG_BFI.pickle", "rb") as file:
    args = pickle.load(file)

args.station_data_path = "C:/Users/paul.scheer/Python/bus_toolbox/eBus-Toolbox/Haltestellen.csv"

ANIMATION_DURATION_MIN = 480


def main():
    pickle_path = None #"vehicle_data_frames.pickle"
    with open(str(args.station_data_path), "r", encoding='utf-8') as f:
        delim = util.get_csv_delim(args.station_data_path)
        reader = csv.DictReader(f, delimiter=delim)
        stations = dict()
        for row in reader:
            elevation = float(row['elevation'])
            name = str(row['Endhaltestelle'])
            lon = float(row['lon'])
            lat = float(row['lat'])
            stations[name] = station(name, lat, lon, elevation)

    start_time = scenario.start_time
    time_step = datetime.timedelta(hours=1 / scenario.stepsPerHour)

    if not pickle_path:
        vehicle_data_frames = pd.DataFrame()
        for v_id, soc_data in scenario.vehicle_socs.items():
            rotations = get_sorted_rotations(v_id, schedule)
            trips = [trip for rot in rotations for trip in rot.trips]
            lats = []
            lons = []
            trip_nr = 0
            for counter, soc in enumerate(soc_data):
                try:
                    current_time = start_time + counter * time_step
                    current_trip, found_trip_index = find_current_trip(trips[trip_nr:],
                                                                       current_time)
                    trip_nr += found_trip_index
                    rel_time_of_trip = get_rel_time_of_trip(current_time, current_trip)
                    try:
                        current_station = stations[current_trip.departure_name]
                    except:
                        current_station = last_station
                    try:
                        next_station = stations[current_trip.arrival_name]
                    except:
                        next_station = current_station
                    last_station=current_station
                    lat, lon = current_station.get_lat_lon(next_station, rel_time_of_trip)
                    lats.append(lat)
                    lons.append(lon)
                except IndexError:
                    pass
            columns = [(v_id, "soc"), (v_id, "lat"), (v_id, "lon")]
            data = np.array([soc_data, lats, lons]).transpose()
            vehicle_data_frames = pd.concat(
                (vehicle_data_frames, pd.DataFrame(data, columns=columns)),
                axis=1)

        vehicle_data_frames.columns = pd.MultiIndex.from_tuples(vehicle_data_frames.columns,
                                                                names=['Vehicle_id', 'Data'])

        vehicle_data_frames.index = np.arange(start_time, start_time +
                                              len(vehicle_data_frames) * time_step, time_step)
    else:
        with open(pickle_path, "rb") as f:
            print("depickeling")
            vehicle_data_frames = pickle.load(f)

    stations_to_annotate = {name: stat for name, stat in stations.items() if
                            name in schedule.stations}

    vehicle_data_frames = vehicle_data_frames.iloc[:24*60,:]
    plot_merge_animate_battery(vehicle_data_frames, station_data=stations_to_annotate,
                               soc_threshold=1,track_black=False, vehicle_black=True,
                               save=False, repeat=True)


class station:
    def __init__(self, name, lat, lon, elevation):
        self.name = name
        self.lat = lat
        self.lon = lon
        self.elevation = elevation

    def get_lat_lon(self, next_station, rel_pos):
        lat_pos = (next_station.lat - self.lat) * rel_pos + self.lat
        lon_pos = (next_station.lon - self.lon) * rel_pos + self.lon
        return lat_pos, lon_pos


def get_sorted_rotations(v_id, schedule):
    rots = []
    for rot in schedule.rotations.values():
        if rot.vehicle_id == v_id:
            rots.append(rot)
    return sorted(rots, key=lambda x: x.departure_time)


def plot_merge_animate_battery(data: pd.DataFrame, station_data=None,
                               save=False, repeat=True, vehicle_black=False, track_black=True,
                               track_method="min", soc_threshold=1.0):

    data.xs("lat", level="Data", axis=1)
    lat_boundary = (data.xs("lat", level="Data", axis=1).min().min(),
                    data.xs("lat", level="Data", axis=1).max().max())
    lon_boundary = (data.xs("lon", level="Data", axis=1).min().min(),
                    data.xs("lon", level="Data", axis=1).max().max())


    v_max = 1
    v_min = 0
    fig = plt.figure(figsize=(14, 8))

    fig.suptitle('Vehicle Tracking', fontsize=14)

    sub2 = fig.add_subplot(111)
    ax = plt.gca()

    # Get all the vehicle ids from data
    vehicles = list(data.columns.levels[0])

    # Filter data for relevant SOCs
    for v_id in vehicles:
            if data[v_id]['soc'].min()>soc_threshold:
                data=data.drop(v_id,axis=1)
    data.columns = data.columns.remove_unused_levels()
    vehicles = list(data.columns.levels[0])

    # Plot track onto canvas via line which doesnt allow soc printing
    # or via scatter plot. With a scatter plot the points can take the color of the soc which drove
    # above the position.


    round_nr = 3
    stacked_array_unique = make_unique_data(vehicles, data, track_method)
    x_org=  stacked_array_unique[0, :].copy()
    y_org= stacked_array_unique[1, :].copy()
    #
    # stacked_array_unique[0, :] = stacked_array_unique[0, :] - min(stacked_array_unique[0, :])
    # stacked_array_unique[1, :] = stacked_array_unique[1, :] - min(stacked_array_unique[1, :])
    # stacked_array_unique[0, :] = stacked_array_unique[0, :]
    # stacked_array_unique[0, :] = stacked_array_unique[0, :] * 10 ** round_nr
    # stacked_array_unique[1, :] = stacked_array_unique[1, :] * 10 ** round_nr
    ##############
    # from scipy.ndimage import gaussian_filter as gauss
    # fig = plt.figure(figsize=(8, 8))
    #
    # fig.suptitle('Vehicle Tracking', fontsize=14)
    # sub2 = fig.add_subplot(111)
    # ax = plt.gca()
    #
    # stacked_array_pic = make_unique_data(vehicles, data, track_black, track_method,
    #                                      round_nr=round_nr)
    #
    # x = (stacked_array_pic[0, :] * 10 ** round_nr).astype(int)
    # X = x - np.min(x)
    # y = (stacked_array_pic[1, :] * 10 ** round_nr).astype(int)
    # Y = y - np.min(y)
    # z = stacked_array_pic[2, :]
    # idx1 = X
    # idx2 = Y
    # grid_data = z
    # grid = np.ones((np.max(Y) + 1, np.max(X) + 1))
    # grid[idx2, idx1] = grid_data
    # ax.imshow(gauss(grid, 10))
    # #############################
    #     fig = plt.figure(figsize=(14, 8))
    #
    #     fig.suptitle('Vehicle Tracking', fontsize=14)
    #
    #     sub2 = fig.add_subplot(111)
    #     ax = plt.gca()
    #
    #     x=stacked_array_unique[0,:]
    #     y=stacked_array_unique[1,:]
    #     z=stacked_array_unique[2,:]
    #
    #     steps=500
    #     X = np.linspace(stacked_array_unique[0, :].min(), stacked_array_unique[0, :].max(), steps)
    #     Y = np.linspace(stacked_array_unique[1, :].min(), stacked_array_unique[1, :].max(), steps)
    #     from scipy.interpolate import griddata
    #     Z = griddata((x, y), z,((X[None, :], Y[:, None])), method='linear')
    #     contour = sub2.contourf(X, Y, Z)
    # ########################
    #     ax.imshow(Z)
    # ########################

    if track_black:
        for v_id in vehicles:
            plot_dict = dict(color=(0.9, 0.9, 0.9), alpha=1)
            lns2_1 = sub2.plot(data[v_id]['lat'], data[v_id]['lon'],
                               **plot_dict, linestyle='-',
                               linewidth=2)
    else:
        plot_dict = dict(c=stacked_array_unique[2, :] * 1, alpha=0.8)
        lns2_1 = sub2.scatter(stacked_array_unique[0, :], stacked_array_unique[1, :],
                              **plot_dict, linestyle='-',
                              linewidth=0, vmin=v_min, vmax=v_max)
        lns2_1.set_clim(vmin=v_min, vmax=v_max)

    # Helper plot so we can plot a colorbar aftwards, not possible with only line plot
    plot_dict = dict()
    lns2_1 = sub2.scatter([], [],
                          **plot_dict, linestyle='-',
                          linewidth=0, vmin=v_min, vmax=v_max)

    cbar = fig.colorbar(lns2_1, ax=ax, pad=0.0)
    cbar.set_alpha(1)
    cbar.draw_all()
    cbar.set_label("SOC")

    # Plot Station Names
    if station_data:
        for station in station_data.values():
            # xy=((station.lat-np.min(x_org))*10**round_nr, (station.lon-np.min(y_org))*10**round_nr)
            xy=(station.lat, station.lon)
            ax.annotate(station.name,
                        xy=xy, xycoords='data', xytext=(0, 5), textcoords='offset points', fontsize=8, ha='center')
            ax.plot(xy[0], xy[1], 'ko', zorder=100)

    artist_objects = []
    for _ in vehicles:
        lns2_1 = sub2.scatter([0, 0], [0, 0], [0, 100])
        artist_objects.append(lns2_1)

    d_bound = 0.005
    sub2.set_ylim((lon_boundary[0] - d_bound, lon_boundary[1] + d_bound))
    sub2.set_xlim((lat_boundary[0] - d_bound, lat_boundary[1] + d_bound))
    sub2.set_ylabel('Position')
    sub2.set_xlabel('Position')

    text_box = TextArea(data.index[0])

    time_annotation_box = [AnnotationBbox(text_box, xy=(0.1, 0.1), xycoords='axes fraction',
                                          fontsize=15)]
    ax.add_artist(time_annotation_box[0])
    hover_texts = [str(v_id) for v_id in vehicles]

    def animate(i, data, counter=[]):
        roll_v = 1
        len_v = 1
        ax = plt.gca()
        vehicles = list(data.columns.levels[0])
        if not counter:
            counter.append(i)
        else:
            counter[0] += 1
        i = counter[0]
        start_index = i * roll_v + 8 * 60
        end_index = start_index + len_v
        try:
            ax.artists.remove(time_annotation_box[0])
        except:
            return artist_objects
        text_box = TextArea(data.index[start_index])
        time_annotation_box[0] = AnnotationBbox(text_box, xy=(0.1, 0.1), xycoords='axes fraction',
                                                fontsize=15)
        ax.add_artist(time_annotation_box[0])
        #
        #

        plot_dict = dict(color=(0, 0, 0))
        artist_objects[0].remove()
        lat_data=data.xs("lat", level="Data", axis=1).iloc[start_index:end_index].values
        lon_data=data.xs("lon", level="Data", axis=1).iloc[start_index:end_index].values
        soc_data=data.xs("soc", level="Data", axis=1).iloc[start_index:end_index].mean().values
        if vehicle_black:
            plot_dict = dict(color=(0, 0, 0))
        else:
            plot_dict = dict(c=soc_data)
        artist_objects[0] = sub2.scatter(lat_data,
                                         lon_data,
                                         **plot_dict,
                                         vmin=0, vmax=1,
                                         linestyle='-',
                                         linewidth=0, zorder=50)

        #
        #
        # for k, v_id in enumerate(vehicles):
        #     artist_objects[k].remove()
        #     mean_soc = data[v_id]["soc"][start_index:end_index].mean()
        #     if vehicle_black:
        #         plot_dict = dict(color=(0, 0, 0))
        #     else:
        #         l = len(data[v_id]["lat"][start_index:end_index])
        #         plot_dict = dict(c=mean_soc * np.ones(l))
        #     artist_objects[k] = sub2.scatter(data[v_id]["lat"][start_index:end_index],
        #                                      data[v_id]["lon"][start_index:end_index],
        #                                      **plot_dict,
        #                                      vmin=0, vmax=1,
        #                                      linestyle='-',
        #                                      linewidth=0, zorder=50)
        #     hover_texts[k] = (v_id, str(round(mean_soc, 3)))

        if end_index >= len(data):
            counter[0] = 0

        return artist_objects

    ax.axis('off')
    # create animation using the animate() function

    my_animation = animation.FuncAnimation(fig,
                                          lambda i: animate(i, data),
                                          frames=ANIMATION_DURATION_MIN,
                                          interval=15, blit=False, repeat=repeat)

    #
    sub2.set_ylabel('Position')
    sub2.set_xlabel('Position')

    fig.canvas.mpl_connect("motion_notify_event", lambda event: hover_for_scatter(
        event, fig, ax, artist_objects, hover_texts))
    # Toggle save  for saving
    if save:
        my_animation.save('SOC_animation.gif', writer='imagemagick')

    fig.tight_layout()

    def on_close(event):
        plt.close()
        print('Closed Figure!')
    fig.canvas.mpl_connect('close_event',lambda event:on_close(event))
    plt.show(block=True)
    print("stopped")

def make_unique_data(vehicles, data, track_method, round_nr=None):
    v_iter = iter(vehicles)
    v_id = next(v_iter)
    data = pd.DataFrame(data)
    if not round_nr:
        r = lambda x: x
    else:
        r = lambda i: round(i, round_nr)
    stacked_array = np.array((r(data[v_id]['lat']), r(data[v_id]['lon']), data[v_id]['soc']))
    for v_id in v_iter:
        stacked_array = np.hstack(
            (stacked_array, (r(data[v_id]['lat']), r(data[v_id]['lon']), data[v_id]['soc'])))
    stacked_array = stacked_array.astype(float)
    # stacked_array has all vehicle socs and geo positions in
    unique_mask = np.unique(stacked_array[0:2, :], axis=1, return_index=True)[1]
    stacked_array_unique = stacked_array[:, unique_mask]

    # find soc data for stacked array depending on method. Since its only needed if the track is not black check this as well

    apply_function = None
    if track_method == "mean":
        apply_function = np.mean
    elif track_method == "max":
        apply_function = np.max
    else:
        apply_function = np.min

    for geo_index in range(len(stacked_array_unique)):
        geo_loc = stacked_array_unique[0:2, geo_index]
        found_positions = np.all([stacked_array[0:2, :].T == geo_loc], axis=0)[:, 0]
        fill_value = apply_function(stacked_array[2, found_positions])
        stacked_array_unique[2, geo_index] = fill_value

        return stacked_array_unique


def hover_for_scatter(event, fig, ax, plot_points: typing.Iterable[matplotlib.lines.Line2D],
                      hover_texts: typing.Iterable[str], annotations=[]):
    """Called when user hovers over plot.
    Checks if user hovers over point. If so, delete old annotation and
    create new one with relevant info from the hover_texts list.
    If user does not hover over point, remove annotation, if any.
    """
    if len(annotations) == 0:
        annotations.append(None)
    for i, points in enumerate(plot_points):
        if points and event.inaxes == ax:
            # results shown, mouse within plot: get event info
            # cont: any points hovered?
            # ind:  list of points hovered
            cont, ind = points.contains(event)

            if cont and "ind" in ind:
                ind = ind["ind"]
                # points hovered
                # get all point coordinates
                xy = points.get_offsets().data
                text = hover_texts[i]

                # # remove old annotation
                if annotations and annotations[0]:
                    annotations[0].remove()
                    annotations[0] = None

                # create new annotation
                annotations[0] = ax.annotate(
                    text,
                    xy=(xy[ind[0]][0], xy[ind[0]][1]),
                    xytext=(-20, 20),
                    textcoords="offset points",
                    bbox=dict(boxstyle="round", fc="w"),
                    arrowprops={'arrowstyle': "-"},
                    annotation_clip=False)
                fig.canvas.draw()


def find_current_trip(trips, current_time):
    for i, trip in enumerate(trips):
        try:
            next_trip = trips[i + 1]
            if next_trip.departure_time > current_time:
                return trip, i
        except IndexError:
            return trip, i
    else:
        raise Exception("no trip found for time: " + str(current_time))


def get_rel_time_of_trip(current_time: datetime.datetime, current_trip: 'ebus_toolbox.trip.Trip'):
    trip_time= (current_trip.arrival_time - current_trip.departure_time)
    rel_time = (current_time - current_trip.departure_time) / \
               max(trip_time,datetime.timedelta(minutes=1))
    return min(1, max(rel_time, 0))


if __name__ == "__main__":
    main()
