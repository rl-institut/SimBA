import json
import warnings
import subprocess


def get_git_revision_hash() -> str:
    return subprocess.check_output(['git', 'rev-parse', 'HEAD']).decode('ascii').strip()


def save_version(file_path):
    with open(file_path, "w", encoding='utf-8') as f:
        f.write("Git Hash eBus-Toolbox:" + get_git_revision_hash())


def get_buffer_time(trip, default=0):
    """ Get buffer time at arrival station of a trip.
        Buffer_time is an abstraction of delays like docking procedures and
        is added to the planned arrival time

    :param trip: The of buffer time of this trips arrival is returned.
    :type trip: ebus_toolbox.Trip
    :param default: Default buffer time if no station specific buffer time is given. [minutes]
    :type default: dict, numeric
    :return: Buffer time
    :rtype: numeric

    NOTE: Buffertime dictionaries map hours of the day to a buffer time.
    Keys are ranges of hours and corresponding values provide buffer time in
    minutes for that time range.
    An entry with key "else" is a must if not all
    hours of the day are covered.
    E.g.
        buffer_time = {
            "10-22": 2,
            "22-6": 3,
            "else": 1
        }
    """
    schedule = trip.rotation.schedule
    buffer_time = schedule.stations.get(trip.arrival_name, {}).get('buffer_time', default)

    # distinct buffer times depending on time of day can be provided
    # in that case buffer time is of type dict instead of int
    if isinstance(buffer_time, dict):
        # sort dict to make sure 'else' key is last key
        buffer_time = {key: buffer_time[key] for key in sorted(buffer_time)}
        current_hour = trip.arrival_time.hour
        for time_range, buffer in buffer_time.items():
            if time_range == 'else':
                buffer_time = buffer
                break
            else:
                start_hour, end_hour = [int(t) for t in time_range.split('-')]
                if end_hour < start_hour:
                    if current_hour >= start_hour or current_hour < end_hour:
                        buffer_time = buffer
                        break
                else:
                    if start_hour <= current_hour < end_hour:
                        buffer_time = buffer
                        break

    return buffer_time


def uncomment_json_file(f, char='//'):
    """
    Remove comments from JSON file.

    Wouldn't it be nice to have comments in JSON data? Now you can!
    Both full-line comments and trailing lines are supported.
    :param f: file to read
    :type f: JSON file handle
    :param char: char sequence used for commenting, defaults to '//'
    :type char: string
    :return: JSON file content
    :rtype: dict
    """
    uncommented_data = ""
    for line in f:
        comment_idx = line.find(char)
        if comment_idx == -1:
            # no comment in line
            uncommented_data += line
        else:
            # remove comment from line
            uncommented_data += line[:comment_idx]
    return json.loads(uncommented_data)


def get_csv_delim(path, other_delims=set()):
    """ Get the delimiter of a character separated file.
     Checks the file for ",", "tabulator" and ";" as well as optional other characters
     or strings.
     In case no clear delimiter is found "," as delimiter is returned

    :param path: Path to file to be checked
    :type path: str
    :param other_delims: Other delimiters besides the default delimiters to be checked
    :type other_delims: set() of str
    :return: delimiter
    :rtype: str
    """

    # create union of default and optional other delimiters
    possible_delims = {",", ";", "\t"}.union(other_delims)
    # create a dict which counts the occurrences of the delimiter per row

    with open(path, "r", encoding='utf-8') as f:
        # count delimiters in first line
        line = f.readline()
        counters = {d: line.count(d) for d in possible_delims if line.count(d) > 0}

        # count delimiters in all following lines, rejecting those with different count than before
        for line_nr, line in enumerate(f):
            # for every delimiter in the dictionary. Casting to set creates new instance
            # needed in case of counter changing during the iteration.
            possible_delims = set(counters.keys())
            for delim in possible_delims:
                # compare the counted amount with the first row values
                amount = line.count(delim)
                # delete the counter if it is different to the first row
                if counters[delim] != amount:
                    del counters[delim]
                # if only one delimiter is remaining
                if len(counters) == 1:
                    # take the last item and return the key
                    return counters.popitem()[0]
                # if not even a single delimiter is remaining
                elif not counters:
                    warnings.warn("Warning: Delimiter could not be found.\n"
                                  "Returning standard Delimiter ','", stacklevel=100)
                    return ","
    #  multiple delimiters are possible. Every row was checked but more than 1 delimiter
    # has the same amount of occurrences (>0) in every row.
    warnings.warn("Warning: Delimiter could not be found.\n"
                  "Returning standard delimiter ','", stacklevel=100)
    return ","


def nd_interp(input_values, lookup_table):
    """ Get the interpolated output value for a given input value of n dimensions and a given
    lookup table with n+1 columns. Input values outside of the lookup table are
    clipped to the bounds.
    :param input_values: tuple with n input values
    :type input_values: tuple(floats)
    :param lookup_table: Table with n+1 columns, with the output values in the (n+1)th column
    :type lookup_table: list() of lists()
    :return: interpolated value
    :rtype: float
    """
    # find all unique values in table per column
    dim_sets = [set() for _ in input_values]
    for row in lookup_table:
        for i, v in enumerate(row[:-1]):
            dim_sets[i].add(v)
    dim_values = [sorted(s) for s in dim_sets]
    # find nearest value(s) per column
    # go through sorted column values until last less / first greater
    lower = [None] * len(input_values)
    upper = [None] * len(input_values)
    for i, v in enumerate(input_values):
        # initialize for out of bound values -> Constant value since lower and upper will both
        # be the same boundary value. Still allows for interpolation in other dimensions
        # forcing lower<upper could be implemented for extrapolation beyond the bounds.
        lower[i] = dim_values[i][0]
        upper[i] = dim_values[i][-1]
        for c in dim_values[i]:
            if v >= c:
                lower[i] = c
            if v <= c:
                upper[i] = c
                break
    # find rows in table made up of only lower or upper values
    points = []
    for row in lookup_table:
        for i, v in enumerate(row[:-1]):
            if lower[i] != v and upper[i] != v:
                break
        else:
            points.append(row)

    # interpolate between points that differ only in current dimension
    for i, x in enumerate(input_values):
        new_points = []
        # find points that differ in just that dimension
        for j, p1 in enumerate(points):
            for p2 in points[j + 1:]:
                for k in range(len(input_values)):
                    if p1[k] != p2[k] and i != k:
                        break
                else:
                    # differing row found
                    x1 = p1[i]
                    y1 = p1[-1]
                    x2 = p2[i]
                    y2 = p2[-1]
                    dx = x2 - x1
                    dy = y2 - y1
                    m = dy / dx
                    n = y1 - m * x1
                    y = m * x + n
                    # generate new point at interpolation
                    p = [v for v in p1]
                    p[i] = x
                    p[-1] = y
                    new_points.append(p)
                    # only couple
                    break
            else:
                # no matching row (singleton dimension?)
                new_points.append(p1)
        points = new_points

    return points[0][-1]
