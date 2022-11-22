"""
run these tests with `pytest tests/test_something.py` or `pytest tests` or simply `pytest`
pytest will look for all files starting with "test_" and run all functions
within this file. For basic example of tests you can look at our workshop
https://github.com/rl-institut/workshop/tree/master/test-driven-development.
Otherwise https://docs.pytest.org/en/latest/ and https://docs.python.org/3/library/unittest.html
are also good support.
"""

import random

from ebus_toolbox.util import nd_interp

random.seed(5)
linear_function = None
TOLERANCE = 0.0001
POINTS_TO_CHECK = 10
DIM_AMOUNT = 4


def test_lookup_tables():
    table = [(1, 2, 3), (2, 2, 4), (1, 5, 3), (2, 5, 10)]
    assert nd_interp((1, 2), table) == 3
    assert nd_interp((100, 100), table) == 10
    assert nd_interp((1.5, 2), table) == 3.5
    table = [(1, 2, 3), (2, 2, 4), ]
    assert nd_interp((1, 2), table) == 3
    table = [(1, 2, 3)]
    assert nd_interp((1, 2), table) == 3
    assert nd_interp((100, 100), table) == 3
    assert nd_interp((-100, -100), table) == 3


def get_boundaries(table):
    idims = len(table[0]) - 1
    # get the boundaries of each dimension
    lower_bounds = [float('inf')] * idims
    upper_bounds = [float('-inf')] * idims
    for row in table:
        for dim in range(0, idims):
            lower_bounds[dim] = min(lower_bounds[dim], row[dim])
            upper_bounds[dim] = max(upper_bounds[dim], row[dim])
    return lower_bounds, upper_bounds


# Test if out of bound dims are correctly calculated. At this stage clipping of values is assumed,
# e.g. constant values of the dimension onwards. Extrapolation could also be possible but isnt
# used at this point
def test_out_of_bounds():
    for num_dim in range(1, DIM_AMOUNT):
        data_table = generate_data_table(num_dimensions=num_dim, random_values=True,
                                         dim_length=[2] * num_dim)
        lower_bounds, upper_bounds = get_boundaries(data_table)
        for out_of_bounds_dim in range(1, num_dim):
            # copy table
            stuffed_table = [p for p in data_table]
            point = get_outer_point(data_table, dims_out_of_bound=out_of_bounds_dim)
            for row in data_table:
                for dim, v in enumerate(point):
                    row = list(row)
                    if (v < lower_bounds[dim] and row[dim] == lower_bounds[dim]) \
                            or (v > upper_bounds[dim] and row[dim] == upper_bounds[dim]):
                        row[dim] = v
                stuffed_table.append(tuple(row))
            assert nd_interp(point, data_table) \
                   == nd_interp(point, stuffed_table)


# select a random point in between each boundary of each dimensions
def get_outer_point(table, dims_out_of_bound=1):
    idims = len(table[0]) - 1
    dims_out_of_bound = min(dims_out_of_bound, idims)
    low_bounds, upper_bounds = get_boundaries(table)

    # define a point with the number of dims_out_of_bound values which are not inside
    # the dimension values
    point = ()
    for dim in range(0, idims):
        offset = random.random()
        # out of bounds
        if dim < dims_out_of_bound:
            point += ((low_bounds[dim] - offset) + (offset > 0.5) *
                      ((-low_bounds[dim] + offset) + (upper_bounds[dim] + offset)),)
        else:
            # inside bounds
            point += (offset * (upper_bounds[dim] - low_bounds[dim]) + low_bounds[dim],)
    return point


# test each point on the grid of an automatically and randomly generated data table.
def test_grid_points():
    for i in range(1, DIM_AMOUNT):
        data_table = generate_data_table(num_dimensions=i, random_values=True)
        for row in data_table:
            assert nd_interp(row[:-1], data_table) == row[-1]


# test 10 randomly generated points in between the table bounds. Values are checked against
# a linear function with a TOLERANCE of 0.0001
def test_random_points():
    global linear_function
    for i in range(1, DIM_AMOUNT):
        data_table = generate_data_table(num_dimensions=i, random_values=False)
        for _ in range(0, POINTS_TO_CHECK):
            point = get_random_inner_point(data_table)
            exact_value = linear_function(point)
            # get signature of value for tolerance (sig=+1 or -1)
            sig = 1 - (exact_value < 0) * 2
            assert ((1 - TOLERANCE * sig) * linear_function(point)) <= nd_interp(point, data_table)
            assert nd_interp(point, data_table) <= ((1 + TOLERANCE * sig) * linear_function(point))


# select a random point in between each boundary of each dimensions
def get_random_inner_point(table):
    idims = len(table[0]) - 1
    lower_bounds, upper_bounds = get_boundaries(table)
    point = ()
    for dim in range(0, idims):
        point += (random.random() * (upper_bounds[dim] - lower_bounds[dim]) + lower_bounds[dim],)
    return point


def test_singleton():
    for i in range(1, DIM_AMOUNT):
        data_table = generate_data_table(num_dimensions=i, random_values=True)
        singleton_value = random.random()
        for ii, row in enumerate(data_table):
            row = list(row)
            row.insert(0, singleton_value)
            data_table[ii] = row
        for row in data_table:
            assert nd_interp(row[:-1], data_table) == row[-1]


# generate a data table with num_dimensions as input plus 1 output dimension
# output values can be random values or follow a linear function
# nr of values per dimension can be specified via dim_length
# dimension values are not sorted and dont follow a constant stride / step size
def generate_data_table(num_dimensions: int, random_values=True, dim_length=None):
    dims = []
    # generate dimensions with random values with random step size
    for _ in range(0, num_dimensions):
        start = random.randint(-100, 100)
        if dim_length is not None:
            steps = dim_length[_]
        else:
            steps = random.randint(1, 5)
        dims.append(list(
            {start + (random.randint(1, 4) * random.randint(-10, 10)) * k for k in
             range(0, steps)}))
    # combine dimensions via list_combinator and append random or linear value as output to each
    # row
    if random_values:
        output_table = [row + (random.randint(1, 4) * random.randint(1, 10),) for row in
                        list_combinator(dims)]
        return output_table
    else:
        def get_linear_function(num_dimensions):
            linear_factors = [random.random() * random.randint(-5, 5) for _ in
                              range(0, num_dimensions)]

            def linear_function(inputs):
                result = 0
                for i, inp in enumerate(inputs):
                    result += inp * linear_factors[i]
                return result

            return linear_function

        global linear_function
        linear_function = get_linear_function(num_dimensions)

        output_table = [row + (linear_function(row),) for row in
                        list_combinator(dims)]
        return output_table


def list_combinator(dimension_list):
    table = [(x,) for x in dimension_list[0]]
    for i, l in enumerate(dimension_list):
        if i == 0:
            continue
        new_table = []
        for v in l:
            for ii, t in enumerate(table):
                new_table.append((*table[ii], v,))
        table = new_table
    return table
