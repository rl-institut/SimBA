import random
import pytest
from simba.util import nd_interp

TOLERANCE = 0.0001
POINTS_TO_CHECK = 10
DIM_AMOUNT = 4


# get the boundaries of a data_table, where the lower boundaries and upper boundaries are given
# as a tuple of two lists
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


def get_outer_point(table, dims_out_of_bound=1):
    # Number of input dimensions
    idims = len(table[0]) - 1

    # Number of dimensions which are out of bounds
    dims_out_of_bound = min(dims_out_of_bound, idims)

    # Lists for boundaries for each dimension, e.g. Dimension 2 has lower boundary of
    # low_bounds[1]
    low_bounds, upper_bounds = get_boundaries(table)

    # define a point where the first dims_out_of_bound values are outside of the boundaries, and
    # the following are inside of the boundaries
    point = ()
    for dim in range(idims):
        offset = random.random()
        if offset == 0:
            offset = 0.1
        if dim < dims_out_of_bound:
            # out of bounds
            # the point is randomly below the lower bounds or above the upper bound, depending
            # on the offset value
            if offset > 0.5:
                out_of_bounds_value = (low_bounds[dim] - offset)
            else:
                out_of_bounds_value = (upper_bounds[dim] + offset)
            point += (out_of_bounds_value,)
        else:
            # inside bounds
            scaling = offset
            # find a random value between boundaries
            inside_bounds_value = scaling * (upper_bounds[dim] - low_bounds[dim]) + low_bounds[dim]
            point += (inside_bounds_value,)
    return point



class TestNdInterpol:
    random.seed(5)
    linear_function = None
    tolerance = TOLERANCE
    points_to_check = POINTS_TO_CHECK
    dim_amount = DIM_AMOUNT

    def test_specific(self):
        nd_interp((1, 2), [(1, 2, 3), (0, 2, 4), (1, 2, 3), (0, 2, 4)])
        nd_interp((5, 10), [(1, 2, 3), (0, 2, 4), (1, 2, 3), (0, 2, 4)])

    # some manual test for a simple data table, where interpolated values can be easily
    # calculated by hand
    def test_lookup_tables(self):
        table = [(1, 2, 3), (2, 2, 4), (1, 5, 3), (2, 5, 10)]

        # perfect grid lookup
        # Looking for value for input (1,2) which is exactly found as first element in the above
        # table. Therefore 3 should be given back
        self.approx(nd_interp((1, 2), table), 3)

        # Looking for out of bounds value for both input dimension. Since the last element in
        # the above table has maximum of first and second input dimension, 10 should be given back
        self.approx(nd_interp((100, 100), table), 10)

        # looking for the value between the first and second element in the above table. Both of
        # them share the same value for the second input dimension. Therefore the interpolation
        # is one dimensional in the first dimension, i.e. f(1) = 3 and f(2)= 4 --> f(1.5) --> 3.5
        self.approx(nd_interp((1.5, 2), table), 3.5)

        # perfect grid lookup
        # making sure the perfect lookup from above works for 2 and 1 dimensions as well
        table = [(1, 2, 3), (2, 2, 4), ]

        self.approx(nd_interp((1, 2), table), 3)
        table = [(1, 2, 3)]
        self.approx(nd_interp((1, 2), table), 3)

        # checking the out of bounds for above and below for both input dimensions at the same time.
        self.approx(nd_interp((100, 100), table), 3)
        self.approx(nd_interp((-100, -100), table), 3)

    # test if out of bound dims are correctly calculated. At this stage clipping of values is
    #  assumed, e.g. constant values of the dimension onwards. Extrapolation could also be
    # possible but is not used at this point
    def test_out_of_bounds(self):
        for num_dim in range(1, self.dim_amount):
            data_table = self.generate_data_table(
                        num_dimensions=num_dim, random_values=True, dim_lengths=[2] * num_dim)
            lower_bounds, upper_bounds = get_boundaries(data_table)
            for out_of_bounds_dim in range(1, num_dim+1):
                # copy table
                stuffed_table = [p for p in data_table]

                # get a point which has one or many input values out of the table boundaries
                point = get_outer_point(data_table, dims_out_of_bound=out_of_bounds_dim)

                for dim, v in enumerate(point):
                    if lower_bounds[dim] < v < upper_bounds[dim]:
                        # if the value is inside of the boundaries, no stuffing at the borders has
                        # to take place. Therefore this dimension can be skipped
                        continue

                    if v < lower_bounds[dim]:
                        # input value is below the lower boundary of the dimension
                        value_below_bounds = True
                    else:
                        # value is not below the lower boundary. With the previous skip condition
                        # this means the value is above the upper bounds of the dimension
                        value_below_bounds = False

                    # Points which will extend the boundaries of the data_table, so that the given
                    # point with out of bound dimensions can be found in the stuffed table.
                    stuffing = []
                    # go through all rows in the data_table
                    for row in stuffed_table:
                        mutated_row = list(row)

                        # copy boundary values to out of bound values.
                        # if a row contains a boundary value which is broken by the input, the
                        # value of this row is copied to a new row, where the boundary value
                        # changes to the input value.
                        # Eg The lower boundary of dimension 0 is '5'. The data_table has a row with
                        # the value '5' for dimension 0, but the input value is '4'.
                        # In this case the row in the data_table including the output value
                        # is copied, but the value for dimension 0 is changed to '4'.
                        # In the case where the input value exceeds the boundary of the dimension,
                        # only the rows of the data_table which have the upper boundary are copied.
                        lower_bound_broken = value_below_bounds and\
                            mutated_row[dim] == lower_bounds[dim]
                        upper_bound_broken = (not value_below_bounds) and\
                            mutated_row[dim] == upper_bounds[dim]
                        if lower_bound_broken or upper_bound_broken:
                            mutated_row[dim] = v
                            stuffing.append(mutated_row)
                    stuffed_table.extend(stuffing)
                self.approx(nd_interp(point, data_table), nd_interp(point, stuffed_table))

    # test each point on the grid of an automatically and randomly generated data table.
    def test_grid_points(self):
        for i in range(1, DIM_AMOUNT):
            data_table = self.generate_data_table(num_dimensions=i, random_values=True)
            for row in data_table:
                self.approx(nd_interp(row[:-1], data_table), row[-1])

    # test 10 randomly generated points in between the table bounds. Values are checked against
    # a linear function with a TOLERANCE of 0.0001
    def test_random_points(self):
        for i in range(1, DIM_AMOUNT):
            data_table = self.generate_data_table(num_dimensions=i, random_values=False)
            for _ in range(0, POINTS_TO_CHECK):
                point = get_random_inner_point(data_table)
                exact_value = self.linear_function(point)
                self.approx(exact_value, nd_interp(point, data_table))
                self.approx(nd_interp(point, data_table), self.linear_function(point))

    # test handling of singleton dimension, e.g [1,1,11], [2,1,12] <-- second dimension only has
    # 1 value
    def test_singleton(self):
        for i in range(1, DIM_AMOUNT):
            data_table = self.generate_data_table(num_dimensions=i, random_values=True)
            singleton_value = random.random()
            for ii, row in enumerate(data_table):
                row = list(row)
                row.insert(0, singleton_value)
                data_table[ii] = row
            for row in data_table:
                self.approx(nd_interp(row[:-1], data_table), row[-1])

    # generate a data table with num_dimensions as input plus 1 output dimension
    # output values can be random values or follow a linear function
    # nr of values per dimension can be specified via dim_length
    # dimension values are not sorted and don't follow a constant stride / step size
    def generate_data_table(self, num_dimensions: int, random_values=True, dim_lengths=None):
        dims = []
        # generate dimensions with random values with random step size
        for dim in range(num_dimensions):
            start = random.randint(-100, 100)
            if dim_lengths is not None:
                steps = dim_lengths[dim]
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
                linear_factors = [random.uniform(-5, 5) for _ in range(num_dimensions)]

                def linear_function(inputs):
                    result = 0
                    for i, inp in enumerate(inputs):
                        result += inp * linear_factors[i]
                    return result

                return linear_function

            self.linear_function = get_linear_function(num_dimensions)

            output_table = [row + (self.linear_function(row),) for row in
                            list_combinator(dims)]
            return output_table

    # Wrapper for pytest.approx with constant tolerance
    def approx(self, value_a, value_b):
        tolerance = self.__class__.tolerance
        assert value_a == pytest.approx(value_b, rel=tolerance, abs=tolerance)


# creates a full matrix, if given a list of lists with unique values per dimension
def list_combinator(dimension_list):
    table = [(x,) for x in dimension_list[0]]
    for unique_values_single_dimension in dimension_list[1:]:
        new_table = []
        for v in unique_values_single_dimension:
            for ii, t in enumerate(table):
                new_table.append((*table[ii], v,))
        table = new_table
    return table


# select a random point in between each boundary of each dimensions
def get_random_inner_point(table):
    # Number of input dimensions
    idims = len(table[0]) - 1
    lower_bounds, upper_bounds = get_boundaries(table)
    point = []
    for dim in range(0, idims):
        point.append(random.random() * (upper_bounds[dim] - lower_bounds[dim]) + lower_bounds[dim],)
    return point
