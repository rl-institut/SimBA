""" Module to generate meaningful output files and/or figures to describe simulation process
"""


def generate(schedule, scenario, args):
    schedule.generate_rotations_overview(scenario, args)
