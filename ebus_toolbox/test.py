from schedule import Schedule
import os

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

schedule = Schedule.from_csv("./data/private_examples/Trips_example.csv")
