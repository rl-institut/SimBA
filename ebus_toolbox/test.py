from schedule import Schedule

schedule = Schedule.from_csv("./data/private_examples/trips_example-bvg.csv")
schedule.calculate_consumption()
print("done")
