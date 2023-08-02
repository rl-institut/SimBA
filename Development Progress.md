# SimBA: Simulation tool for Bus Applications

## Idea

This toolbox helps to do bus feasibility studies and answer research questions

## Tools

For overview and interfaces of tools see here:
https://miro.com/app/board/o9J_lu8coPI=/

### scedule data preparation

* [X] from VIP to RLI scedule data format
* [ ] from eFlips to RLI scedule data format

### scedule data analysis

* [ ] Generate table of rotation distances
* [ ] Ranking by accumulated break time

### Consumption analysis
* [ ] air conditioning consumption based on weather data

### SOC analysis
Assumptions: Leave depot fully loaded

SOC trend for each rotation based on:
* [X] charging locations (OC)
* [X] number of available charging infrastructure (OC)
* [ ] charging characteristics (OC)

### Charge demand analysis
Generate load profile for each charging location based on SOC analysis

### Schedule adjustment
* [X] split up rotations with negative SOC
* [X] recombine splitted rotations

### Charging Location Tool