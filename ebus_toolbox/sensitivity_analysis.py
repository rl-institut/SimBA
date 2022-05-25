""" Sensitivity analysis based on monte carlo simulation
"""
from ebus_toolbox.util import read_arguments
import pandas as pd

# Störgrößen definieren

# 1. Verspätung

verspaetung = pd.read_csv('data/monte_carlo/Fitted_Delay_at_terminus_over_weekdays.csv', sep=';', decimal='.')

# configs erstellen


def no_optimization():
    args = read_arguments()
    return "converged"


if __name__ == '__main__':
    no_optimization()

# config manipulieren/ergänzen

## gc_power_opps, gc_power_deps (Netzauslastung)
## cs_power_opps, cs_power_deps_depb, cs_power_deps_oppb (Technik)

## ergänzung: temp, verspätung, aufall hpc, mehrkilometer, (batteriealterung), zeitl. abläufe depot

# neue config speichern
