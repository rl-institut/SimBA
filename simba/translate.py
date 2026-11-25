language = None
translation_tables = {
    "de": {
        "Active Rotations": "Aktive Umläufe",
        "bat. stored energy [kWh]": "stat. Speicher Energieinhalt in kWh",
        "battery power [kW]": "stat. Speicher Leistung in kW",
        "Block": "Umlauf",
        "Charging type": "Ladetyp",
        "Depot": "Depot",
        "Distance [km]": "Distanz in km",
        "Distribution of energy consumption of rotations per vehicle type":
            "Verteilung des Energiebedarfs über Umläufe und Fahrzeugtypen",
        "Distribution of rotation length per vehicle type":
            "Verteilung der Umlauflänge je Fahrzeugtyp",
        "Energy consumption [kWh]": "Energiebedarf in kWh",
        "Feasibility of rotations per charging type": "Machbarkeit der Umläufe je Ladetyp",
        "fixed load [kW]": "Feste Last in kW",
        "grid supply [kW]": "Netzanschlusspunkt in kW",
        "Inside window": "Innerhalb des HLZ",
        "local generation [kW]": "lokale Erzeugung in kW",
        "negative rotations": "Negative Umläufe",
        "Number of active rotations": "Anzahl aktiver Umläufe",
        "Number of rotations": "Anzahl Umläufe",
        "Opportunity": "Gelegenheit",
        "Outside window": "Außerhalb des HLZ",
        "Power": "Leistung",
        "Power [kW]": "Leistung in kW",
        "price [ct/kWh]": "Strompreis in ct/kWh",
        "stored battery energy [kWh]": "Batterie Energieinhalt in kWh",
        "successful rotations": "Erfolgreiche Umläufe",
        "sum CS power [kW]": "Ladestationen in kW",
        "Vehicle ID": "Fahrzeug ID",
    },
}


def set_language(new_language):
    global language
    language = new_language


def translate(text):
    if language in translation_tables:
        # language set and known
        table = translation_tables[language]
        # find text in translation table, otherwise return text unchanged
        return table.get(text, text)
    # language not set or not known
    return text
