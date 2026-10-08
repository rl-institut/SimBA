# Konzept: Kostenkorrektur für Speicherfüllstand-Änderungen

Status: **Implementiert auf Branch `feature/rebase_el_cost` (noch nicht committet/gemerged)**
Betrifft: `simba/costs.py` (`Costs.set_storage_level_cost_correction` und Helfer),
Konfiguration (`simba/util.py::get_args`, Flag `clean_costs_storage_level`),
`data/examples/configs/extensive.cfg`, `docs/source/simba_features.rst`

## 0. Umsetzung (Abweichungen vom Konzept unten)

Die Implementierung folgt Stufe C wie unten beschrieben, mit einer wichtigen praktischen
Korrektur gegenüber dem ursprünglichen Konzept: Für Stationen mit einer `fixed_*`
`cost_calculation_method` liefert SpiceEVs Echtzeit-Preissignal (`scenario.prices`) **keinen**
brauchbaren Preis – jeder GC bekommt beim Szenario-Aufbau einen hartkodierten
Platzhalterpreis (`{"type": "fixed", "value": 0.3}`, siehe `simba/schedule.py:1009`), der nur
der Lade-Strategie zur Echtzeit-Entscheidung dient, nicht der eigentlichen (nachträglichen)
Preisblatt-basierten Abrechnung. Für `fixed_*`-Methoden wird deshalb statt des Echtzeitsignals
ein aus den bereits berechneten (Preisblatt-basierten) Energiekosten dieses GC abgeleiteter
Flatrate-Preis verwendet (`annual_cost / annual_kwh`, beides bereits annualisiert). Für
`variable_*`/`balanced_market`/`flex_window` wird weiterhin das echte Preissignal
(CSV-Preise oder `scenario.prices`) verwendet. Diese Fallunterscheidung ist notwendig, damit
die Korrektur überhaupt einen realistischen (nicht nahe Null liegenden) Preis verwendet.

An einem konkreten Testlauf verifiziert (`simba_w_bat.cfg`/`simba_wo_bat.cfg`,
`fixed_w_plw`, nach dem in Abschnitt 0.2 beschriebenen Effizienz-Fix): Die jährlichen
Stromkosten des Batterie-Szenarios sinken durch die Korrektur von 837.700 €/a auf 127.766 €/a
(vorher fälschlich höher als das batterielose Szenario mit 205.583 €/a; nach Korrektur beider
Szenarien liegt das Batterie-Szenario bei 127.766 €/a gegenüber 185.954 €/a ohne Batterie – ein
plausibles Bild, das für eine funktionierende Peak-Shaving-Batterie zu erwarten ist).

### 0.1 Automatisierte Tests

In `tests/test_cost_calculation.py` (drei Ebenen, orientiert am bestehenden Testaufbau des
Moduls):

- `TestStorageCorrectionUnit`: reine, deterministische Tests für
  `Costs._storage_correction_eur` ohne jede Simulation – Regressionstest für den
  Effizienz-Bug (0.2), Defizit-Fall, mehrstufiger Überschuss mit wechselndem Preis (vergleicht
  gleitenden Durchschnitt explizit gegen einen naiven Blend-Durchschnitt), Randfälle (nie
  geladen, ΔE=0, führende/nachgestellte `None`-Werte, fehlender Preis an einem Ladeschritt).
- `TestStorageCorrectionHelpers`: `_last_rotation_gc_per_vehicle` und `_vehicle_price_series`
  mit leichten Stand-ins (`SimpleNamespace`) statt echter Simulation.
- `TestStorageCorrectionIntegration`: mit echtem `(schedule, scenario)` aus
  `BasicSchedule().basic_run()` (öffentliche Beispieldaten). Da die dort verwendete
  `balanced`-Strategie die stationäre Batterie nie zyklisiert, wird für diese Tests ein
  bekannter Überschuss direkt in `scenario.batteryLevels` injiziert, um Verkabelung (Flag an/aus,
  `c_el_annual = vorher + Korrektur`), CSV-Zeile und Log-Meldung deterministisch zu prüfen.

### 0.2 Gefundener und behobener Bug: Wirkungsgrad im Überschuss-Zweig doppelt angewendet

Beim ersten Testlauf traten unplausibel große negative Korrekturen auf (z. B. Station A:
182.234 €/a → -6.103 €/a), obwohl der Tarif nachweislich exakt fix war (0,255220 €/kWh bei
allen vier Stationen identisch). Ursache: `avg_price` im gleitenden Durchschnitt (Überschuss-Fall)
konvergiert korrekt auf einen Preis **pro gespeicherter Batterie-kWh** (das ist bereits
`Netzpreis / η`, da 1 kWh im Speicher `1/η` kWh Netzbezug kostet). Die finale Zeile
`-(delta_e / efficiency) * avg_price` hat diesen Wirkungsgrad aber ein zweites Mal angewendet,
weil `delta_e` dort nochmal auf Netz-Äquivalent umgerechnet wurde – macht ihn effektiv `Preis / η²`
statt `Preis / η`. Fix: `avg_price` wird jetzt konsistent als "Preis pro Batterie-kWh" geführt
(Startwert `own_avg_price / η` statt `own_avg_price`) und am Ende ohne nochmalige
Effizienz-Division verwendet (`-delta_e * avg_price`). Der Defizit-Zweig war davon nicht
betroffen (dort ist `own_avg_price` durchgehend ein Preis pro Netz-kWh, konsistent mit der
Multiplikation mit `-delta_e / η`).

Nach dem Fix (Testlauf `simba_w_bat.cfg`): Station A liegt korrigiert bei +3.314 €/a (statt
-6.103 €/a), alle vier Stationen sind nach Korrektur positiv, kumuliert 837.700 €/a → 127.766 €/a
(statt fälschlich 90.098 €/a vor dem Fix). Regressionstest:
`TestStorageCorrectionUnit.test_surplus_single_charge_efficiency_not_double_counted`.

## 1. Problem

Beim Vergleich zweier sonst identischer Szenarien (mit vs. ohne stationäre Batterie,
`data/private_examples/VHH_Speicher/HLZ/simba_w_bat.cfg` vs. `simba_wo_bat.cfg`) stiegen die
jährlichen Stromkosten mit Batterie unerwartet stark an (837.700 €/a vs. 205.583 €/a), obwohl eine
Peak-Load-Batterie tendenziell Kosten senken sollte.

### Ursache (im konkreten Fall bereits identifiziert)

1. SpiceEVs `peak_load_window`-Strategie lädt stationäre Batterien außerhalb der Lastfenster immer
   vollständig auf (`spice_ev/strategies/peak_load_window.py`), unabhängig davon, ob das für die
   Spitzenkappung im nächsten Fenster überhaupt nötig ist.
2. Die Batterien starten in `electrified_stations_hlz_w_bat.json` leer (`"soc": 0`) mit
   `"capacity": 5000` kWh und `"efficiency": 0.95`. Das einmalige Vollladen zieht
   `5000 / 0.95 ≈ 5263 kWh` Netzenergie **pro Station, einmalig**.
3. Der Testzeitraum (`trips_example_small.csv`) umfasst nur ca. 3 Tage. Der einmalige
   Ladevorgang macht dadurch einen unverhältnismäßig großen Anteil der gemessenen Energie aus.
4. `simba/costs.py` skaliert Energiekosten linear aufs Jahr
   (`annual_factor = (365*24*3600) / scenario_duration_s`, siehe `costs.py:500-503`). Bei 3
   simulierten Tagen ist `annual_factor ≈ 122`. Der einmalige Ladevorgang wird dadurch behandelt,
   als träte er ~122-mal im Jahr auf.

### Verallgemeinerung

Das Problem ist nicht auf diesen Testfall beschränkt: **Jede** Simulation, bei der der
Energieinhalt der Speicher (stationäre Batterien **und** Fahrzeugbatterien) am Ende der
Simulation vom Anfangszustand abweicht, verzerrt die hochskalierten (annualisierten)
Energiekosten – besonders bei kurzen Simulationszeiträumen. Das erschwert den fairen Vergleich
zweier Szenarien, wenn sich deren Speicher unterschiedlich stark auf- oder entladen.

## 2. Ziel

Eine optionale Nachbehandlung der Kostenrechnung, die den Effekt von Speicherfüllstand-Änderungen
(SoC-Drift) auf die **energieproportionalen** Kostenanteile herausrechnet, sodass Simulationen mit
unterschiedlichem Speicherverhalten am Ende und Anfang der Simulation fair vergleichbar bleiben –
unabhängig von der Simulationsdauer.

Vorschlag für die Konfiguration:

```
clean_costs_storage_level = true
```

als zusätzliches Flag neben `cost_calculation` in der SimBA-Config.

## 3. Was wird korrigiert – und was nicht

Kosten lassen sich in zwei Kategorien teilen:

| Kategorie | Beispiele | Betroffen von SoC-Drift? |
|---|---|---|
| Energieproportional (€/kWh) | Strombeschaffung, energiebezogene Netzentgelte, Stromsteuer, Einspeisevergütung | **Ja** – wird korrigiert |
| Leistungs-/kapazitätsproportional (€/kW oder €/Jahr fix) | Leistungspreis (Peak Power), Investition/Wartung von GC, Ladepunkten, Speicherkapazität, Fahrzeugen | Nein – bleibt unverändert |

Die Korrektur greift ausschließlich in die energieproportionalen Kostenanteile ein
(`c_el_procurement_annual`, `c_el_energy_price_annual`, `c_el_taxes_annual`,
`c_el_feed_in_remuneration_annual` in `simba/costs.py`), nicht in `c_el_power_price_annual` oder
irgendeine Investitions-/Wartungsposition.

## 4. Grundprinzip der Korrektur

Für jeden Speicher (jede stationäre Batterie, jede Fahrzeugbatterie) wird die Differenz des
Energieinhalts zwischen Simulationsende und -anfang bestimmt:

```
ΔE = (SoC_Ende − SoC_Anfang) × Kapazität
```

- **ΔE > 0** (am Ende mehr gespeichert): Ein Teil der bezogenen Energie wurde nicht für den
  Fahrbetrieb "verbraucht", sondern steckt noch im Speicher. Die dafür nötige Netzenergie
  (`ΔE / η_charge`, unter Berücksichtigung des Ladewirkungsgrads) wird virtuell **nicht gekauft**
  → Kosten werden reduziert.
- **ΔE < 0** (am Ende weniger gespeichert): Der Speicher hat Energie aus einem vorherigen,
  nicht simulierten Zyklus verbraucht. Um den Ursprungszustand wiederherzustellen, müsste
  `|ΔE| / η_charge` virtuell **nachgekauft** werden → Kosten werden erhöht.

Wichtig: Es wird in beiden Fällen mit dem **Einkaufspreis** gerechnet, nicht mit der
Einspeisevergütung – ein Überschuss wird nicht wie eine echte Netzeinspeisung behandelt, sondern
wie ein vermiedener Kauf. Die Vorzeichen sind symmetrisch, die Preisbasis ist immer dieselbe.

## 5. Ausbaustufen (von einfach zu genau)

### Stufe A – Systemweiter Durchschnittspreis (verworfen)

Ein einziger, mengengewichteter Durchschnittspreis über alle GCs. Einfach, aber ungenau: jeder GC
kann an einem anderen Netzbetreiber mit eigenem Tarif hängen, und das widerspricht der ohnehin pro
GC ausgewiesenen Kostenstruktur in `summary_vehicles_costs.csv`.

### Stufe B – Durchschnittspreis je GC (verworfen als Endzustand, Vorstufe zu C)

Pro Grid Connector wird ein mengengewichteter Durchschnittspreis aus allen energieproportionalen
Kostenkomponenten gebildet (`(c_el_procurement + c_el_energy_price + c_el_taxes) / bezogene_kWh`),
gewichtet nach tatsächlich bezogener Menge, nicht nach Zeit, und auf alle Geräte an diesem GC
gleichermaßen angewendet.

**Bekannte Schwäche:** Ein einzelner Blend-Preis über *alle* an einem GC bezogene Energie
(Fahrzeuge **und** stationäre Batterie gemeinsam) verzerrt das Ergebnis, wenn ein Teil der
Ladevorgänge (typischerweise die flexible/optionale Batterieladung) gezielt in
Niedrigpreisfenstern stattfindet. Der "Überschuss" am Ende stammt dann überproportional aus
günstigen Zeiten, wird aber mit dem (höheren) Blend-Durchschnitt bewertet → Korrektur schreibt zu
viel gut. Bei `fixed_w_plw`/`fixed_wo_plw` (siehe unten: zeitlich konstanter Preis) ist dieser
Fehler irrelevant; bei marktbasierten Preisen (`variable_*`, `balanced_market`) potenziell nicht
mehr, da eine Batterie ja gerade *deshalb* eingesetzt wird, um Preisvarianz gezielt auszunutzen.
Zusätzlich war die Attribution von Fahrzeugbatterien auf mehrere GCs (anteilig nach geladener
kWh) hier noch ungelöst komplex. Beides wird in Stufe C aufgelöst.

### Stufe C – Konkretes, umsetzbares Konzept (empfohlen)

Stufe C behandelt Überschuss und Defizit **asymmetrisch**, weil es sich um zwei unterschiedliche
Fragestellungen handelt, und trennt außerdem **Preisberechnung** (welcher Preis gilt für die
ΔE-Menge eines Geräts) von **Attribution** (welchem GC wird der resultierende €-Betrag in der
Kostenausgabe zugeschlagen).

#### 5.C.1 Asymmetrische Bewertung von Überschuss und Defizit

| Fall | Frage | Methode |
|---|---|---|
| ΔE > 0 (Überschuss) | "Was ist die Energie wert, die aktuell noch im Speicher steckt?" | **Gleitender Durchschnittspreis** (Lagerbewertung) über die eigene Lade-/Entladehistorie des Geräts |
| ΔE < 0 (Defizit) | "Was kostet es typischerweise, dieses Gerät nachzuladen?" | **Einfacher mengengewichteter Durchschnitt** aller eigenen Ladevorgänge des Geräts über die gesamte Simulation |

Der Überschussfall ist eine Lagerbewertungsfrage (was ist aktuell noch vorhanden – dafür ist der
gleitende Durchschnitt das richtige Werkzeug, weil er zwischenzeitlich bereits wieder verbrauchte
billige Energie korrekt herausrechnet). Der Defizitfall ist dagegen eine Wiederbeschaffungsfrage
("was zahlt dieses Gerät normalerweise beim Laden") – dafür ist ein gleitender Durchschnitt nicht
sinnvoll, weil nichts "im Lager" bewertet werden muss; der schlichte Durchschnitt über die eigene
Ladehistorie ist hier die aussagekräftigere Schätzung.

Dieselbe Logik gilt für stationäre Batterien wie für Fahrzeugbatterien. Bei der stationären
Batterie ist praktisch die *gesamte* Ladehistorie "flexibel/optional" – gerade deshalb ist der
gleitende Durchschnitt dort besonders wichtig, sobald sie über die Simulation mehrfach lädt und
entlädt (zyklisiert).

**Gleitender Durchschnitt – Algorithmus je Gerät (Stationär-Batterie oder Fahrzeug):**

1. Startbestand = `SoC_Start × Kapazität` (kWh).
2. **Startwert-Konvention (entschieden):** Der Startbestand wird mit dem eigenen
   mengengewichteten Durchschnittspreis dieses Geräts (siehe Defizit-Spalte oben, berechnet über
   die gesamte Simulation) bewertet. Das ist konsistent (keine zusätzliche freie Annahme) und die
   neutralste verfügbare Schätzung, da der tatsächliche Preis vor Simulationsbeginn nicht bekannt
   ist.
3. Für jeden Zeitschritt, in dem das Gerät **lädt**: zugeführte Energiemenge (Netzenergie,
   also `ΔSoC × Kapazität / η_charge`) zum zu diesem Zeitpunkt am zugehörigen GC geltenden Preis
   einbuchen; neuer Bestandsdurchschnitt = gewichtetes Mittel aus Altbestand und Zukauf.
4. Für jeden Zeitschritt, in dem das Gerät **entlädt**: Menge vom Bestand abziehen, Durchschnitt
   bleibt unverändert.
5. Am Simulationsende: `ΔE = Endbestand − Startbestand` (kWh, energieseitig am Speicher).
   - Wenn `ΔE > 0`: Korrektur = `− (ΔE / η_charge) × Bestandsdurchschnittspreis_Ende`
     (vermiedener Kauf).
   - Wenn `ΔE < 0`: Korrektur = `+ (|ΔE| / η_charge) × Ø-Preis_eigene_Ladehistorie`
     (angenommener Nachkauf).

**Warum ohne zweite Simulation:** Alle benötigten Zeitreihen liegen aus dem einen bereits
gelaufenen Simulationslauf vor (siehe 5.C.3). Die Buchführung ist ein zusätzlicher
Nachbearbeitungsschritt (ein Durchlauf über vorhandene Zeitreihen je Gerät), keine erneute
Simulation.

**Verbleibende Näherung:** Der gleitende Durchschnitt vermischt bei mehreren Zyklen leicht ältere
und neuere Ladevorgänge (statt exakter FIFO/LIFO-Zuordnung) – das ist der in der Praxis übliche,
für diesen Zweck ausreichend genaue Kompromiss.

#### 5.C.2 Attributionsregel (Vorschlag Julian, übernommen)

Die **Preisberechnung** (5.C.1) läuft je Gerät unabhängig davon, an welchem GC tatsächlich
geladen wurde. Für die **Zuordnung des resultierenden €-Betrags** zu einer Spalte in
`summary_vehicles_costs.csv` gilt:

- **Stationäre Batterie:** eindeutig einem GC zugeordnet → trivial.
- **Fahrzeugbatterie:** Der gesamte Korrekturbetrag eines Fahrzeugs wird dem GC zugeschlagen, an
  dem der **zuletzt gestartete Umlauf** dieses Fahrzeugs begonnen hat. Zwischenladungen an
  anderen Stationen werden für die Attribution nicht berücksichtigt (fließen aber sehr wohl in
  die Preisberechnung nach 5.C.1 ein, da dort die gesamte Ladehistorie zählt).

  Konkret in den vorhandenen Datenstrukturen: `Schedule.rotations` (Dict `rotation_id ->
  Rotation`) nach `vehicle_id` filtern, die Rotation mit dem größten `departure_time` wählen,
  deren `departure_name` ist der Ziel-GC (jede Rotation beginnt und endet laut
  Konsistenzannahme an derselben Station, siehe `docs/source/simba_features.rst`).

  Da Fahrzeuge in SimBA nur für tatsächlich zugewiesene Rotationen instanziiert werden
  (`Schedule.assign_vehicles`), hat praktisch jedes Fahrzeug mindestens einen Umlauf im
  Simulationszeitraum – der Fall "Fahrzeug ohne jeden Umlauf" sollte nicht auftreten.
- Die bestehende Nicht-Doppelzählungs-Regel für die Spalte "cumulated" (Busse mit mehreren
  Depots werden dort nur einmal gezählt, siehe `docs/source/simba_features.rst`) gilt
  unverändert auch für die Korrekturbeträge.

#### 5.C.3 Datenquellen (bereits vorhanden, keine SpiceEV-Änderung nötig)

- **Stationäre Batterie – Energieinhalt je Zeitschritt:** `scenario.batteryLevels[batName]`
  (kWh, siehe `spice_ev/scenario.py:178`). Differenz zwischen aufeinanderfolgenden Zeitschritten
  = Lade-/Entladeleistung dieses Schritts (Batterien laden oder entladen nie beides gleichzeitig
  im selben Zeitschritt in diesem Modell, daher ist die Differenz eindeutig einer Richtung
  zuordenbar).
- **Fahrzeugbatterie – SoC-Verlauf:** aus `scenario.socs` (Werte während Verbindung an CS) und
  `scenario.disconnect` (interpolierte Werte während Fahrt) zusammensetzen – exakt die Logik,
  die bereits für die `vehicle_SoCs.csv`-Ausgabe in `report.py` existiert und wiederverwendet
  werden sollte, statt neu implementiert zu werden.
- **Zuordnung Fahrzeug → GC beim Laden:** `scenario.connected[t]` (Dict `vehicle_id -> cs_id`)
  kombiniert mit `charging_stations[cs_id].parent` (GC-ID) und `scenario.connChargeByTS[gcID][t]`
  (Ladeleistung je `cs_id`).
- **Preis je Zeitschritt und GC:** `getattr(scenario, f"{gcID}_timeseries")["price [EUR/kWh]"]`
  (siehe `simba/costs.py:325-328`) – dieselbe Preis-Zeitreihe, die auch für die reguläre
  Kostenberechnung verwendet wird, keine neue Datenquelle nötig.
- **Wichtige Randbedingung:** Bei `fixed_wo_plw`/`fixed_w_plw` kommt der Preis aus einem
  statischen Preisblatt (`spice_ev/costs.py`, `price_sheet`) und ist über die gesamte Simulation
  **konstant**. Für diese (aktuell verwendeten) Cost-Calculation-Methoden ist der gleitende
  Durchschnitt daher rechnerisch identisch mit einem einfachen Durchschnitt – die Sophistication
  von Stufe C zahlt sich erst bei `variable_*`/`balanced_market`/`flex_window` wirklich aus. Der
  Algorithmus muss deshalb aber nicht speziell für den Fixpreis-Fall verzweigt werden: er deckt
  ihn automatisch als (trivialen) Spezialfall ab.

### Stufe D – Iterative Einschwing-Simulation (nicht empfohlen)

Statt die Kosten nachträglich zu korrigieren, könnte man vor dem eigentlichen Lauf den Start-SoC
aller Speicher iterativ an den End-SoC angleichen ("Rebasing"), sodass `ΔE ≈ 0` gilt und sich die
Vergleichbarkeit ohne nachträgliche Korrektur ergibt. Das ist physikalisch am saubersten, aber
deutlich teurer (mehrere zusätzliche vollständige Simulationsläufe) und für den hier verfolgten
Zweck (reine Kostenvergleichbarkeit) unverhältnismäßig aufwendig gegenüber Stufe C.

## 6. Architekturentscheidung: SimBA oder SpiceEV?

**Empfehlung: SimBA**, konkret in `simba/costs.py`, aus folgenden Gründen:

- Die Annualisierung (`annual_factor`) und die Trennung der Kostenarten
  (energie- vs. leistungsproportional) passiert bereits dort – die Korrektur muss vor genau dieser
  Skalierung ansetzen.
- Alle nötigen Rohdaten (Batterie-SoC-Zeitreihen der stationären Speicher, Fahrzeug-SoC je
  Rotation, Kapazitäten, Preis-Zeitreihen je GC) sind über das von SpiceEV gelieferte
  `scenario`-Objekt bereits zugänglich, ohne dass SpiceEV etwas zusätzlich exportieren müsste.
- Der Konfigurations-Flag gehört ohnehin zu SimBA (`simba/util.py::get_args`, analog zu
  `cost_calculation`).

Ein Eingriff in SpiceEV selbst wäre nur nötig, wenn die Korrektur direkt in die
Kostenberechnung von SpiceEV (`spice_ev.costs.calculate_costs`) integriert werden sollte – dafür
besteht aber kein Grund, da SimBA bereits Zugriff auf alle notwendigen Zeitreihen hat.

## 7. Offene Punkte

Geklärt und umgesetzt (siehe 5.C und Abschnitt 0): Datenquellen, Attributionsregel für
Fahrzeugbatterien, asymmetrische Bewertung Überschuss/Defizit, Startwert-Konvention des
gleitenden Durchschnitts, Preisquelle je nach `cost_calculation_method` (fixed vs. variable),
Reihenfolge der Korrektur relativ zur `annual_factor`-Skalierung (Korrektur wird auf Basis der
simulierten Periode berechnet und dann mit demselben `annual_factor` skaliert wie die übrigen
Energiekostenposten), Fahrzeug-SoC-Rekonstruktion (`spice_ev.report.generate_soc_timeseries`
wird bei Bedarf selbst angestoßen, falls `scenario.vehicle_socs` zum Zeitpunkt der
Kostenberechnung noch nicht existiert – siehe Abschnitt 0).

Erledigt: Automatisierte Tests in `tests/test_cost_calculation.py` (siehe Abschnitt 0.1).

Verbleibend:

- Negative "corrected" `c_el_annual`-Werte pro GC sind grundsätzlich möglich, wenn die
  herausgerechnete Anomalie größer ist als die ursprünglich berechneten Kosten dieses GC. Das
  wäre mathematisch korrekt (eine Vergleichsgröße, keine reale Rechnung), sollte aber in der
  Doku/UI klar so kommuniziert bleiben, damit es nicht als Fehler missverstanden wird.
- Flatrate-Preis-Fallback (`annual_cost / annual_kwh` aus `sum(scenario.totalLoad)`) kann bei
  GCs mit nennenswerter lokaler Einspeisung (PV/Feed-in) unerwartete Werte liefern, da
  `totalLoad` dort durch negative Lastanteile verkleinert bzw. verzerrt wird. Beobachtet im
  öffentlichen Test-Fixture (`data/examples/electrified_stations/electrified_stations.json`,
  Station-0 mit PV-Einspeisung): die Korrektur kam dort mit unerwartet großem positivem Vorzeichen
  heraus. Die Verkabelung (`c_el_annual = vorher + Korrektur`) bleibt davon unberührt und ist
  getestet; die *Höhe* der Korrektur bei GCs mit Einspeisung ist aber noch nicht separat
  verifiziert/diskutiert.
- Kein Testfall mit `variable_*`/`balanced_market` auf Integrationsebene (nur die
  Fixpreis-Fallback-Logik wird dort durchlaufen); die gleitende-Durchschnitt-Logik selbst ist
  aber bereits durch die reinen Unit-Tests in Abschnitt 0.1 (unabhängig von der
  Cost-Calculation-Methode) abgedeckt.

## 8. Nicht Teil dieses Konzepts

- Änderungen an der eigentlichen Ladephysik/-strategie in SpiceEV.
- Eine Automatisierung, die den Start-SoC von Speichern realistisch vorbelegt (das würde das
  ursprüngliche Symptom mindern, aber nicht das allgemeinere Vergleichbarkeitsproblem lösen).
