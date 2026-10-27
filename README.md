# Stadium venue library

Builds a GIS layout file for every FBS, FCS, NFL, MLB and MLS stadium in the U.S. (338 venues), which the Venue Safety Plan dashboard loads from its **Stadium library** search.

## What's in it

The venue list (`venue_pipeline/venues.csv`) comes from Wikipedia's current stadium lists for each league: name, city, teams, capacity, roof type and coordinates. Edit the CSV to add, fix or remove venues. Leave `lat`/`lon` blank and the pipeline looks the venue up by name.

For each venue the pipeline pulls OpenStreetMap data and writes `venues/<slug>.geojson`:

| Layer | Where it comes from | How much to trust it |
|---|---|---|
| Stadium outline, playing field, grandstands | OpenStreetMap | Mapped, but unverified |
| Gates/entrances, elevators | OpenStreetMap | Only where volunteers mapped them. Many venues have none |
| Nearby buildings as shelter candidates | OpenStreetMap, within 350 m, at least 400 m², homes excluded | Mapped footprints. Whether a building is open and usable during events is unknown |
| Seating zones (8 per venue) | **Modeled** from the gap between field and outer outline. Capacity is spread by area | Estimate only |
| Concourse segments | **Modeled** as the outer band of the outline | Estimate only. Real concourses are usually under the stands and larger |

Every modeled feature carries `modeled: yes` and "(modeled)" in its name. The dashboard badges these features and adds a warning to the plan.

It also writes `venues/index.json`, the list the dashboard searches, and `venues/coverage.csv`, which shows what was found for each venue. Check `coverage.csv` first after a run.

### Model assumptions (edit at the top of `build_venues.py`)

- Seats start 8 m from the field edge. The concourse band is 12 m deep, or a quarter of the field-to-wall depth if that's smaller.
- Exit width per zone assumes each zone clears in 8 minutes at the Green Guide stepped-surface rate (66 people per metre per minute ≈ 20 per foot). Every zone's exit-queue time is therefore 8 minutes until you enter real widths.
- No field mapped: the field is assumed to be the outline scaled to 45%. No outline mapped: nothing is modeled.
- Domed venues get enclosed concourses and a note that lightning sheltering mainly applies to plazas, lots and tailgates.

## Running it

**Actions → Build venue library → Run workflow**. A full run takes roughly 30–60 minutes because it pauses between OpenStreetMap requests to stay within their usage policy. After that it rebuilds every Sunday at 2:17 AM Central. To rebuild a few stadiums, enter their slugs in the run box, e.g. `kroger-field-at-c-m-newton-grounds`.

If OpenStreetMap's server is down, the run stops after 8 failures in a row, keeps the existing files and shows red. Run it again later.

## Replacing modeled zones with real ones

1. In the dashboard, pick the stadium from **Stadium library**.
2. Trace the real seating sections and concourses over the venue map, fix shelter areas and exit widths, and mark what venue ops confirms as verified. Delete modeled features you don't want.
3. Click **Export GeoJSON** and upload the file to `venues/edits/` with the same name (`<slug>.geojson`).
4. On the next run, the pipeline keeps your edits. Traced seating zones replace the modeled ones, features you deleted stay deleted, and edits you made to OSM features win over fresh OSM data.

## License note

The venue files are derived from OpenStreetMap, so they're under the Open Database License (ODbL). Attribution is "© OpenStreetMap contributors". If you publish the dataset, it has to stay under ODbL. Keep restricted areas (service corridors, security posts, control rooms) out of this public repo.
