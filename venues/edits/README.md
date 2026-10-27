# Hand edits

Put corrected or traced layouts here as `<venue-slug>.geojson` (export them from the Venue Safety Plan dashboard, or from RustGIS in WGS84).

On every pipeline run:
- a feature here replaces the generated feature with the same id
- ids listed in the file's `removed_ids` are dropped from the generated layout
- if you traced real seating zones, the modeled seating zones are dropped; if you traced real concourses, the modeled concourses are dropped

Nothing in this folder is ever overwritten by the pipeline.
