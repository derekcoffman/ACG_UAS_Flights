# UAS_Flights

Framework for processing ACG UAS flights. 

## Requires: 
- holoviz environment

## Configuration
Edit config/settings.py to set base data folder. This is generally for envDataSystem data and actual paths can be overridden in notebook

## Use:
If multiple people are using, it would be best to create a new git branch for each Payload/Project. E.g., "keywest2022_CloudySky" would allow separate version control for all Cloudysky flights. Run helper script to create flight folders for the payload/project. This will also create a settings file specific to the flight. Also, the appropriate jupyter notebook will be created from a template using paths from settings files.

1. Create/check out payload/project branch
2. From UAS_Flights directory run uasflights in bin directory: <br>
```python bin/uasflights.py create_flight --payload="CloudySky" --project="KeyWest2022" --platform="FVR-55" --flight_id="Flight_01"```
3. Run notebook server if necessary
4. Open jupyter notebook to load/process/analyze data. Optionally, save to igor text file.
