# Earthquake Signal GUI

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![GUI](https://img.shields.io/badge/GUI-DearPyGui-ff69b4.svg)](https://github.com/hoffstadt/DearPyGui)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

Interactive GUI for earthquake signal analysis from MATLAB (.mat) files. Provides P/S picking, polarity, azimuth, S–P distance, epicenter estimation on a Turkey map, and Excel export.

---

- [Features](#features)
- [Quick Start](#quick-start)
- [Screenshot](#screenshot)
- [Data Format](#data-format)
- [Original vs Calculated](#original-vs-calculated)
- [Directory Structure](#directory-structure)
- [Troubleshooting](#troubleshooting)
- [License](#license)

## Features

- Multi-channel visualization (North, East, Up)
- P/S wave picking with synchronized drag-lines
- P-wave polarity detection (up/down)
- Azimuth and back-azimuth, event azimuth
- S–P distance and epicenter estimation on map
- File search/filter with P/S status indicators
- Bulk Excel export

## Quick Start

1) Install dependencies
```bash
pip install -r requirements.txt
```

2) Add your `.mat` files to `sample_data/`
- At least 5 recordings are recommended
- The app defaults to `sample_data/` as the data folder

3) Launch the app
```bash
python earthquake_signal_gui.py
```

## Screenshot

![GUI Screenshot](assets/screenshot.png)

If the image does not render on GitHub, ensure the file exists at `assets/screenshot.png` and is committed.

## Data Format

The GUI expects one of the following .mat structures (either is fine):

- Under `EQ.anEQ`:
  - `Accel` (N×3 array for North, East, Up)
  - `Ptime` (seconds), `Stime` (seconds)
  - `epicenter` [lat, lon], `statco` [lat, lon]
- Or directly under `anEQ` the same fields as above

See `sample_data/README.md` for a full field list.

## Original vs Calculated

The middle panel shows paired values: Original (from metadata) vs Calculated (from signal and picks). Briefly:

- Original Distance (km)
  - Computed from file metadata epicenter and station coordinates using the haversine formula.

- Original Azimuth / Back Azimuth (deg)
  - Geographic azimuths between station and epicenter from metadata coordinates.
  - Azimuth = direction station → epicenter, Back Azimuth = epicenter → station.

- Calculated Azimuth (deg)
  - From the signal at the P pick using atan2 on the horizontal components:
    - az_calc = degrees(atan2(East, North)) mapped to [0, 360).
    - For consistency the GUI uses detrended values shown in the P Y boxes.

- Calculated Back Azimuth (deg)
  - back_az_calc = (az_calc + 180) mod 360.

- Event Azimuth (Calc.) (deg)
  - Applies P polarity correction to the calculated azimuth:
    - If P polarity = down (compression): event_az = az_calc
    - If P polarity = up (dilation): event_az = back_az_calc
    - If unknown: shown as “-”.

- S–P Distance (km)
  - Derived from S–P time and user-specified velocities Vp and Vs:
    - d = (Stime − Ptime) × (Vp × Vs) / (Vp − Vs)

Notes:
- Changing Vp/Vs immediately updates S–P Distance and the map epicenter prediction.
- Moving P/S picks updates the calculated azimuths and P polarity.

## Directory Structure

```
.
├── earthquake_signal_gui.py      # Main GUI
├── requirements.txt              # Dependencies
├── README.md                     # This page
├── LICENSE                       # MIT License
├── .gitignore                    # Git ignore rules
├── assets/
│   └── screenshot.png            # Screenshot (you added this)
├── sample_data/
│   └── README.md                 # Data format notes
└── (earthquake_analysis_results.xlsx)  # Output, ignored by git
```

## Troubleshooting

- Map not showing: verify `cartopy` install and system libs (GEOS/PROJ). On macOS:
  ```bash
  brew install geos proj
  ```
- No files listed: `.mat` files must be under `sample_data/` with the expected fields.
- Excel button produces no file: check write permissions and `openpyxl` is installed.

## License

MIT. See [LICENSE](LICENSE).
