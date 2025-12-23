# RS808 LV2

RS808 is an open-source, MIDI-triggered, 808-style drum synthesizer
implemented as an LV2 plugin.

It is a SuperCollider-accurate translation of classic 808-style
analog drum synthesis, featuring:

- Dual-mono main outputs (L = R)
- Individual mono outputs per instrument
- MIDI note triggering (no internal sequencer)
- Per-parameter MIDI CC mapping with learn mode
- Soft takeover (pickup) for MIDI CCs
- PipeWire / JACK / ALSA compatible via LV2 hosts

## Features

- Instruments: BD, SD, LT, MT, HT, CH, OH, Clap, Cowbell, Cymbal,
  Rimshot, Maracas, Claves
- One knob per synthesis parameter
- Per-instrument level controls
- Global accent and master level
- Configurable MIDI CC smoothing
- GTK3-based UI

## Outputs

- `main_l`, `main_r` (dual mono)
- Individual outputs for each instrument

Use the individual outputs for proper mixing.

## Build Dependencies

- C++17 compiler
- meson
- ninja
- lv2-dev
- gtk+-3.0-dev
- cairo-dev

## Building the plugin
On Ubuntu Ubuntu Studio:

Install build dependencies:
sudo apt update
sudo apt install -y build-essential meson ninja-build pkg-config \
  lv2-dev lilv-utils g++ \
  libc6-dev libgtk-3-dev libcairo2-dev


How to build:
cd rs808lv2
rm -rf build
meson setup build --buildtype=release --prefix=$HOME/.local
ninja -C build
ninja -C build install

Then the bundle ends up in:
~/.local/lib/lv2/rs808.lv2/

You can copy it to your favourite lv2 folder, mine happend to be ~/.lv2

Carla / Ardour

RS808 appears as an LV2 instrument plugin.

License

RS808 is licensed under the GNU General Public License v3.0 or later.

See the LICENSE file for details.

Trademark Notice

Roland and TR-808 are trademarks of Roland Corporation.
This project is not affiliated with or endorsed by Roland.

