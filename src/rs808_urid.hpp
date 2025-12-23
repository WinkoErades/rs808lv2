/*
 * RS808 LV2 Drum Synth
 * Copyright (C) 2025 Winko Erades van den Berg
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once
#include <lv2/urid/urid.h>
#include <lv2/atom/atom.h>
#include <lv2/midi/midi.h>

struct RS808_URIs {
  LV2_URID atom_Sequence;
  LV2_URID atom_Float;
  LV2_URID atom_Int;
  LV2_URID midi_MidiEvent;

  LV2_URID rs808_LearnMsg;
  LV2_URID rs808_paramIndex;
  LV2_URID rs808_ccNumber;
};

inline void map_uris(LV2_URID_Map* map, RS808_URIs* uris) {
  uris->atom_Sequence  = map->map(map->handle, LV2_ATOM__Sequence);
  uris->atom_Float     = map->map(map->handle, LV2_ATOM__Float);
  uris->atom_Int       = map->map(map->handle, LV2_ATOM__Int);
  uris->midi_MidiEvent = map->map(map->handle, LV2_MIDI__MidiEvent);

  uris->rs808_LearnMsg   = map->map(map->handle, "https://example.org/lv2/rs808#LearnMsg");
  uris->rs808_paramIndex = map->map(map->handle, "https://example.org/lv2/rs808#paramIndex");
  uris->rs808_ccNumber   = map->map(map->handle, "https://example.org/lv2/rs808#ccNumber");
}
