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
#include <cstdint>
#include <vector>

struct RS808_State {
  static constexpr uint32_t kNoCC = 255;
  std::vector<uint8_t> cc_map;

  void init(size_t n_params) {
    cc_map.assign(n_params, (uint8_t)kNoCC);
  }
};
