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
#include <gtk/gtk.h>
#include <functional>
#include <string>

struct KnobMappingView {
  int cc = -1;       // -1 none
  int ch = 0;        // 0 omni, 1..16
  bool range_enabled = false;
  float range_min = 0.0f;
  float range_max = 1.0f;
  bool pickup = true;
  bool conflict = false;
};

struct Knob {
  GtkWidget* box=nullptr;
  GtkWidget* area=nullptr;
  GtkWidget* label=nullptr;
  GtkWidget* warn=nullptr;

  float value=0.0f, minv=0.0f, maxv=1.0f;
  bool dragging=false;
  double drag_y=0.0;
  float drag_start=0.0f;

  std::string param_name;
  KnobMappingView map;

  std::function<void(float)> on_change;

  // context actions (implemented by UI)
  std::function<void()> on_learn;
  std::function<void()> on_unlearn;
  std::function<void(int cc, int ch)> on_set_cc;
  std::function<void(bool en, float mn, float mx)> on_set_range;
  std::function<void(bool pickup)> on_set_pickup;

  Knob(const char* name, float minv_, float maxv_, float defv);
  GtkWidget* widget() const { return box; }

  void set_value(float v);
  void set_mapping(const KnobMappingView& m);

  void update_tooltip();
};
