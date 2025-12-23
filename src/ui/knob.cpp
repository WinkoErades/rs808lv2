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
#include "knob.hpp"
#include <cairo.h>
#include <cmath>
#include <algorithm>

static gboolean on_draw(GtkWidget* w, cairo_t* cr, gpointer user){
  Knob* k=(Knob*)user;
  GtkAllocation a; gtk_widget_get_allocation(w, &a);
  double cx=a.width/2.0, cy=a.height/2.0;
  double r=std::min(a.width,a.height)*0.38;

  cairo_set_source_rgb(cr, 0.12,0.12,0.12);
  cairo_paint(cr);

  // knob body
  cairo_set_source_rgb(cr, 0.22,0.22,0.22);
  cairo_arc(cr, cx, cy, r, 0, 2*M_PI);
  cairo_fill(cr);

  // pointer
  float t = (k->value - k->minv) / (k->maxv - k->minv);
  t = std::clamp(t, 0.0f, 1.0f);

  const double ang0 = 3.0*M_PI/4.0;
  const double ang1 = 9.0*M_PI/4.0;
  const double ang = ang0 + (ang1-ang0)*t;

  cairo_set_line_width(cr, 3.0);
  cairo_set_source_rgb(cr, 0.85,0.85,0.85);
  cairo_move_to(cr, cx, cy);
  cairo_line_to(cr, cx + std::cos(ang)*r*0.85, cy + std::sin(ang)*r*0.85);
  cairo_stroke(cr);

  return FALSE;
}

static gboolean on_button(GtkWidget*, GdkEventButton* ev, gpointer user){
  Knob* k=(Knob*)user;
  if(ev->button==1){
    k->dragging=true;
    k->drag_y=ev->y_root;
    k->drag_start=k->value;
    return TRUE;
  }
  if(ev->button==3){
    // right-click menu
    GtkWidget* menu = gtk_menu_new();

    auto add_item=[&](const char* txt, std::function<void()> fn){
      GtkWidget* mi = gtk_menu_item_new_with_label(txt);
      g_signal_connect(G_OBJECT(mi), "activate", G_CALLBACK(+[](GtkMenuItem*, gpointer u){
        auto* f = (std::function<void()>*)u;
        (*f)();
        delete f;
      }), new std::function<void()>(std::move(fn)));
      gtk_menu_shell_append(GTK_MENU_SHELL(menu), mi);
    };

    add_item("Learn",   [k](){ if(k->on_learn) k->on_learn(); });
    add_item("Unlearn", [k](){ if(k->on_unlearn) k->on_unlearn(); });

    add_item("Set CC manually...", [k](){
      GtkWidget* d = gtk_dialog_new_with_buttons("Set CC",
                                                 nullptr, GTK_DIALOG_MODAL,
                                                 "_Cancel", GTK_RESPONSE_CANCEL,
                                                 "_OK", GTK_RESPONSE_OK,
                                                 nullptr);

      GtkWidget* box = gtk_dialog_get_content_area(GTK_DIALOG(d));

      GtkWidget* grid = gtk_grid_new();
      gtk_grid_set_row_spacing(GTK_GRID(grid), 6);
      gtk_grid_set_column_spacing(GTK_GRID(grid), 6);

      GtkWidget* l1 = gtk_label_new("CC (0-127, -1=None):");
      GtkWidget* cc_spin = gtk_spin_button_new_with_range(-1,127,1);
      gtk_spin_button_set_value(GTK_SPIN_BUTTON(cc_spin), k->map.cc);

      GtkWidget* l2 = gtk_label_new("Channel (0=Omni, 1-16):");
      GtkWidget* ch_spin = gtk_spin_button_new_with_range(0,16,1);
      gtk_spin_button_set_value(GTK_SPIN_BUTTON(ch_spin), k->map.ch);

      gtk_grid_attach(GTK_GRID(grid), l1, 0,0,1,1);
      gtk_grid_attach(GTK_GRID(grid), cc_spin, 1,0,1,1);
      gtk_grid_attach(GTK_GRID(grid), l2, 0,1,1,1);
      gtk_grid_attach(GTK_GRID(grid), ch_spin, 1,1,1,1);

      gtk_container_add(GTK_CONTAINER(box), grid);
      gtk_widget_show_all(d);

      if(gtk_dialog_run(GTK_DIALOG(d)) == GTK_RESPONSE_OK){
        int cc = (int)gtk_spin_button_get_value(GTK_SPIN_BUTTON(cc_spin));
        int ch = (int)gtk_spin_button_get_value(GTK_SPIN_BUTTON(ch_spin));
        if(k->on_set_cc) k->on_set_cc(cc, ch);
      }
      gtk_widget_destroy(d);
    });

    add_item("Set Range...", [k](){
      GtkWidget* d = gtk_dialog_new_with_buttons("Range",
                                                 nullptr, GTK_DIALOG_MODAL,
                                                 "_Cancel", GTK_RESPONSE_CANCEL,
                                                 "_OK", GTK_RESPONSE_OK,
                                                 nullptr);

      GtkWidget* box = gtk_dialog_get_content_area(GTK_DIALOG(d));
      GtkWidget* grid = gtk_grid_new();
      gtk_grid_set_row_spacing(GTK_GRID(grid), 6);
      gtk_grid_set_column_spacing(GTK_GRID(grid), 6);

      GtkWidget* en = gtk_check_button_new_with_label("Enable custom range");
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(en), k->map.range_enabled);

      GtkWidget* mn = gtk_spin_button_new_with_range(k->minv, k->maxv, (k->maxv-k->minv)/500.0);
      GtkWidget* mx = gtk_spin_button_new_with_range(k->minv, k->maxv, (k->maxv-k->minv)/500.0);
      gtk_spin_button_set_value(GTK_SPIN_BUTTON(mn), k->map.range_min);
      gtk_spin_button_set_value(GTK_SPIN_BUTTON(mx), k->map.range_max);

      gtk_grid_attach(GTK_GRID(grid), en, 0,0,2,1);
      gtk_grid_attach(GTK_GRID(grid), gtk_label_new("Min:"), 0,1,1,1);
      gtk_grid_attach(GTK_GRID(grid), mn, 1,1,1,1);
      gtk_grid_attach(GTK_GRID(grid), gtk_label_new("Max:"), 0,2,1,1);
      gtk_grid_attach(GTK_GRID(grid), mx, 1,2,1,1);

      gtk_container_add(GTK_CONTAINER(box), grid);
      gtk_widget_show_all(d);

      if(gtk_dialog_run(GTK_DIALOG(d)) == GTK_RESPONSE_OK){
        bool enabled = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(en));
        float rmn = (float)gtk_spin_button_get_value(GTK_SPIN_BUTTON(mn));
        float rmx = (float)gtk_spin_button_get_value(GTK_SPIN_BUTTON(mx));
        if(k->on_set_range) k->on_set_range(enabled, rmn, rmx);
      }
      gtk_widget_destroy(d);
    });

    add_item("Pickup (Soft Takeover)", [k](){
      bool newv = !k->map.pickup;
      if(k->on_set_pickup) k->on_set_pickup(newv);
    });

      gtk_widget_show_all(menu);
      gtk_menu_popup_at_pointer(GTK_MENU(menu), (GdkEvent*)ev);
      return TRUE;
  }
  return FALSE;
}

static gboolean on_release(GtkWidget*, GdkEventButton* ev, gpointer user){
  Knob* k=(Knob*)user;
  if(ev->button==1){ k->dragging=false; return TRUE; }
  return FALSE;
}

static gboolean on_motion(GtkWidget*, GdkEventMotion* ev, gpointer user){
  Knob* k=(Knob*)user;
  if(!k->dragging) return FALSE;
  double dy = (k->drag_y - ev->y_root);
  float range = (k->maxv - k->minv);
  float delta = float(dy) * (range / 240.0f);
  float v = k->drag_start + delta;
  k->set_value(v);
  if(k->on_change) k->on_change(k->value);
  return TRUE;
}

Knob::Knob(const char* name, float minv_, float maxv_, float defv)
: minv(minv_), maxv(maxv_), param_name(name)
{
  box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 2);

  GtkWidget* top = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 2);
  label = gtk_label_new(name);
  gtk_widget_set_halign(label, GTK_ALIGN_START);

  warn = gtk_label_new("⚠");
  gtk_widget_set_halign(warn, GTK_ALIGN_END);
  gtk_widget_set_opacity(warn, 0.0); // hidden by default

  gtk_box_pack_start(GTK_BOX(top), label, TRUE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(top), warn, FALSE, FALSE, 0);

  area = gtk_drawing_area_new();
  gtk_widget_set_size_request(area, 56, 56);

  gtk_box_pack_start(GTK_BOX(box), top, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(box), area, FALSE, FALSE, 0);

  gtk_widget_add_events(area, GDK_BUTTON_PRESS_MASK|GDK_BUTTON_RELEASE_MASK|GDK_POINTER_MOTION_MASK);
  g_signal_connect(G_OBJECT(area), "draw", G_CALLBACK(on_draw), this);
  g_signal_connect(G_OBJECT(area), "button-press-event", G_CALLBACK(on_button), this);
  g_signal_connect(G_OBJECT(area), "button-release-event", G_CALLBACK(on_release), this);
  g_signal_connect(G_OBJECT(area), "motion-notify-event", G_CALLBACK(on_motion), this);

  set_value(defv);
  update_tooltip();
}

void Knob::set_value(float v){
  value = std::clamp(v, minv, maxv);
  gtk_widget_queue_draw(area);
  update_tooltip();
}

void Knob::set_mapping(const KnobMappingView& m){
  map = m;
  gtk_widget_set_opacity(warn, map.conflict ? 1.0 : 0.0);
  update_tooltip();
}

void Knob::update_tooltip(){
  std::string ccstr = "None";
  if(map.cc >= 0){
    if(map.ch == 0) ccstr = "CC" + std::to_string(map.cc) + " Omni";
    else ccstr = "CC" + std::to_string(map.cc) + " Ch" + std::to_string(map.ch);
  }

  std::string rstr = "Normal";
  if(map.range_enabled){
    rstr = "Range " + std::to_string(map.range_min) + " .. " + std::to_string(map.range_max);
  }

  std::string pstr = map.pickup ? "Pickup: On" : "Pickup: Off";

  char buf[512];
  std::snprintf(buf, sizeof(buf),
                "%s\nValue: %.4f\nCC: %s\n%s\n%s",
                param_name.c_str(), value, ccstr.c_str(), rstr.c_str(), pstr.c_str());

  gtk_widget_set_tooltip_text(area, buf);
}
