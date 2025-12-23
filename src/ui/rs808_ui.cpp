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

#include <cstring>
#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <gtk/gtk.h>
#include <lv2/ui/ui.h>
#include <lv2/urid/urid.h>
#include <lv2/atom/atom.h>
#include <lv2/atom/forge.h>
#include <lv2/atom/util.h>
#include "../rs808_urid.hpp"
#include "knob.hpp"

#define RS808_UI_URI "https://example.org/lv2/rs808#ui"

enum PortIndex : uint32_t {
  P_MAIN_L=0, P_MAIN_R,
  P_BD, P_SD, P_CH, P_OH, P_CP, P_CB, P_CY, P_RIM, P_MA, P_CLV, P_LT, P_MT, P_HT,
  P_MIDI_IN, P_NOTIFY_OUT,
  P_MASTER_LEVEL, P_ACCENT_AMOUNT, P_CC_SMOOTH_MS,
  P_BD_LEVEL, P_BD_TUNE, P_BD_DECAY, P_BD_CLICK,
  P_SD_LEVEL, P_SD_TUNE, P_SD_DECAY, P_SD_TONE, P_SD_SNAPPY,
  P_LT_LEVEL, P_LT_TUNE, P_LT_DECAY,
  P_MT_LEVEL, P_MT_TUNE, P_MT_DECAY,
  P_HT_LEVEL, P_HT_TUNE, P_HT_DECAY,
  P_CH_LEVEL, P_CH_TONE, P_CH_DECAY,
  P_OH_LEVEL, P_OH_TONE, P_OH_DECAY,
  P_CP_LEVEL, P_CP_TONE, P_CP_DECAY,
  P_CB_LEVEL, P_CB_TUNE, P_CB_DECAY, P_CB_TONE,
  P_CY_LEVEL, P_CY_TONE, P_CY_DECAY,
  P_RIM_LEVEL, P_RIM_TONE,
  P_MA_LEVEL, P_MA_TONE,
  P_CLV_LEVEL
};

struct UI {
  LV2UI_Write_Function write;
  LV2UI_Controller controller;
  LV2_URID_Map* map=nullptr;
  RS808_URIs uris{};

  GtkWidget* root=nullptr;
  GtkWidget* status=nullptr;

  std::map<uint32_t, Knob*> knobs;
  std::vector<uint32_t> knob_ports; // same order as plugin (we rebuild locally)

  // current mapping cache for tooltip/warnings (key: knob index)
  std::vector<KnobMappingView> mapping;

  UI(LV2UI_Write_Function w, LV2UI_Controller c, LV2_URID_Map* m)
  : write(w), controller(c), map(m)
  {
    map_uris(map, &uris);

    auto add=[&](uint32_t p){ knob_ports.push_back(p); };
    add(P_MASTER_LEVEL);
    add(P_ACCENT_AMOUNT);

    add(P_BD_LEVEL); add(P_BD_TUNE); add(P_BD_DECAY); add(P_BD_CLICK);
    add(P_SD_LEVEL); add(P_SD_TUNE); add(P_SD_DECAY); add(P_SD_TONE); add(P_SD_SNAPPY);

    add(P_LT_LEVEL); add(P_LT_TUNE); add(P_LT_DECAY);
    add(P_MT_LEVEL); add(P_MT_TUNE); add(P_MT_DECAY);
    add(P_HT_LEVEL); add(P_HT_TUNE); add(P_HT_DECAY);

    add(P_CH_LEVEL); add(P_CH_TONE); add(P_CH_DECAY);
    add(P_OH_LEVEL); add(P_OH_TONE); add(P_OH_DECAY);

    add(P_CP_LEVEL); add(P_CP_TONE); add(P_CP_DECAY);

    add(P_CB_LEVEL); add(P_CB_TUNE); add(P_CB_DECAY); add(P_CB_TONE);

    add(P_CY_LEVEL); add(P_CY_TONE); add(P_CY_DECAY);

    add(P_RIM_LEVEL); add(P_RIM_TONE);
    add(P_MA_LEVEL); add(P_MA_TONE);
    add(P_CLV_LEVEL);

    mapping.assign(knob_ports.size(), KnobMappingView{});
    for(auto& mvv : mapping){
      mvv.cc = -1;
      mvv.ch = 0;
      mvv.range_enabled = false;
      mvv.range_min = 0.0f;
      mvv.range_max = 1.0f;
      mvv.pickup = true;
      mvv.conflict = false;
    }
  }

  void write_float(uint32_t port, float v){
    write(controller, port, sizeof(float), 0, &v);
  }

  int knob_index_for_port(uint32_t port){
    for(size_t i=0;i<knob_ports.size();++i) if(knob_ports[i]==port) return (int)i;
    return -1;
  }

  void set_status(const std::string& s){
    gtk_label_set_text(GTK_LABEL(status), s.c_str());
  }

  // We use “learn mode” by sending a special “arm learn” request via Notify port:
  // In this simplified version, plugin learns on next CC automatically when armed.
  // To avoid adding another Atom input port, we encode arm requests by writing negative values
  // into the smoothing port (harmless) is ugly — so we do it cleanly:
  // We'll use the existing NOTIFY port *from UI to plugin* (allowed by LV2 UI: we can write to atom ports).
  // Host must connect it; jalv does.
  void send_arm_learn(int knobIndex){
    // We'll send a LearnMsg object on the notify port index (P_NOTIFY_OUT) as UI->plugin event.
    // Many hosts allow this; jalv does. If a host doesn't, you can still “manual set”.
    // (If your host blocks UI->plugin atom writes, tell me; we’ll add a dedicated atom input port.)
    LV2_Atom_Forge f;
    lv2_atom_forge_init(&f, map);
    uint8_t buf[256];
    lv2_atom_forge_set_buffer(&f, buf, sizeof(buf));

    LV2_Atom_Forge_Frame seq;
    lv2_atom_forge_sequence_head(&f, &seq, 0);

    LV2_Atom_Forge_Frame obj;
    lv2_atom_forge_frame_time(&f, 0);
    lv2_atom_forge_object(&f, &obj, 0, uris.rs808_LearnMsg);
    lv2_atom_forge_key(&f, uris.rs808_paramIndex);
    lv2_atom_forge_int(&f, knobIndex);
    lv2_atom_forge_key(&f, uris.rs808_ccNumber);
    lv2_atom_forge_int(&f, -1); // special: arm learn
    lv2_atom_forge_pop(&f, &obj);
    lv2_atom_forge_pop(&f, &seq);

    write(controller, P_NOTIFY_OUT, lv2_atom_total_size((LV2_Atom*)buf), uris.atom_Sequence, buf);
    set_status("Learn: move a CC now...");
  }

  void compute_conflicts(){
    // duplicates allowed but show warning when same CC+ch used by >1 knobs (excluding None)
    for(auto& mvv : mapping) mvv.conflict = false;

    for(size_t i=0;i<mapping.size();++i){
      if(mapping[i].cc < 0) continue;
      for(size_t j=i+1;j<mapping.size();++j){
        if(mapping[j].cc < 0) continue;
        if(mapping[i].cc == mapping[j].cc && mapping[i].ch == mapping[j].ch){
          mapping[i].conflict = true;
          mapping[j].conflict = true;
        }
      }
    }

    for(size_t i=0;i<knob_ports.size();++i){
      auto it = knobs.find(knob_ports[i]);
      if(it!=knobs.end()) it->second->set_mapping(mapping[i]);
    }
  }

  GtkWidget* section(const char* title, const std::vector<Knob*>& ks, int cols=6){
    GtkWidget* frame=gtk_frame_new(title);
    GtkWidget* grid=gtk_grid_new();
    gtk_grid_set_row_spacing(GTK_GRID(grid), 4);
    gtk_grid_set_column_spacing(GTK_GRID(grid), 6);
    gtk_container_add(GTK_CONTAINER(frame), grid);

    int col=0,row=0;
    for(auto* k: ks){
      gtk_grid_attach(GTK_GRID(grid), k->widget(), col, row, 1, 1);
      col++;
      if(col>=cols){ col=0; row++; }
    }
    return frame;
  }

  Knob* add_knob(uint32_t port, const char* name, float mn, float mx, float defv){
    Knob* k = new Knob(name, mn, mx, defv);

    k->on_change = [this,port](float v){ write_float(port, v); };

    k->on_learn = [this,port](){
      int idx = knob_index_for_port(port);
      if(idx>=0) send_arm_learn(idx);
    };

      k->on_unlearn = [this, port, k]() {
        int idx = knob_index_for_port(port);
        if (idx >= 0) {
          mapping[idx].cc = -1;
          mapping[idx].ch = 0;
          mapping[idx].range_enabled = false;
          mapping[idx].pickup = true;
          compute_conflicts();
          set_status("Unlearned.");
          // We also send “manual set none”
          if (k->on_set_cc) {
            k->on_set_cc(-1, 0);
          }
        }
      };


      k->on_set_cc = [this,port](int cc, int ch){
        int idx = knob_index_for_port(port);
        if(idx<0) return;
        mapping[idx].cc = cc;
        mapping[idx].ch = ch;
        compute_conflicts();

        // Send mapping to plugin via notify atom message (cc packed with channel)
        LV2_Atom_Forge f;
        lv2_atom_forge_init(&f, map);
        uint8_t buf[256];
        lv2_atom_forge_set_buffer(&f, buf, sizeof(buf));

        LV2_Atom_Forge_Frame seq;
        lv2_atom_forge_sequence_head(&f, &seq, 0);

        LV2_Atom_Forge_Frame obj;
        lv2_atom_forge_frame_time(&f, 0);
        lv2_atom_forge_object(&f, &obj, 0, uris.rs808_LearnMsg);
        lv2_atom_forge_key(&f, uris.rs808_paramIndex);
        lv2_atom_forge_int(&f, idx);
        lv2_atom_forge_key(&f, uris.rs808_ccNumber);
        // pack: (ch<<8)|cc, with cc=-1 meaning none
        int packed = (cc < 0) ? -1 : ((ch<<8) | (cc & 127));
        lv2_atom_forge_int(&f, packed);
        lv2_atom_forge_pop(&f, &obj);
        lv2_atom_forge_pop(&f, &seq);

        write(controller, P_NOTIFY_OUT, lv2_atom_total_size((LV2_Atom*)buf), uris.atom_Sequence, buf);

        if(cc < 0) set_status("CC cleared.");
        else {
          std::string s = "Assigned CC" + std::to_string(cc) + (ch==0? " Omni" : (" Ch"+std::to_string(ch)));
          set_status(s);
        }
      };

      k->on_set_range = [this,port](bool en, float rmn, float rmx){
        int idx = knob_index_for_port(port);
        if(idx<0) return;
        mapping[idx].range_enabled = en;
        mapping[idx].range_min = rmn;
        mapping[idx].range_max = rmx;
        compute_conflicts();
        set_status(en ? "Range set." : "Range cleared.");

        // NOTE: range is UI-only in this minimal patch (tooltip + future extension).
        // If you want range enforced in DSP mapping, say so and I’ll add the atom message + DSP handling.
      };

      k->on_set_pickup = [this,port](bool pickup){
        int idx = knob_index_for_port(port);
        if(idx<0) return;
        mapping[idx].pickup = pickup;
        compute_conflicts();
        set_status(pickup ? "Pickup enabled." : "Pickup disabled.");
        // Same note as range: UI-only in this minimal patch.
      };

      knobs[port]=k;
      return k;
  }

  void build(){
    root = gtk_box_new(GTK_ORIENTATION_VERTICAL, 6);
    status = gtk_label_new("Right-click a knob → Learn / Unlearn / Set CC.");
    gtk_box_pack_start(GTK_BOX(root), status, FALSE, FALSE, 2);

    std::vector<Knob*> global;
    global.push_back(add_knob(P_MASTER_LEVEL,"Master Level",0,1,0.8f));
    global.push_back(add_knob(P_ACCENT_AMOUNT,"Accent Amount",0,2,1.0f));
    global.push_back(add_knob(P_CC_SMOOTH_MS,"CC Smooth (ms)",0,50,12.0f));

    std::vector<Knob*> bd;
    bd.push_back(add_knob(P_BD_LEVEL,"BD Level",0,1.5f,0.8f));
    bd.push_back(add_knob(P_BD_TUNE,"BD Tune",20,120,56.0f));
    bd.push_back(add_knob(P_BD_DECAY,"BD Decay",1,60,30.0f));
    bd.push_back(add_knob(P_BD_CLICK,"BD Attack/Click",0,1,0.5f));

    std::vector<Knob*> sd;
    sd.push_back(add_knob(P_SD_LEVEL,"SD Level",0,1.5f,0.7f));
    sd.push_back(add_knob(P_SD_TUNE,"SD Tune (semi)",-12,12,0.0f));
    sd.push_back(add_knob(P_SD_DECAY,"SD Decay",0.2f,3.0f,1.0f));
    sd.push_back(add_knob(P_SD_TONE,"SD Tone",80,1200,340.0f));
    sd.push_back(add_knob(P_SD_SNAPPY,"SD Snappy",0,1,0.3f));

    std::vector<Knob*> toms;
    toms.push_back(add_knob(P_LT_LEVEL,"LT Level",0,1.5f,0.6f));
    toms.push_back(add_knob(P_LT_TUNE,"LT Tune",40,200,80.0f));
    toms.push_back(add_knob(P_LT_DECAY,"LT Decay",1,30,20.0f));

    toms.push_back(add_knob(P_MT_LEVEL,"MT Level",0,1.5f,0.6f));
    toms.push_back(add_knob(P_MT_TUNE,"MT Tune",60,260,120.0f));
    toms.push_back(add_knob(P_MT_DECAY,"MT Decay",1,25,15.0f));

    toms.push_back(add_knob(P_HT_LEVEL,"HT Level",0,1.5f,0.6f));
    toms.push_back(add_knob(P_HT_TUNE,"HT Tune",80,400,165.0f));
    toms.push_back(add_knob(P_HT_DECAY,"HT Decay",1,20,10.0f));

    std::vector<Knob*> hats;
    hats.push_back(add_knob(P_CH_LEVEL,"CH Level",0,1.5f,0.55f));
    hats.push_back(add_knob(P_CH_TONE,"CH Tone",0,1,0.5f));
    hats.push_back(add_knob(P_CH_DECAY,"CH Decay",0.02f,2.0f,0.42f));

    hats.push_back(add_knob(P_OH_LEVEL,"OH Level",0,1.5f,0.65f));
    hats.push_back(add_knob(P_OH_TONE,"OH Tone",0,1,0.5f));
    hats.push_back(add_knob(P_OH_DECAY,"OH Decay",0.05f,8.0f,0.5f));

    std::vector<Knob*> clap;
    clap.push_back(add_knob(P_CP_LEVEL,"Clap Level",0,1.5f,0.6f));
    clap.push_back(add_knob(P_CP_TONE,"Clap Tone",0,1,0.5f));
    clap.push_back(add_knob(P_CP_DECAY,"Clap Decay",0.1f,3.0f,0.5f));

    std::vector<Knob*> cow;
    cow.push_back(add_knob(P_CB_LEVEL,"Cowbell Level",0,1.5f,0.55f));
    cow.push_back(add_knob(P_CB_TUNE,"Cowbell Tune",0.5f,2.0f,1.0f));
    cow.push_back(add_knob(P_CB_DECAY,"Cowbell Decay",0.2f,3.0f,1.0f));
    cow.push_back(add_knob(P_CB_TONE,"Cowbell Tone",0,1,0.5f));

    std::vector<Knob*> cym;
    cym.push_back(add_knob(P_CY_LEVEL,"Cymbal Level",0,1.5f,0.55f));
    cym.push_back(add_knob(P_CY_TONE,"Cymbal Tone",0,1,0.5f));
    cym.push_back(add_knob(P_CY_DECAY,"Cymbal Decay",0.2f,8.0f,2.0f));

    std::vector<Knob*> perc;
    perc.push_back(add_knob(P_RIM_LEVEL,"Rimshot Level",0,1.5f,0.55f));
    perc.push_back(add_knob(P_RIM_TONE,"Rimshot Tone",0,1,0.5f));
    perc.push_back(add_knob(P_MA_LEVEL,"Maracas Level",0,1.5f,0.5f));
    perc.push_back(add_knob(P_MA_TONE,"Maracas Tone",0,1,0.5f));
    perc.push_back(add_knob(P_CLV_LEVEL,"Claves Level",0,1.5f,0.5f));

    GtkWidget* grid=gtk_grid_new();
    gtk_grid_set_row_spacing(GTK_GRID(grid), 6);
    gtk_grid_set_column_spacing(GTK_GRID(grid), 6);

    gtk_grid_attach(GTK_GRID(grid), section("Global", global, 3), 0,0,2,1);
    gtk_grid_attach(GTK_GRID(grid), section("BD", bd, 4), 0,1,2,1);
    gtk_grid_attach(GTK_GRID(grid), section("SD", sd, 5), 0,2,2,1);
    gtk_grid_attach(GTK_GRID(grid), section("Toms (LT/MT/HT)", toms, 6), 0,3,2,1);
    gtk_grid_attach(GTK_GRID(grid), section("Hats (CH/OH)", hats, 6), 0,4,2,1);
    gtk_grid_attach(GTK_GRID(grid), section("Clap", clap, 3), 0,5,1,1);
    gtk_grid_attach(GTK_GRID(grid), section("Cowbell", cow, 4), 1,5,1,1);
    gtk_grid_attach(GTK_GRID(grid), section("Cymbal", cym, 3), 0,6,1,1);
    gtk_grid_attach(GTK_GRID(grid), section("Rim/Maracas/Claves", perc, 5), 1,6,1,1);

    gtk_box_pack_start(GTK_BOX(root), grid, TRUE, TRUE, 2);
    gtk_widget_show_all(root);

    compute_conflicts();
  }

  void port_event(uint32_t port, uint32_t size, uint32_t format, const void* buffer){
    (void)format;
    // Plugin -> UI notify (learn assigned)
    if(port == P_NOTIFY_OUT && buffer && size >= sizeof(LV2_Atom)){
      const LV2_Atom* atom=(const LV2_Atom*)buffer;
      if(atom->type == uris.atom_Sequence){
        const LV2_Atom_Sequence* seq=(const LV2_Atom_Sequence*)buffer;
        LV2_ATOM_SEQUENCE_FOREACH(seq, ev){
          const LV2_Atom* a=&ev->body;
          if(a->type == uris.rs808_LearnMsg){
            const LV2_Atom_Object* obj=(const LV2_Atom_Object*)a;
            int pidx=-1, packed=-1;
            LV2_ATOM_OBJECT_FOREACH(obj, prop){
              if(prop->key == uris.rs808_paramIndex && prop->value.type == uris.atom_Int)
                pidx=((const LV2_Atom_Int*)&prop->value)->body;
              else if(prop->key == uris.rs808_ccNumber && prop->value.type == uris.atom_Int)
                packed=((const LV2_Atom_Int*)&prop->value)->body;
            }
            if(pidx>=0 && pidx < (int)mapping.size() && packed>=0){
              int cc = packed & 127;
              int ch = (packed >> 8) & 255;
              mapping[pidx].cc = cc;
              mapping[pidx].ch = ch;
              compute_conflicts();
              std::string s = "Assigned CC" + std::to_string(cc) + " Ch" + std::to_string(ch);
              set_status(s);
            }
          }
        }
      }
      return;
    }

    // regular param updates from host automation
    auto it = knobs.find(port);
    if(it!=knobs.end() && buffer && size==sizeof(float)){
      it->second->set_value(*(const float*)buffer);
    }
  }
};

static LV2UI_Handle
ui_instantiate(const LV2UI_Descriptor*, const char*, const char*,
               LV2UI_Write_Function write_function,
               LV2UI_Controller controller,
               LV2UI_Widget* widget,
               const LV2_Feature* const* features)
{
  LV2_URID_Map* map=nullptr;
  for(int i=0; features && features[i]; ++i){
    if(!std::strcmp(features[i]->URI, LV2_URID__map))
      map=(LV2_URID_Map*)features[i]->data;
  }
  if(!map) return nullptr;

  UI* ui=new UI(write_function, controller, map);
  ui->build();
  *widget = ui->root;
  return (LV2UI_Handle)ui;
}

static void ui_cleanup(LV2UI_Handle handle){ delete (UI*)handle; }

static void ui_port_event(LV2UI_Handle handle, uint32_t port_index,
                          uint32_t buffer_size, uint32_t format, const void* buffer){
  ((UI*)handle)->port_event(port_index, buffer_size, format, buffer);
                          }

                          static const void* ui_extension_data(const char*){ return nullptr; }

                          static const LV2UI_Descriptor ui_descriptor = {
                            RS808_UI_URI,
                            ui_instantiate,
                            ui_cleanup,
                            ui_port_event,
                            ui_extension_data
                          };

                          extern "C" LV2_SYMBOL_EXPORT
                          const LV2UI_Descriptor* lv2ui_descriptor(uint32_t index){
                            return index==0 ? &ui_descriptor : nullptr;
                          }
