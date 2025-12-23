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
#include <vector>
#include <string>
#include <lv2/core/lv2.h>
#include <lv2/atom/atom.h>
#include <lv2/atom/forge.h>
#include <lv2/state/state.h>
#include <lv2/urid/urid.h>
#include <lv2/midi/midi.h>

#include "rs808_sc_dsp.hpp"
#include "rs808_urid.hpp"

#define RS808_URI "https://example.org/lv2/rs808"

enum PortIndex : uint32_t {
  // Audio outs (dual mono main + per instrument)
  P_MAIN_L=0, P_MAIN_R,
  P_BD, P_SD, P_CH, P_OH, P_CP, P_CB, P_CY, P_RIM, P_MA, P_CLV, P_LT, P_MT, P_HT,

  // Atom
  P_MIDI_IN,
  P_NOTIFY_OUT,

  // Global controls
  P_MASTER_LEVEL,
  P_ACCENT_AMOUNT,
  P_CC_SMOOTH_MS,

  // BD
  P_BD_LEVEL, P_BD_TUNE, P_BD_DECAY, P_BD_CLICK,

  // SD
  P_SD_LEVEL, P_SD_TUNE, P_SD_DECAY, P_SD_TONE, P_SD_SNAPPY,

  // LT/MT/HT
  P_LT_LEVEL, P_LT_TUNE, P_LT_DECAY,
  P_MT_LEVEL, P_MT_TUNE, P_MT_DECAY,
  P_HT_LEVEL, P_HT_TUNE, P_HT_DECAY,

  // CH/OH
  P_CH_LEVEL, P_CH_TONE, P_CH_DECAY,
  P_OH_LEVEL, P_OH_TONE, P_OH_DECAY,

  // Clap
  P_CP_LEVEL, P_CP_TONE, P_CP_DECAY,

  // Cowbell
  P_CB_LEVEL, P_CB_TUNE, P_CB_DECAY, P_CB_TONE,

  // Cymbal
  P_CY_LEVEL, P_CY_TONE, P_CY_DECAY,

  // Rimshot
  P_RIM_LEVEL, P_RIM_TONE,

  // Maracas
  P_MA_LEVEL, P_MA_TONE,

  // Claves
  P_CLV_LEVEL,

  P_COUNT
};

static inline float cc_to_01(uint8_t cc) { return float(cc) / 127.0f; }

enum class ChanMode : uint8_t {
  Omni = 0,
  Ch1  = 1, // ... Ch16 = 16
};

struct CCMapEntry {
  uint8_t cc = 255;          // 0..127 or 255 none
  uint8_t chan = 0;          // 0 omni, 1..16 specific
  float scale_min = 0.0f;    // if scale_enabled==0, ignored
  float scale_max = 1.0f;
  uint8_t scale_enabled = 0; // 0 normal, 1 custom
  uint8_t pickup = 1;        // soft takeover
  // internal pickup latch
  uint8_t picked_up = 0;
};

struct SmoothValue {
  float y=0.0f;
  float a=0.0f; // 0..1 (closer to 0 = faster)
  void set_time_ms(float sr, float ms){
    ms = clampf(ms, 0.0f, 200.0f);
    float t = ms * 0.001f;
    if(t <= 0.0f){
      a = 0.0f;
    } else {
      a = std::exp(-1.0f / (t * sr));
    }
  }
  float process(float x){
    y = (1.0f - a) * x + a * y;
    return y;
  }
};

struct RS808Plugin {
  // Audio
  float* out_main_l=nullptr; float* out_main_r=nullptr;
  float* out_bd=nullptr; float* out_sd=nullptr; float* out_ch=nullptr; float* out_oh=nullptr;
  float* out_cp=nullptr; float* out_cb=nullptr; float* out_cy=nullptr;
  float* out_rim=nullptr; float* out_ma=nullptr; float* out_clv=nullptr;
  float* out_lt=nullptr; float* out_mt=nullptr; float* out_ht=nullptr;

  const LV2_Atom_Sequence* midi_in=nullptr;
  LV2_Atom_Sequence* notify_out=nullptr;

  const float* ctrl[P_COUNT]{};

  LV2_URID_Map* map=nullptr;
  RS808_URIs uris{};
  LV2_Atom_Forge forge{};
  LV2_Atom_Forge_Frame notify_frame{};

  RS808 dsp{};
  RS808_Params params{};

  // Mapping for each continuous port we expose as a knob
  std::vector<uint32_t> knob_ports;
  std::vector<CCMapEntry> ccmap;
  std::vector<SmoothValue> smooth;

  int learn_knob_index = -1;

  double sr=48000.0;

  RS808Plugin(double rate, LV2_URID_Map* m)
  : map(m), sr(rate)
  {
    dsp.setSR((float)rate);
    map_uris(map, &uris);
    lv2_atom_forge_init(&forge, map);

    auto add_knob=[&](uint32_t p){ knob_ports.push_back(p); };

    add_knob(P_MASTER_LEVEL);
    add_knob(P_ACCENT_AMOUNT);

    add_knob(P_BD_LEVEL); add_knob(P_BD_TUNE); add_knob(P_BD_DECAY); add_knob(P_BD_CLICK);

    add_knob(P_SD_LEVEL); add_knob(P_SD_TUNE); add_knob(P_SD_DECAY); add_knob(P_SD_TONE); add_knob(P_SD_SNAPPY);

    add_knob(P_LT_LEVEL); add_knob(P_LT_TUNE); add_knob(P_LT_DECAY);
    add_knob(P_MT_LEVEL); add_knob(P_MT_TUNE); add_knob(P_MT_DECAY);
    add_knob(P_HT_LEVEL); add_knob(P_HT_TUNE); add_knob(P_HT_DECAY);

    add_knob(P_CH_LEVEL); add_knob(P_CH_TONE); add_knob(P_CH_DECAY);
    add_knob(P_OH_LEVEL); add_knob(P_OH_TONE); add_knob(P_OH_DECAY);

    add_knob(P_CP_LEVEL); add_knob(P_CP_TONE); add_knob(P_CP_DECAY);

    add_knob(P_CB_LEVEL); add_knob(P_CB_TUNE); add_knob(P_CB_DECAY); add_knob(P_CB_TONE);

    add_knob(P_CY_LEVEL); add_knob(P_CY_TONE); add_knob(P_CY_DECAY);

    add_knob(P_RIM_LEVEL); add_knob(P_RIM_TONE);
    add_knob(P_MA_LEVEL); add_knob(P_MA_TONE);
    add_knob(P_CLV_LEVEL);

    ccmap.assign(knob_ports.size(), CCMapEntry{});
    smooth.assign(knob_ports.size(), SmoothValue{});
    for(auto& s : smooth) s.y = 0.0f;
  }

  static Inst note_to_inst(int note){
    switch(note){
      case 36: return I_BD;
      case 38: return I_SD;
      case 39: return I_CP;
      case 37: return I_RIM;
      case 70: return I_MA;
      case 75: return I_CLV;
      case 56: return I_CB;
      case 42: return I_CH;
      case 46: return I_OH;
      case 49: return I_CY;
      case 43: return I_LT;
      case 47: return I_MT;
      case 50: return I_HT;
      default: return I_COUNT;
    }
  }

  float getf(uint32_t p, float defv){
    return (ctrl[p] ? *ctrl[p] : defv);
  }

  static void port_range(uint32_t port, float& mn, float& mx){
    // Ranges match the “knob feel”, not raw SC internals
    mn=0.0f; mx=1.0f;
    switch(port){
      case P_MASTER_LEVEL: mn=0; mx=1; break;
      case P_ACCENT_AMOUNT: mn=0; mx=2; break;

      case P_CC_SMOOTH_MS: mn=0; mx=50; break;

      case P_BD_LEVEL: mn=0; mx=1.5f; break;
      case P_BD_TUNE:  mn=20; mx=120; break;
      case P_BD_DECAY: mn=1; mx=60; break;
      case P_BD_CLICK: mn=0; mx=1; break;

      case P_SD_LEVEL: mn=0; mx=1.5f; break;
      case P_SD_TUNE:  mn=-12; mx=12; break;       // semitone shift
      case P_SD_DECAY: mn=0.2f; mx=3.0f; break;    // multiplier
      case P_SD_TONE:  mn=80; mx=1200; break;
      case P_SD_SNAPPY:mn=0; mx=1; break;

      case P_LT_LEVEL: case P_MT_LEVEL: case P_HT_LEVEL:
      case P_CH_LEVEL: case P_OH_LEVEL:
      case P_CP_LEVEL:
      case P_CB_LEVEL:
      case P_CY_LEVEL:
      case P_RIM_LEVEL:
      case P_MA_LEVEL:
      case P_CLV_LEVEL:
        mn=0; mx=1.5f; break;

      case P_LT_TUNE: mn=40; mx=200; break;
      case P_MT_TUNE: mn=60; mx=260; break;
      case P_HT_TUNE: mn=80; mx=400; break;

      case P_LT_DECAY: mn=1; mx=30; break;
      case P_MT_DECAY: mn=1; mx=25; break;
      case P_HT_DECAY: mn=1; mx=20; break;

      case P_CH_TONE: case P_OH_TONE: case P_CP_TONE:
      case P_CB_TONE: case P_CY_TONE:
      case P_RIM_TONE: case P_MA_TONE:
        mn=0; mx=1; break;

      case P_CH_DECAY: mn=0.02f; mx=2.0f; break;
      case P_OH_DECAY: mn=0.05f; mx=8.0f; break;
      case P_CP_DECAY: mn=0.1f; mx=3.0f; break;
      case P_CB_TUNE:  mn=0.5f; mx=2.0f; break;
      case P_CB_DECAY: mn=0.2f; mx=3.0f; break;

      case P_CY_DECAY: mn=0.2f; mx=8.0f; break;

      default: break;
    }
  }

  int knob_index_for_port(uint32_t port){
    for(size_t i=0;i<knob_ports.size();++i) if(knob_ports[i]==port) return (int)i;
    return -1;
  }

  void send_notify_learn(int knobIndex, int cc, int ch){
    if(!notify_out) return;

    lv2_atom_forge_set_buffer(&forge, (uint8_t*)notify_out, 4096);
    lv2_atom_forge_sequence_head(&forge, &notify_frame, 0);

    LV2_Atom_Forge_Frame obj;
    lv2_atom_forge_frame_time(&forge, 0);
    lv2_atom_forge_object(&forge, &obj, 0, uris.rs808_LearnMsg);

    lv2_atom_forge_key(&forge, uris.rs808_paramIndex);
    lv2_atom_forge_int(&forge, knobIndex);

    lv2_atom_forge_key(&forge, uris.rs808_ccNumber);
    // pack cc + channel into one int for simplicity: (ch<<8)|cc
    lv2_atom_forge_int(&forge, (ch<<8) | cc);

    lv2_atom_forge_pop(&forge, &obj);
    lv2_atom_forge_pop(&forge, &notify_frame);
  }

  void update_params_from_ports(){
    params.master = getf(P_MASTER_LEVEL, 0.8f);
    params.accent = getf(P_ACCENT_AMOUNT, 1.0f);

    params.bd_level = getf(P_BD_LEVEL, 0.8f);
    params.bd_tune  = getf(P_BD_TUNE, 56.0f);
    params.bd_decay = getf(P_BD_DECAY, 30.0f);
    params.bd_click = getf(P_BD_CLICK, 0.5f);

    params.sd_level  = getf(P_SD_LEVEL, 0.7f);
    params.sd_tune   = getf(P_SD_TUNE, 0.0f);
    params.sd_decay  = getf(P_SD_DECAY, 1.0f);
    params.sd_tone   = getf(P_SD_TONE, 340.0f);
    params.sd_snappy = getf(P_SD_SNAPPY, 0.3f);

    params.lt_level = getf(P_LT_LEVEL, 0.6f);
    params.lt_tune  = getf(P_LT_TUNE, 80.0f);
    params.lt_decay = getf(P_LT_DECAY, 20.0f);

    params.mt_level = getf(P_MT_LEVEL, 0.6f);
    params.mt_tune  = getf(P_MT_TUNE, 120.0f);
    params.mt_decay = getf(P_MT_DECAY, 15.0f);

    params.ht_level = getf(P_HT_LEVEL, 0.6f);
    params.ht_tune  = getf(P_HT_TUNE, 165.0f);
    params.ht_decay = getf(P_HT_DECAY, 10.0f);

    params.ch_level = getf(P_CH_LEVEL, 0.55f);
    params.ch_tone  = getf(P_CH_TONE, 0.5f);
    params.ch_decay = getf(P_CH_DECAY, 0.42f);

    params.oh_level = getf(P_OH_LEVEL, 0.65f);
    params.oh_tone  = getf(P_OH_TONE, 0.5f);
    params.oh_decay = getf(P_OH_DECAY, 0.5f);

    params.cp_level = getf(P_CP_LEVEL, 0.6f);
    params.cp_tone  = getf(P_CP_TONE, 0.5f);
    params.cp_decay = getf(P_CP_DECAY, 0.5f);

    params.cb_level = getf(P_CB_LEVEL, 0.55f);
    params.cb_tune  = getf(P_CB_TUNE, 1.0f);
    params.cb_decay = getf(P_CB_DECAY, 1.0f);
    params.cb_tone  = getf(P_CB_TONE, 0.5f);

    params.cy_level = getf(P_CY_LEVEL, 0.55f);
    params.cy_tone  = getf(P_CY_TONE, 0.5f);
    params.cy_decay = getf(P_CY_DECAY, 2.0f);

    params.rim_level= getf(P_RIM_LEVEL, 0.55f);
    params.rim_tone = getf(P_RIM_TONE, 0.5f);

    params.ma_level = getf(P_MA_LEVEL, 0.5f);
    params.ma_tone  = getf(P_MA_TONE, 0.5f);

    params.clv_level= getf(P_CLV_LEVEL, 0.5f);

    // smoothing time for CC updates
    float ms = getf(P_CC_SMOOTH_MS, 12.0f);
    for(auto& s : smooth) s.set_time_ms((float)sr, ms);
  }

  // Apply CC to mapped knobs with pickup + smoothing.
  // Host ports remain authoritative when automated; CC acts like “performance control”.
  void apply_cc(uint8_t status, uint8_t cc, uint8_t val){
    uint8_t ch = (status & 0x0F) + 1; // 1..16
    float v01 = cc_to_01(val);

    for(size_t i=0;i<ccmap.size();++i){
      auto& m = ccmap[i];
      if(m.cc == 255 || m.cc != cc) continue;
      if(m.chan != 0 && m.chan != ch) continue;

      uint32_t port = knob_ports[i];
      float mn, mx; port_range(port, mn, mx);

      float target = mn + v01*(mx-mn);
      if(m.scale_enabled){
        float smn = m.scale_min;
        float smx = m.scale_max;
        // scale_min/max are in “port units”
        target = clampf(target, std::min(smn,smx), std::max(smn,smx));
      }

      // pickup: wait until controller crosses current value
      float current = get_port_value(port);
      if(m.pickup){
        if(!m.picked_up){
          // If within 1% of range, pick up
          float rng = std::max(1e-6f, mx-mn);
          if(std::fabs(target - current) <= 0.01f*rng){
            m.picked_up = 1;
          } else {
            // not picked up, ignore
            continue;
          }
        }
      } else {
        m.picked_up = 1;
      }

      float smoothed = smooth[i].process(target);
      set_param_runtime(port, smoothed);
    }
  }

  float get_port_value(uint32_t port){
    // return current “effective” param from params
    switch(port){
      case P_MASTER_LEVEL: return params.master;
      case P_ACCENT_AMOUNT: return params.accent;

      case P_BD_LEVEL: return params.bd_level;
      case P_BD_TUNE: return params.bd_tune;
      case P_BD_DECAY: return params.bd_decay;
      case P_BD_CLICK: return params.bd_click;

      case P_SD_LEVEL: return params.sd_level;
      case P_SD_TUNE: return params.sd_tune;
      case P_SD_DECAY: return params.sd_decay;
      case P_SD_TONE: return params.sd_tone;
      case P_SD_SNAPPY: return params.sd_snappy;

      case P_LT_LEVEL: return params.lt_level;
      case P_LT_TUNE: return params.lt_tune;
      case P_LT_DECAY: return params.lt_decay;

      case P_MT_LEVEL: return params.mt_level;
      case P_MT_TUNE: return params.mt_tune;
      case P_MT_DECAY: return params.mt_decay;

      case P_HT_LEVEL: return params.ht_level;
      case P_HT_TUNE: return params.ht_tune;
      case P_HT_DECAY: return params.ht_decay;

      case P_CH_LEVEL: return params.ch_level;
      case P_CH_TONE: return params.ch_tone;
      case P_CH_DECAY: return params.ch_decay;

      case P_OH_LEVEL: return params.oh_level;
      case P_OH_TONE: return params.oh_tone;
      case P_OH_DECAY: return params.oh_decay;

      case P_CP_LEVEL: return params.cp_level;
      case P_CP_TONE: return params.cp_tone;
      case P_CP_DECAY: return params.cp_decay;

      case P_CB_LEVEL: return params.cb_level;
      case P_CB_TUNE: return params.cb_tune;
      case P_CB_DECAY: return params.cb_decay;
      case P_CB_TONE: return params.cb_tone;

      case P_CY_LEVEL: return params.cy_level;
      case P_CY_TONE: return params.cy_tone;
      case P_CY_DECAY: return params.cy_decay;

      case P_RIM_LEVEL: return params.rim_level;
      case P_RIM_TONE: return params.rim_tone;

      case P_MA_LEVEL: return params.ma_level;
      case P_MA_TONE: return params.ma_tone;

      case P_CLV_LEVEL: return params.clv_level;

      default: return 0.0f;
    }
  }

  void set_param_runtime(uint32_t port, float v){
    // apply CC result into params (runtime)
    switch(port){
      case P_MASTER_LEVEL: params.master=v; break;
      case P_ACCENT_AMOUNT: params.accent=v; break;

      case P_BD_LEVEL: params.bd_level=v; break;
      case P_BD_TUNE: params.bd_tune=v; break;
      case P_BD_DECAY: params.bd_decay=v; break;
      case P_BD_CLICK: params.bd_click=v; break;

      case P_SD_LEVEL: params.sd_level=v; break;
      case P_SD_TUNE: params.sd_tune=v; break;
      case P_SD_DECAY: params.sd_decay=v; break;
      case P_SD_TONE: params.sd_tone=v; break;
      case P_SD_SNAPPY: params.sd_snappy=v; break;

      case P_LT_LEVEL: params.lt_level=v; break;
      case P_LT_TUNE: params.lt_tune=v; break;
      case P_LT_DECAY: params.lt_decay=v; break;

      case P_MT_LEVEL: params.mt_level=v; break;
      case P_MT_TUNE: params.mt_tune=v; break;
      case P_MT_DECAY: params.mt_decay=v; break;

      case P_HT_LEVEL: params.ht_level=v; break;
      case P_HT_TUNE: params.ht_tune=v; break;
      case P_HT_DECAY: params.ht_decay=v; break;

      case P_CH_LEVEL: params.ch_level=v; break;
      case P_CH_TONE: params.ch_tone=v; break;
      case P_CH_DECAY: params.ch_decay=v; break;

      case P_OH_LEVEL: params.oh_level=v; break;
      case P_OH_TONE: params.oh_tone=v; break;
      case P_OH_DECAY: params.oh_decay=v; break;

      case P_CP_LEVEL: params.cp_level=v; break;
      case P_CP_TONE: params.cp_tone=v; break;
      case P_CP_DECAY: params.cp_decay=v; break;

      case P_CB_LEVEL: params.cb_level=v; break;
      case P_CB_TUNE: params.cb_tune=v; break;
      case P_CB_DECAY: params.cb_decay=v; break;
      case P_CB_TONE: params.cb_tone=v; break;

      case P_CY_LEVEL: params.cy_level=v; break;
      case P_CY_TONE: params.cy_tone=v; break;
      case P_CY_DECAY: params.cy_decay=v; break;

      case P_RIM_LEVEL: params.rim_level=v; break;
      case P_RIM_TONE: params.rim_tone=v; break;

      case P_MA_LEVEL: params.ma_level=v; break;
      case P_MA_TONE: params.ma_tone=v; break;

      case P_CLV_LEVEL: params.clv_level=v; break;

      default: break;
    }
  }

  void arm_learn(int knobIndex){
    if(knobIndex < 0 || knobIndex >= (int)ccmap.size()) return;
    learn_knob_index = knobIndex;
    ccmap[knobIndex].picked_up = 0; // reset pickup latch on new mapping
  }

  void unlearn(int knobIndex){
    if(knobIndex < 0 || knobIndex >= (int)ccmap.size()) return;
    ccmap[knobIndex] = CCMapEntry{};
  }

  void handle_midi(const uint8_t* msg, uint32_t size){
    if(size < 1) return;
    uint8_t st = msg[0];
    uint8_t hi = st & 0xF0;

    if(hi == 0x90 && size >= 3){
      uint8_t note = msg[1] & 0x7F;
      uint8_t vel7 = msg[2] & 0x7F;
      if(vel7==0) return;

      float v = float(vel7)/127.0f;
      v *= clampf(params.accent, 0.0f, 2.0f);
      v = clampf(v, 0.0f, 1.0f);

      Inst inst = note_to_inst(note);
      if(inst != I_COUNT){
        dsp.trigger(inst, params);
        lastVel[inst] = v;
      }
    } else if(hi == 0xB0 && size >= 3){
      uint8_t cc = msg[1] & 0x7F;
      uint8_t val = msg[2] & 0x7F;

      if(learn_knob_index >= 0){
        // assign CC+channel from this event
        uint8_t ch = (st & 0x0F) + 1;
        ccmap[(size_t)learn_knob_index].cc = cc;
        ccmap[(size_t)learn_knob_index].chan = ch; // learn specific channel
        ccmap[(size_t)learn_knob_index].picked_up = 0;
        send_notify_learn(learn_knob_index, cc, ch);
        learn_knob_index = -1;
      } else {
        apply_cc(st, cc, val);
      }
    }
  }

  float lastVel[I_COUNT]{};

  void run(uint32_t n_samples){
    update_params_from_ports();

    auto clear=[&](float* b){ if(b) std::memset(b,0,sizeof(float)*n_samples); };
    clear(out_main_l); clear(out_main_r);
    clear(out_bd); clear(out_sd); clear(out_ch); clear(out_oh);
    clear(out_cp); clear(out_cb); clear(out_cy);
    clear(out_rim); clear(out_ma); clear(out_clv);
    clear(out_lt); clear(out_mt); clear(out_ht);

    uint32_t frame=0;
    if(midi_in){
      LV2_ATOM_SEQUENCE_FOREACH(midi_in, ev){
        uint32_t ev_frame = ev->time.frames;
        if(ev_frame > frame){
          render_block(frame, ev_frame, n_samples);
          frame = ev_frame;
        }
        if(ev->body.type == uris.midi_MidiEvent){
          const uint8_t* msg = (const uint8_t*)(ev+1);
          handle_midi(msg, ev->body.size);
        }
      }
    }
    if(frame < n_samples) render_block(frame, n_samples, n_samples);

    float m = clampf(params.master, 0.0f, 1.0f);
    for(uint32_t i=0;i<n_samples;i++){
      if(out_main_l) out_main_l[i] *= m;
      if(out_main_r) out_main_r[i] *= m;
    }
  }

  void render_block(uint32_t start, uint32_t end, uint32_t){
    for(uint32_t i=start;i<end;i++){
      float bd  = dsp.render_one(I_BD, params) * lastVel[I_BD];
      float sd  = dsp.render_one(I_SD, params) * lastVel[I_SD];
      float cp  = dsp.render_one(I_CP, params) * lastVel[I_CP];
      float rim = dsp.render_one(I_RIM,params) * lastVel[I_RIM];
      float ma  = dsp.render_one(I_MA, params) * lastVel[I_MA];
      float clv = dsp.render_one(I_CLV,params) * lastVel[I_CLV];
      float cb  = dsp.render_one(I_CB, params) * lastVel[I_CB];
      float ch  = dsp.render_one(I_CH, params) * lastVel[I_CH];
      float oh  = dsp.render_one(I_OH, params) * lastVel[I_OH];
      float cy  = dsp.render_one(I_CY, params) * lastVel[I_CY];
      float lt  = dsp.render_one(I_LT, params) * lastVel[I_LT];
      float mt  = dsp.render_one(I_MT, params) * lastVel[I_MT];
      float ht  = dsp.render_one(I_HT, params) * lastVel[I_HT];

      if(out_bd)  out_bd[i]=bd;
      if(out_sd)  out_sd[i]=sd;
      if(out_cp)  out_cp[i]=cp;
      if(out_rim) out_rim[i]=rim;
      if(out_ma)  out_ma[i]=ma;
      if(out_clv) out_clv[i]=clv;
      if(out_cb)  out_cb[i]=cb;
      if(out_ch)  out_ch[i]=ch;
      if(out_oh)  out_oh[i]=oh;
      if(out_cy)  out_cy[i]=cy;
      if(out_lt)  out_lt[i]=lt;
      if(out_mt)  out_mt[i]=mt;
      if(out_ht)  out_ht[i]=ht;

      float sum = bd+sd+cp+rim+ma+clv+cb+ch+oh+cy+lt+mt+ht;

      // dual mono main outs
      if(out_main_l) out_main_l[i]=sum;
      if(out_main_r) out_main_r[i]=sum;
    }
  }
};

// ===== LV2 State save/restore for CC mappings =====
static const char* KEY_CCMAP = "https://example.org/lv2/rs808#ccmap_v2";

static LV2_State_Status
state_save(LV2_Handle instance, LV2_State_Store_Function store, LV2_State_Handle handle, uint32_t, const LV2_Feature* const*) {
  auto* self = (RS808Plugin*)instance;
  LV2_URID key = self->map->map(self->map->handle, KEY_CCMAP);
  store(handle, key,
        self->ccmap.data(),
        (uint32_t)(self->ccmap.size() * sizeof(CCMapEntry)),
        self->uris.atom_Int,
        LV2_STATE_IS_POD | LV2_STATE_IS_PORTABLE);
  return LV2_STATE_SUCCESS;
}

static LV2_State_Status
state_restore(LV2_Handle instance, LV2_State_Retrieve_Function retrieve, LV2_State_Handle handle, uint32_t, const LV2_Feature* const*) {
  auto* self = (RS808Plugin*)instance;
  LV2_URID key = self->map->map(self->map->handle, KEY_CCMAP);
  size_t size=0; uint32_t type=0; uint32_t flags=0;
  const void* data = retrieve(handle, key, &size, &type, &flags);
  if(data){
    size_t need = self->ccmap.size() * sizeof(CCMapEntry);
    if(size == need){
      std::memcpy(self->ccmap.data(), data, size);
    }
  }
  return LV2_STATE_SUCCESS;
}

static const LV2_State_Interface state_iface = { state_save, state_restore };

static LV2_Handle
instantiate(const LV2_Descriptor*, double rate, const char*, const LV2_Feature* const* features){
  LV2_URID_Map* map=nullptr;
  for(int i=0; features && features[i]; ++i){
    if(!std::strcmp(features[i]->URI, LV2_URID__map)) map=(LV2_URID_Map*)features[i]->data;
  }
  if(!map) return nullptr;
  return (LV2_Handle)new RS808Plugin(rate, map);
}

static void connect_port(LV2_Handle instance, uint32_t port, void* data){
  auto* self=(RS808Plugin*)instance;
  switch(port){
    case P_MAIN_L: self->out_main_l=(float*)data; break;
    case P_MAIN_R: self->out_main_r=(float*)data; break;

    case P_BD: self->out_bd=(float*)data; break;
    case P_SD: self->out_sd=(float*)data; break;
    case P_CH: self->out_ch=(float*)data; break;
    case P_OH: self->out_oh=(float*)data; break;
    case P_CP: self->out_cp=(float*)data; break;
    case P_CB: self->out_cb=(float*)data; break;
    case P_CY: self->out_cy=(float*)data; break;
    case P_RIM: self->out_rim=(float*)data; break;
    case P_MA: self->out_ma=(float*)data; break;
    case P_CLV: self->out_clv=(float*)data; break;
    case P_LT: self->out_lt=(float*)data; break;
    case P_MT: self->out_mt=(float*)data; break;
    case P_HT: self->out_ht=(float*)data; break;

    case P_MIDI_IN: self->midi_in=(const LV2_Atom_Sequence*)data; break;
    case P_NOTIFY_OUT: self->notify_out=(LV2_Atom_Sequence*)data; break;

    default:
      if(port < P_COUNT) self->ctrl[port]=(const float*)data;
      break;
  }
}

static void run_cb(LV2_Handle instance, uint32_t n_samples){
  ((RS808Plugin*)instance)->run(n_samples);
}

static void cleanup(LV2_Handle instance){
  delete (RS808Plugin*)instance;
}

static const void* extension_data(const char* uri){
  if(!std::strcmp(uri, LV2_STATE__interface)) return &state_iface;
  return nullptr;
}

static const LV2_Descriptor descriptor = {
  RS808_URI,
  instantiate,
  connect_port,
  nullptr,
  run_cb,
  nullptr,
  cleanup,
  extension_data
};

extern "C" LV2_SYMBOL_EXPORT
const LV2_Descriptor* lv2_descriptor(uint32_t index){
  return index==0 ? &descriptor : nullptr;
}
