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
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <cstring>

static inline float clampf(float v, float lo, float hi) {
  return std::max(lo, std::min(hi, v));
}

// ===== Noise (WhiteNoise.ar) =====
struct WhiteNoise {
  uint32_t s = 22222u;
  float tick() {
    s = s * 1664525u + 1013904223u;
    uint32_t bits = (s >> 9) | 0x3f800000u;
    float f;
    std::memcpy(&f, &bits, sizeof(float));
    f = f - 1.5f;
    return f * 2.0f;
  }
};

// ===== Oscillators =====
struct SinOsc {
  float ph=0.0f, sr=48000.0f;
  void setSR(float s){ sr=s; }
  float tick(float hz, float phaseOffset=0.0f){
    float inc = (2.0f * float(M_PI) * hz) / sr;
    ph += inc;
    if (ph > 2.0f*float(M_PI)) ph -= 2.0f*float(M_PI);
    return std::sin(ph + phaseOffset);
  }
};

struct LFTri {
  float ph=0.0f, sr=48000.0f;
  void setSR(float s){ sr=s; }
  float tick(float hz, float iphase=0.0f){
    float inc = hz / sr;
    ph += inc;
    ph -= std::floor(ph);
    float p = ph + (iphase/(2.0f*float(M_PI)));
    p -= std::floor(p);
    float v = 4.0f * std::fabs(p - 0.5f) - 1.0f;
    return -v;
  }
};

struct LFPulse {
  float ph=0.0f, sr=48000.0f;
  void setSR(float s){ sr=s; }
  float tick(float hz, float width=0.5f, float iphase=0.0f){
    float inc = hz / sr;
    ph += inc;
    ph -= std::floor(ph);
    float p = ph + (iphase/(2.0f*float(M_PI)));
    p -= std::floor(p);
    return (p < width) ? 1.0f : 0.0f;
  }
};

// ===== Env curves =====
static inline float env_curve_interp(float a, float b, float t01, float curve) {
  t01 = clampf(t01, 0.0f, 1.0f);
  if (std::fabs(curve) < 1e-6f) return a + (b-a)*t01;
  float denom = 1.0f - std::exp(curve);
  if (std::fabs(denom) < 1e-12f) return a + (b-a)*t01;
  float num = 1.0f - std::exp(curve * t01);
  float shaped = num / denom;
  return a + (b-a) * shaped;
}

struct EnvSeg2 {
  float sr=48000.0f;
  float l0=0, l1=1, l2=0;
  float t1=0.001f, t2=0.001f;
  float curve=0.0f;
  float pos=0.0f;
  bool active=false;

  void setSR(float s){ sr=s; }

  void trigger(float _l0,float _l1,float _l2,float _t1,float _t2,float _curve){
    l0=_l0; l1=_l1; l2=_l2;
    t1=std::max(0.0f,_t1);
    t2=std::max(0.0f,_t2);
    curve=_curve;
    pos=0.0f;
    active=true;
  }

  float tick(){
    if(!active) return 0.0f;
    float out=0.0f;
    float total = t1 + t2;
    if (total <= 0.0f) { active=false; return l2; }

    if (pos < t1) {
      float t01 = (t1 <= 0.0f) ? 1.0f : (pos / t1);
      out = env_curve_interp(l0,l1,t01,curve);
    } else {
      float p2 = pos - t1;
      float t01 = (t2 <= 0.0f) ? 1.0f : (p2 / t2);
      out = env_curve_interp(l1,l2,t01,curve);
    }

    pos += 1.0f/sr;
    if (pos >= total) { active=false; out=l2; }
    return out;
  }

  bool isActive() const { return active; }
};

struct EnvPerc {
  EnvSeg2 e;
  void setSR(float s){ e.setSR(s); }
  void trigger(float atk, float rel, float level, float curve){
    e.trigger(0.0f, level, 0.0f, atk, rel, curve);
  }
  float tick(){ return e.tick(); }
  bool isActive() const { return e.isActive(); }
};

// ===== Filters =====
struct OnePoleLPF {
  float a=0.0f, y=0.0f, sr=48000.0f;
  void setSR(float s){ sr=s; }
  void setFreq(float hz){
    hz = std::max(5.0f, hz);
    float x = std::exp(-2.0f * float(M_PI) * hz / sr);
    a = x;
  }
  float process(float in){
    y = (1.0f - a) * in + a * y;
    return y;
  }
};

struct OnePoleHPF {
  float a=0.0f, y=0.0f, x1=0.0f, sr=48000.0f;
  void setSR(float s){ sr=s; }
  void setFreq(float hz){
    hz = std::max(5.0f, hz);
    float x = std::exp(-2.0f * float(M_PI) * hz / sr);
    a = x;
  }
  float process(float in){
    y = (1.0f - a) * (y + in - x1);
    x1 = in;
    return y;
  }
};

struct Biquad {
  float b0=1,b1=0,b2=0,a1=0,a2=0;
  float z1=0,z2=0;
  float process(float x){
    float y = b0*x + z1;
    z1 = b1*x - a1*y + z2;
    z2 = b2*x - a2*y;
    return y;
  }
  void setCoeffs(float _b0,float _b1,float _b2,float _a0,float _a1,float _a2){
    b0=_b0/_a0; b1=_b1/_a0; b2=_b2/_a0;
    a1=_a1/_a0; a2=_a2/_a0;
  }
};

static inline void biquad_bandpass(Biquad& bq, float sr, float freq, float rq){
  freq = clampf(freq, 20.0f, 20000.0f);
  rq = std::max(1e-4f, rq);
  float Q = std::max(0.05f, 1.0f / rq);
  float w0=2.0f*float(M_PI)*freq/sr;
  float cw=std::cos(w0), sw=std::sin(w0);
  float alpha=sw/(2.0f*Q);
  float b0 = alpha, b1 = 0.0f, b2 = -alpha;
  float a0 = 1.0f + alpha;
  float a1 = -2.0f * cw;
  float a2 = 1.0f - alpha;
  bq.setCoeffs(b0,b1,b2,a0,a1,a2);
}

static inline void biquad_highpass(Biquad& bq, float sr, float freq, float rq){
  freq = clampf(freq, 20.0f, 20000.0f);
  rq = std::max(1e-4f, rq);
  float Q = std::max(0.05f, 1.0f / rq);
  float w0=2.0f*float(M_PI)*freq/sr;
  float cw=std::cos(w0), sw=std::sin(w0);
  float alpha=sw/(2.0f*Q);
  float b0=(1+cw)/2, b1=-(1+cw), b2=(1+cw)/2;
  float a0=1+alpha, a1=-2*cw, a2=1-alpha;
  bq.setCoeffs(b0,b1,b2,a0,a1,a2);
}

static inline void biquad_peakeq(Biquad& bq, float sr, float freq, float rq, float dbGain){
  freq = clampf(freq, 20.0f, 20000.0f);
  rq = std::max(1e-4f, rq);
  float Q = std::max(0.05f, 1.0f / rq);
  float A = std::pow(10.0f, dbGain/40.0f);
  float w0=2.0f*float(M_PI)*freq/sr;
  float cw=std::cos(w0), sw=std::sin(w0);
  float alpha=sw/(2.0f*Q);

  float b0 = 1 + alpha*A;
  float b1 = -2*cw;
  float b2 = 1 - alpha*A;
  float a0 = 1 + alpha/A;
  float a1 = -2*cw;
  float a2 = 1 - alpha/A;
  bq.setCoeffs(b0,b1,b2,a0,a1,a2);
}

struct Hipass4 {
  Biquad b1, b2;
  float sr=48000.0f;
  void setSR(float s){ sr=s; }
  void set(float freq, float rq){
    biquad_highpass(b1, sr, freq, rq);
    biquad_highpass(b2, sr, freq, rq);
  }
  float process(float x){ return b2.process(b1.process(x)); }
};

static inline float limiter_like(float x, float level){
  float y = std::tanh(x);
  return clampf(y, -level, level);
}

// ===== Instruments =====
enum Inst : int {
  I_BD=0, I_SD, I_CP, I_RIM, I_CLV, I_MA, I_CB, I_CH, I_OH, I_CY,
  I_LT, I_MT, I_HT,
  I_COUNT
};

struct RS808_Params {
  // Global
  float master = 0.8f;
  float accent = 1.0f;

  // BD
  float bd_level=0.8f;
  float bd_tune=56.0f;   // Hz-ish base in original patch
  float bd_decay=30.0f;
  float bd_click=0.5f;   // new: scales punch/click

  // SD
  float sd_level=0.7f;
  float sd_tune=220.0f;  // new: shifts both osc freqs
  float sd_decay=1.0f;   // new: scales env releases
  float sd_tone=340.0f;  // osc2 freq target (original tone)
  float sd_snappy=0.3f;

  // LT/MT/HT
  float lt_level=0.6f, lt_tune=80.0f,  lt_decay=20.0f;
  float mt_level=0.6f, mt_tune=120.0f, mt_decay=15.0f;
  float ht_level=0.6f, ht_tune=165.0f, ht_decay=10.0f;

  // CH/OH
  float ch_level=0.55f, ch_tone=0.5f, ch_decay=0.42f;
  float oh_level=0.65f, oh_tone=0.5f, oh_decay=0.5f;

  // Clap
  float cp_level=0.6f, cp_tone=0.5f, cp_decay=0.5f;

  // Cowbell
  float cb_level=0.55f, cb_tune=1.0f, cb_decay=1.0f, cb_tone=0.5f;

  // Cymbal
  float cy_level=0.55f, cy_tone=0.5f, cy_decay=2.0f;

  // Rimshot/Maracas/Claves
  float rim_level=0.55f, rim_tone=0.5f;
  float ma_level=0.5f, ma_tone=0.5f;
  float clv_level=0.5f;
};

// BD: based on prior translation, add click scaler and tune as “tone”
struct VoiceBD {
  float sr=48000.0f;
  EnvSeg2 env, trienv, fenv, pfenv;
  SinOsc sin1, sin2;
  LFTri tri;
  OnePoleHPF punchHPF;
  WhiteNoise clickNoise;
  EnvPerc clickEnv;
  OnePoleHPF clickHPF;

  void setSR(float s){
    sr=s;
    env.setSR(s); trienv.setSR(s); fenv.setSR(s); pfenv.setSR(s);
    sin1.setSR(s); sin2.setSR(s); tri.setSR(s);
    punchHPF.setSR(s); punchHPF.setFreq(350.0f);
    clickEnv.setSR(s);
    clickHPF.setSR(s); clickHPF.setFreq(2500.0f);
  }

  void trigger(float decay, float tune){
    env.trigger(0.11f, 1.0f, 0.0f, 0.0f, decay, -225.0f);
    trienv.trigger(0.11f, 0.6f, 0.0f, 0.0f, decay, -230.0f);
    fenv.trigger(tune*7.0f, tune*1.35f, tune, 0.05f, 0.6f, -14.0f);
    pfenv.trigger(tune*7.0f, tune*1.35f, tune, 0.03f, 0.6f, -10.0f);
    clickEnv.trigger(0.001f, 0.015f, 1.0f, -60.0f);
  }

  float process(float level, float click){
    if(!env.isActive() && !trienv.isActive() && !clickEnv.isActive()) return 0.0f;

    float e  = env.tick();
    float te = trienv.tick();
    float fe = fenv.tick();
    float pfe= pfenv.tick();

    float sig = sin1.tick(fe, float(M_PI)/2.0f) * e;
    float sub = tri.tick(fe, float(M_PI)/2.0f) * te * 0.05f;

    float punch = sin2.tick(pfe, float(M_PI)/2.0f) * e * 2.0f;
    punch = punchHPF.process(punch);

    float clk = clickHPF.process(clickNoise.tick()) * clickEnv.tick() * 0.2f;

    float sum = (sig + sub + punch * (0.6f + 1.4f*click) + clk * click) * 2.5f;
    sum = limiter_like(sum, 0.5f) * level;
    return sum;
  }
};

// SD: add decay scaler + tune shifting
struct VoiceSD {
  float sr=48000.0f;
  EnvPerc noiseEnv;
  EnvPerc atkEnv;
  WhiteNoise wn;
  OnePoleHPF nHPF;
  OnePoleLPF nLPF;
  SinOsc osc1, osc2;
  OnePoleHPF outHPF;

  void setSR(float s){
    sr=s;
    noiseEnv.setSR(s); atkEnv.setSR(s);
    nHPF.setSR(s); nLPF.setSR(s); outHPF.setSR(s);
    nHPF.setFreq(1800.0f);
    nLPF.setFreq(8850.0f);
    outHPF.setFreq(340.0f);
    osc1.setSR(s); osc2.setSR(s);
  }

  void trigger(float decayMul){
    noiseEnv.trigger(0.001f, 4.2f*decayMul, 1.0f, -115.0f);
    atkEnv.trigger(0.001f, 0.8f*decayMul, 1.0f, -95.0f);
  }

  float process(float level, float tuneShift, float tone, float snappy){
    if(!noiseEnv.isActive() && !atkEnv.isActive()) return 0.0f;

    float ne = noiseEnv.tick();
    float ae = atkEnv.tick();

    float noise = wn.tick();
    noise = nHPF.process(noise);
    noise = nLPF.process(noise);
    noise = noise * ne * snappy;

    // tuneShift shifts both osc freqs musically (simple ratio around 1.0)
    float ratio = std::pow(2.0f, tuneShift/12.0f); // tuneShift in semis (-12..+12)
    float tone2 = (tone * 0.555f) * ratio; // keep relation similar to old tone2
    float tone1 = tone * ratio;

    float o1 = osc1.tick(tone2, float(M_PI)/2.0f) * 0.6f;
    float o2 = osc2.tick(tone1, float(M_PI)/2.0f) * 0.7f;
    float sum = (o1 + o2) * ae;

    float out = (noise + sum) * level * 2.5f;
    out = outHPF.process(out);
    return out;
  }
};

// Clap: add tone (bandpass freq) + decay scaling
struct VoiceClap {
  float sr=48000.0f;
  EnvSeg2 atkenv;
  float dadsr_pos=0.0f;
  bool dadsr_active=false;
  WhiteNoise wn1, wn2;
  OnePoleHPF hpf;
  Biquad bpf;

  EnvPerc revgen;
  WhiteNoise wnrev;
  OnePoleHPF revHPF;
  OnePoleLPF revLPF;

  void setSR(float s){
    sr=s;
    atkenv.setSR(s);
    hpf.setSR(s); hpf.setFreq(500.0f);
    revgen.setSR(s);
    revHPF.setSR(s); revHPF.setFreq(500.0f);
    revLPF.setSR(s); revLPF.setFreq(1000.0f);
  }

  void trigger(float decayMul){
    float d = clampf(decayMul, 0.05f, 4.0f);
    atkenv.trigger(0.5f, 1.0f, 0.0f, 0.0f, 0.3f*d, -160.0f);
    dadsr_pos=0.0f;
    dadsr_active=true;
    revgen.trigger(0.1f, 4.0f*d, 1.0f, -9.0f);
  }

  float dadsr_tick(float decayMul){
    if(!dadsr_active) return 0.0f;
    float t = dadsr_pos;
    dadsr_pos += 1.0f/sr;

    float delay=0.026f, attack=0.0f, decay=6.0f*decayMul;
    float sustain=0.0f, peak=1.0f, curve=-157.0f;

    if (t < delay) return 0.0f;
    t -= delay;

    if (t < attack) {
      float t01 = (attack<=0.0f)?1.0f:(t/attack);
      return env_curve_interp(0.0f, peak, t01, curve);
    }
    t -= attack;

    if (t < decay) {
      float t01 = (decay<=0.0f)?1.0f:(t/decay);
      return env_curve_interp(peak, sustain, t01, curve);
    }
    dadsr_active=false;
    return sustain;
  }

  float process(float level, float tone01, float decayMul){
    float atk = atkenv.tick();
    float denv = dadsr_tick(decayMul);

    float atk_sig = wn1.tick() * atk * 1.4f;
    float dec_sig = wn2.tick() * denv;
    float sum = atk_sig + dec_sig * level;
    sum = hpf.process(sum);

    // tone maps 0..1 -> bandpass 800..3000
    float f = 800.0f + 2200.0f*clampf(tone01,0,1);
    biquad_bandpass(bpf, sr, f, 0.5f);
    sum = bpf.process(sum);
    float raw = sum * 1.5f;

    float rg = revgen.tick();
    float rev = wnrev.tick() * rg * 0.02f;
    rev = revHPF.process(rev);
    rev = revLPF.process(rev);
    float rv = rev * level;

    return (raw + rv) * level;
  }
};

struct VoiceToneEnvSine {
  float sr=48000.0f;
  EnvSeg2 env;
  EnvSeg2 fenv;
  SinOsc osc;
  void setSR(float s){ sr=s; env.setSR(s); fenv.setSR(s); osc.setSR(s); }

  void trigger(float ampAtk, float ampRel, float curveEnv,
               float baseFreq){
    env.trigger(0.0f, 1.0f, 0.0f, 0.0f, ampRel, curveEnv);
    fenv.trigger(baseFreq*1.25f, baseFreq*1.125f, baseFreq, 0.1f, 0.5f, -4.0f);
               }

               float process(float level){
                 if(!env.isActive()) return 0.0f;
                 float e = env.tick();
                 float f = fenv.tick();
                 return osc.tick(f, float(M_PI)/2.0f) * e * level;
               }
};

struct VoiceRimshot {
  float sr=48000.0f;
  EnvSeg2 env;
  LFTri tri1;
  LFPulse pul;
  WhiteNoise wn;
  Biquad peak;
  OnePoleHPF hpf;
  OnePoleLPF lpf;

  void setSR(float s){
    sr=s; env.setSR(s);
    tri1.setSR(s); pul.setSR(s);
    hpf.setSR(s); lpf.setSR(s);
    hpf.setFreq(315.0f); lpf.setFreq(7300.0f);
  }

  void trigger(){ env.trigger(1.0f,1.0f,0.0f,0.00272f,0.07f,-42.0f); }

  float process(float level, float tone01){
    if(!env.isActive()) return 0.0f;
    float e = env.tick();

    float t1 = tri1.tick(1667.0f*1.1f, 1.0f) * e;
    float t2 = pul.tick(455.0f*1.1f, 0.8f) * e;
    float punch = wn.tick() * e * 0.46f;
    float sig = t1 + t2 + punch;

    // tone shifts peakEQ center 350..900
    float f = 350.0f + 550.0f*clampf(tone01,0,1);
    biquad_peakeq(peak, sr, f, 0.44f, 8.0f);
    sig = peak.process(sig);

    sig = hpf.process(sig);
    sig = lpf.process(sig);
    return sig * level;
  }
};

struct VoiceClaves {
  float sr=48000.0f;
  EnvSeg2 env;
  SinOsc osc;
  void setSR(float s){ sr=s; env.setSR(s); osc.setSR(s); }
  void trigger(){ env.trigger(1.0f,1.0f,0.0f,0.0f,0.1f,-20.0f); }
  float process(float level){
    if(!env.isActive()) return 0.0f;
    float e = env.tick();
    return osc.tick(2500.0f, float(M_PI)/2.0f) * e * level;
  }
};

struct VoiceMaracas {
  float sr=48000.0f;
  EnvSeg2 env;
  WhiteNoise wn;
  OnePoleHPF hpf;
  void setSR(float s){ sr=s; env.setSR(s); hpf.setSR(s); }
  void trigger(){ env.trigger(0.3f,1.0f,0.0f,0.027f,0.07f,-250.0f); }
  float process(float level, float tone01){
    if(!env.isActive()) return 0.0f;
    // tone shifts HPF 2500..9000
    float f = 2500.0f + 6500.0f*clampf(tone01,0,1);
    hpf.setFreq(f);
    float e = env.tick();
    float sig = wn.tick() * e * level;
    return hpf.process(sig);
  }
};

struct VoiceCowbell {
  float sr=48000.0f;
  EnvPerc atkenv;
  EnvPerc env;
  LFPulse p1, p2;
  OnePoleHPF hpf;
  OnePoleLPF lpf;

  void setSR(float s){
    sr=s;
    atkenv.setSR(s); env.setSR(s);
    p1.setSR(s); p2.setSR(s);
    hpf.setSR(s); lpf.setSR(s);
  }

  void trigger(float decayMul){
    float d=clampf(decayMul,0.1f,6.0f);
    atkenv.trigger(0.0f, 1.0f*d, 1.0f, -215.0f);
    env.trigger(0.01f, 9.5f*d, 1.0f, -90.0f);
  }

  float process(float level, float tuneMul, float tone01){
    if(!env.isActive() && !atkenv.isActive()) return 0.0f;

    float ae = atkenv.tick();
    float e  = env.tick();

    float f1 = 811.16f * tuneMul;
    float f2 = 538.75f * tuneMul;

    float pul1 = p1.tick(f1);
    float pul2 = p2.tick(f2);

    float atk = (pul1 + pul2) * ae * 6.0f;
    float datk = (pul1 + pul2) * e;

    float sig = (atk + datk) * level;

    // tone shifts HPF/LPF
    float hp = 150.0f + 500.0f*clampf(tone01,0,1);
    float lp = 2500.0f + 4500.0f*clampf(tone01,0,1);
    hpf.setFreq(hp);
    lpf.setFreq(lp);

    sig = hpf.process(sig);
    sig = lpf.process(sig);
    return sig;
  }
};

struct VoiceHatClosed {
  float sr=48000.0f;
  EnvPerc env;
  LFPulse o1,o2,o3,o4,o5,o6;
  Biquad bpf_hi;
  OnePoleHPF hpf_hi;
  Biquad bpf_lo;
  Biquad hpf_lo;
  Biquad peq;

  void setSR(float s){
    sr=s;
    env.setSR(s);
    o1.setSR(s);o2.setSR(s);o3.setSR(s);o4.setSR(s);o5.setSR(s);o6.setSR(s);
    hpf_hi.setSR(s);
  }

  void trigger(float decay){
    env.trigger(0.005f, decay, 1.0f, -30.0f);
  }

  float process(float level, float tone01){
    if(!env.isActive()) return 0.0f;
    float e = env.tick();

    float s = o1.tick(203.52f) + o2.tick(366.31f) + o3.tick(301.77f) +
    o4.tick(518.19f) + o5.tick(811.16f) + o6.tick(538.75f);

    // tone shifts hat band 6k..12k
    float band = 6000.0f + 6000.0f*clampf(tone01,0,1);

    float sighi = s;
    float siglow = s;

    biquad_bandpass(bpf_hi, sr, band, 1.0f);
    sighi = bpf_hi.process(sighi);
    hpf_hi.setFreq(band + 300.0f);
    sighi = hpf_hi.process(sighi);

    biquad_bandpass(bpf_lo, sr, band, 0.8f);
    siglow = bpf_lo.process(siglow);

    biquad_highpass(hpf_lo, sr, band + 300.0f, 0.3f);
    siglow = hpf_lo.process(siglow);

    float sum = siglow + sighi;
    biquad_peakeq(peq, sr, band + 800.0f, 0.8f, 0.7f);
    sum = peq.process(sum);

    return sum * e * level;
  }
};

struct VoiceHatOpen {
  float sr=48000.0f;
  EnvPerc env1;
  EnvSeg2 env2;
  LFPulse o1,o2,o3,o4,o5,o6;

  Biquad bpf;
  Biquad peq;
  Hipass4 hp4;
  OnePoleLPF lpf;

  void setSR(float s){
    sr=s;
    env1.setSR(s); env2.setSR(s);
    o1.setSR(s);o2.setSR(s);o3.setSR(s);o4.setSR(s);o5.setSR(s);o6.setSR(s);
    hp4.setSR(s);
    lpf.setSR(s);
  }

  void trigger(float decay){
    env1.trigger(0.1f, decay, 1.0f, -3.0f);
    env2.trigger(0.0f, 1.0f, 0.0f, 0.0f, decay*5.0f, -150.0f);
  }

  float process(float level, float tone01){
    if(!env1.isActive() && !env2.isActive()) return 0.0f;
    float e1 = env1.tick();
    float e2 = env2.tick();

    float sig = (o1.tick(203.52f)*0.6f + o2.tick(366.31f)*0.6f + o3.tick(301.77f)*0.6f +
    o4.tick(518.19f)*0.6f + o5.tick(811.16f)*0.6f + o6.tick(538.75f)*0.6f);

    float band = 6000.0f + 6000.0f*clampf(tone01,0,1);
    biquad_bandpass(bpf, sr, band, 1.0f);
    sig = bpf.process(sig);

    biquad_peakeq(peq, sr, band + 500.0f, 0.5f, 5.0f);
    sig = peq.process(sig);

    hp4.set(band + 400.0f, 0.7f);
    sig = hp4.process(sig);

    lpf.setFreq(4000.0f + 6000.0f*clampf(tone01,0,1));
    sig = lpf.process(sig);

    float sum = sig * (e1*0.6f + e2);
    return sum * level * 2.0f;
  }
};

struct VoiceCymbal {
  float sr=48000.0f;
  EnvPerc env1;
  EnvSeg2 env2;
  EnvSeg2 env3;
  LFPulse o1,o2,o3,o4,o5,o6;

  Biquad bp1, peq1;
  Hipass4 hp1;
  OnePoleLPF lpf;

  float det_phase=0.0f;
  float det=1.0f;

  void setSR(float s){
    sr=s;
    env1.setSR(s); env2.setSR(s); env3.setSR(s);
    o1.setSR(s);o2.setSR(s);o3.setSR(s);o4.setSR(s);o5.setSR(s);o6.setSR(s);
    hp1.setSR(s);
    lpf.setSR(s);
  }

  void trigger(float decay){
    env1.trigger(0.3f, decay, 1.0f, -3.0f);
    env2.trigger(0.0f,0.6f,0.0f,0.1f, decay*0.7f, -5.0f);
    env3.trigger(0.0f,1.0f,0.0f,0.0f, decay*5.0f, -150.0f);
    det_phase=0.0f; det=1.0f;
  }

  float detune_tick(float rate){
    rate = std::max(0.0001f, rate);
    det_phase += rate / sr;
    if(det_phase >= 1.0f){
      det_phase -= 1.0f;
      uint32_t s = 1234567u ^ (uint32_t)(det_phase*1e9);
      s = s * 1664525u + 1013904223u;
      float r = float((s >> 9) & 0x7FFFFF) / float(0x7FFFFF);
      det = 0.99f + r*0.02f;
    }
    return det;
  }

  float process(float level, float tone01, float decay){
    (void)decay;
    if(!env1.isActive() && !env2.isActive() && !env3.isActive()) return 0.0f;

    float e1=env1.tick();
    float e2=env2.tick();
    float e3=env3.tick();

    // tone controls brightness + detune rate (musical mapping)
    float bright = clampf(tone01,0,1);
    float detRate = 0.0005f + bright*0.05f;
    float d = detune_tick(detRate);

    float sig = (o1.tick(203.52f*d)*0.6f + o2.tick(366.31f*d)*0.6f + o3.tick(301.77f*d)*0.6f +
    o4.tick(518.19f*d)*0.6f + o5.tick(811.16f*d)*0.6f + o6.tick(538.75f*d)*0.6f);

    float band = 2500.0f + 8500.0f*bright;
    biquad_bandpass(bp1, sr, band, 1.0f);
    float x = bp1.process(sig);

    biquad_peakeq(peq1, sr, band+800.0f, 0.6f, 6.0f*bright);
    x = peq1.process(x);

    hp1.set(6500.0f + 5000.0f*bright, 0.8f);
    x = hp1.process(x);

    lpf.setFreq(4000.0f + 8000.0f*bright);
    x = lpf.process(x);

    return x * (e1 + 0.6f*e2 + 0.4f*e3) * level;
  }
};

struct RS808 {
  float sr=48000.0f;

  VoiceBD bd;
  VoiceSD sd;
  VoiceClap cp;
  VoiceRimshot rim;
  VoiceMaracas ma;
  VoiceClaves clv;
  VoiceCowbell cb;
  VoiceHatClosed ch;
  VoiceHatOpen oh;
  VoiceCymbal cy;

  VoiceToneEnvSine lt, mt, ht;

  void setSR(float s){
    sr=s;
    bd.setSR(s); sd.setSR(s); cp.setSR(s);
    rim.setSR(s); ma.setSR(s); clv.setSR(s); cb.setSR(s);
    ch.setSR(s); oh.setSR(s); cy.setSR(s);
    lt.setSR(s); mt.setSR(s); ht.setSR(s);
  }

  void trigger(Inst i, const RS808_Params& p){
    switch(i){
      case I_BD: bd.trigger(p.bd_decay, p.bd_tune); break;
      case I_SD: sd.trigger(p.sd_decay); break;
      case I_CP: cp.trigger(p.cp_decay); break;
      case I_RIM: rim.trigger(); break;
      case I_MA: ma.trigger(); break;
      case I_CLV: clv.trigger(); break;
      case I_CB: cb.trigger(p.cb_decay); break;
      case I_CH: ch.trigger(p.ch_decay); break;
      case I_OH: oh.trigger(p.oh_decay); break;
      case I_CY: cy.trigger(p.cy_decay); break;

      case I_LT: lt.trigger(0, p.lt_decay, -250.0f, p.lt_tune); break;
      case I_MT: mt.trigger(0, p.mt_decay, -250.0f, p.mt_tune); break;
      case I_HT: ht.trigger(0, p.ht_decay, -250.0f, p.ht_tune); break;
      default: break;
    }
  }

  float render_one(Inst i, const RS808_Params& p){
    switch(i){
      case I_BD: return bd.process(p.bd_level, p.bd_click);
      case I_SD: return sd.process(p.sd_level, p.sd_tune, p.sd_tone, p.sd_snappy);
      case I_CP: return cp.process(p.cp_level, p.cp_tone, p.cp_decay);
      case I_RIM: return rim.process(p.rim_level, p.rim_tone);
      case I_MA: return ma.process(p.ma_level, p.ma_tone);
      case I_CLV: return clv.process(p.clv_level);
      case I_CB: return cb.process(p.cb_level, p.cb_tune, p.cb_tone);
      case I_CH: return ch.process(p.ch_level, p.ch_tone);
      case I_OH: return oh.process(p.oh_level, p.oh_tone);
      case I_CY: return cy.process(p.cy_level, p.cy_tone, p.cy_decay);
      case I_LT: return lt.process(p.lt_level * 2.5f);
      case I_MT: return mt.process(p.mt_level * 2.0f);
      case I_HT: return ht.process(p.ht_level * 2.0f);
      default: return 0.0f;
    }
  }
};
