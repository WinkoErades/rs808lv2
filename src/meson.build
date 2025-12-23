lv2_dep = dependency('lv2', required: true)
gtk_dep = dependency('gtk+-3.0', required: true)
cairo_dep = dependency('cairo', required: true)

cpp_args = ['-fvisibility=hidden']

bundle_dir = join_paths(get_option('prefix'), 'lib', 'lv2', 'rs808.lv2')

dsp_sources = [
  'rs808_lv2.cpp',
  'rs808_sc_dsp.cpp',
  'rs808_urid.cpp',
  'rs808_state.cpp',
]

ui_sources = [
  'ui/rs808_ui.cpp',
  'ui/knob.cpp',
]

# DSP binary
shared_library('rs808',
  dsp_sources,
  dependencies: [lv2_dep],
  cpp_args: cpp_args,
  install: true,
  install_dir: bundle_dir,
  name_prefix: '',
  name_suffix: 'so',
)

# UI binary
shared_library('rs808_ui',
  ui_sources,
  dependencies: [lv2_dep, gtk_dep, cairo_dep],
  cpp_args: cpp_args,
  install: true,
  install_dir: bundle_dir,
  name_prefix: '',
  name_suffix: 'so',
)

# TTL files (plain data install)
install_data(
  ['../ttl/manifest.ttl', '../ttl/rs808.ttl', '../ttl/ui.ttl'],
  install_dir: bundle_dir
)
