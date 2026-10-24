# The tests CI leaves out (label `extended`; CI runs `ctest -LE extended`).
#
# CI is there to catch regressions, not to validate accuracy, so it runs a subset: the goldens of the
# default path and one pin per model, one windowed-versus-marched pair per gain kind, one degenerate
# case per model, the refusals and warnings, and the mutation-guarded tests.  Left out are the checks
# against reference solutions, the second and third variants of one check (point, plate, table,
# scattered receivers), convergence and resolution checks, and runs whose cost buys no bug a kept
# test would miss.  `ctest` with no label filter still runs everything; run it before a merge.
#
# The Catch2 cases tagged [extended] are the CTest unit_test_extended, listed here too.

set(EXTENDED_TESTS
  unit_test_extended
  # Urban: the whole-pass sampling check, the wave, plain and GeoJSON building tables (the diffuse
  # table, a superset at s = 0.2, stays), the canyon table's second exactness pairs (the bump pair
  # stays), and option runs whose check would pass with the option ignored.
  buildings_every_state_run buildings_sampled_run buildings_states_sampled urban_table_wave_build
  urban_table_wave urban_table_wave_marched urban_table_wave_windowed_is_exact
  urban_table_wave_golden urban_buildings_build urban_buildings urban_buildings_marched
  urban_buildings_plate urban_buildings_plate_marched urban_buildings_windowed_is_exact
  urban_buildings_plate_windowed_is_exact urban_buildings_golden urban_buildings_plate_golden
  urban_buildings_geojson_build urban_buildings_geojson urban_buildings_geojson_is_the_csv
  urban_table_plate_reference urban_table_plate_open urban_table_plate_open_is_open_ground
  urban_table_canyon_marched urban_table_windowed_is_exact urban_table_plate_marched
  urban_table_plate_windowed_is_exact buildings_scattering_out_of_reach_run
  buildings_scattering_out_of_reach_is_specular buildings_scattering_spectrum_run
  buildings_no_profile_pe_run buildings_no_profile_pe_reaches_no_row
  buildings_no_ground_effect_run buildings_no_ground_effect_is_the_default
  buildings_path_radius_run buildings_path_radius_is_close buildings_turbulence_run
  buildings_turbulence_is_small_on_the_demo buildings_pitched_runs buildings_relief_flat_run
  buildings_relief_flat_is_plain buildings_far_classes_run buildings_far_classes_reference
  buildings_far_classes_are_open point_field_sources_exact_runs point_field_shared_faces_sampled
  fdtd_field_order2_runs
  # Refraction, terrain, ground, tables, atmosphere: second variants of a check (plate, scattered,
  # table, classes) where the point one stays, the stride and resolution checks, marched pairs of a
  # gain whose pair stays elsewhere, and the per-track identities beyond one.
  air_clock_alone air_clock_is_each_tracks_own air_clock_per_track air_first_wind_is_its_own
  air_second_wind_is_its_own air_wind_alone_a air_wind_alone_b air_wind_per_track
  atmosphere_homogeneous_scattered atmosphere_uniform_scattered
  atmosphere_uniform_scattered_is_homogeneous build_atmosphere_manifest_records_provenance
  build_atmosphere_stratified build_atmosphere_stratified_marched
  build_atmosphere_stratified_windowed_is_exact ground_class_scalar_plate
  ground_class_single_is_scalar_plate ground_class_single_is_scalar_terrain
  ground_class_single_plate ground_classes_terrain_impedance_is_the_observers
  ground_classes_terrain_windowed_is_exact propagation_table_carries_harmonoise
  propagation_table_carries_turbulence propagation_table_harmonoise_reference
  propagation_table_harmonoise_run propagation_table_is_the_ground_model
  propagation_table_midpoint_validation propagation_table_plate
  propagation_table_plate_is_the_plate_model propagation_table_plate_marched
  propagation_table_plate_reference propagation_table_plate_windowed_is_exact
  propagation_table_reference propagation_table_run propagation_table_turbulence_matters
  propagation_table_turbulent_reference propagation_table_turbulent_run
  refraction_plate_windowed_is_exact refraction_point_windowed_is_exact
  refraction_uniform_plate_is_straight refraction_uniform_plate_scattered refraction_windy_plate
  refraction_windy_plate_marched refraction_windy_point refraction_windy_point_marched
  scattering_is_not_off straight_plate_scattered stratified_calm_floor_plate
  stratified_calm_floor_plate_is_floorless stratified_duct_plate stratified_duct_plate_golden
  stratified_duct_plate_marched stratified_duct_plate_windowed_is_exact
  stratified_duct_point_exact stratified_duct_point_exact_golden stratified_duct_point_marched
  stratified_duct_point_windowed_is_exact stratified_resolved_is_not_legacy
  stratified_resolved_plate stratified_resolved_plate_marched
  stratified_resolved_plate_windowed_is_exact stratified_resolved_point
  stratified_resolved_point_marched stratified_resolved_point_windowed_is_exact
  stratified_resolved_straight_point stratified_resolved_straight_point_marched
  stratified_resolved_straight_point_windowed_is_exact stratified_stride_duct_is_close
  stratified_stride_is_not_exact stratified_stride_windy_is_close stratified_windy_plate
  stratified_windy_point_exact stratified_windy_point_exact_golden terrain_class_scalar
  terrain_class_single terrain_classes terrain_classes_marched terrain_flat_plate
  terrain_flat_plate_is_flat terrain_flat_plate_reference terrain_flat_plate_screened
  terrain_flat_plate_screened_is_flat terrain_flat_scattered terrain_flat_scattered_is_flat
  terrain_flat_scattered_reference terrain_flat_scattered_screened
  terrain_flat_scattered_screened_is_flat terrain_flat_table_plate
  terrain_flat_table_plate_is_flat terrain_flat_table_plate_reference terrain_flat_table_point
  terrain_flat_table_point_is_flat terrain_flat_table_point_reference terrain_flat_table_scattered
  terrain_flat_table_scattered_is_flat terrain_flat_table_scattered_reference
  terrain_hills_plate_marched terrain_hills_plate_windowed_is_exact terrain_hills_point_marched
  terrain_hills_point_screened terrain_hills_point_windowed_is_exact terrain_hills_scattered
  terrain_hills_scattered_placed terrain_hills_table_bump_plate
  terrain_hills_table_bump_plate_marched terrain_hills_table_bump_plate_windowed_is_exact
  terrain_hills_table_plate terrain_hills_table_plate_is_the_plate_model
  terrain_hills_table_plate_marched terrain_hills_table_plate_windowed_is_exact
  terrain_hills_table_point terrain_hills_table_point_marched
  terrain_hills_table_point_windowed_is_exact terrain_hills_unblocked_is_unscreened
  terrain_reaches_the_ground_model terrain_ridge_screened_classes
  terrain_ridge_screened_classes_marched terrain_ridge_screened_classes_windowed_is_exact
  terrain_ridge_screened_noground terrain_ridge_screened_noground_marched
  terrain_ridge_screened_noground_windowed_is_exact terrain_ridge_screened_plate_marched
  terrain_ridge_screened_plate_windowed_is_exact terrain_ridge_screened_point_marched
  terrain_ridge_screened_point_windowed_is_exact terrain_ridge_screened_table
  terrain_ridge_screened_table_marched terrain_ridge_screened_table_windowed_is_exact
  terrain_tilted_point terrain_tilted_point_is_unscreened terrain_tilted_point_screened wind_built
  wind_is_discretization_independent wind_is_its_written_profile wind_refined wind_replayed
  # The rest: blends and nose headings (their unit tests stay), fleets, local time and periods
  # beyond the metric goldens, histories and contours, provenance, and duplicate refusals.
  blend_barycentric blend_nearest blend_aam blend_barycentric_is_golden
  blend_nearest_moves_footprint blend_aam_moves_footprint blend_aam_is_not_nearest
  blend_aam_spectral blend_aam_spectral_marched blend_nearest_spectral
  blend_nearest_spectral_marched load_factor_two_sided blend_aam_windowed_is_exact
  blend_nearest_windowed_is_exact load_factor_two_sided_moves_footprint nose_default_on_track
  nose_heading_on_track nose_heading_crabbed nose_database_heading nose_database_heading_explicit
  nose_database_velocity nose_heading_on_track_is_velocity nose_heading_crabbed_differs
  nose_database_velocity_differs nose_database_heading_is_default nose_hover_crabbed
  nose_hover_straight nose_hover_takes_the_frame nose_heading_needs_a_frame
  nose_heading_needs_a_heading geodetic_heading_default geodetic_heading_on_track
  geodetic_heading_crabbed geodetic_heading_on_track_is_velocity geodetic_heading_crabbed_differs
  operations_fleet_run operations_fleet_repeated_t operations_fleet_mix
  operations_fleet_logs_its_paths operations_fleet_loads_each_once operations_count_one_run
  operations_count_one_is_bit_identical operations_offset_zero_run
  operations_offset_zero_is_bit_identical operations_offset_zero_geodetic_run
  operations_offset_zero_geodetic_t operations_offset_zero_geodetic_is_bit_identical
  operations_local_noon operations_local_night operations_local_hold
  operations_local_time_replaces_clock operations_local_time_night_adds_10
  operations_local_time_dates_a_static_hold static_hold_dnl_with_start_utc
  operations_period_clockless_runs operations_period_clockless_night_run
  operations_period_clockless_is_the_night operations_clockless_needs_local_time
  operations_manifest_run operations_manifest_records_the_file operations_offset_needs_travel
  operations_offset_on_a_static_hold period_noon period_noon_day_window period_night
  period_evening period_lnight_is_laeq_23_07 period_levening_over_4_hours
  period_lday_over_12_hours period_lday_is_laeq_07_19 period_lday_empty_at_night
  period_laeq_whole_day_is_lden_by_day period_lnight_over_8_hours epnl_multi_trajectory_run
  epnl_multi_trajectory_names band_history_straight band_history_refracted band_history_windowed
  band_history_leaves_footprint band_history_refracted_is_straight band_history_needs_spectral
  contour_run contour_dnl contour_writes_geojson contour_dnl_levels procedure_runs
  mask_suspect_cells_clean mask_suspect_cells_clean_is_golden full_history_runs
  nax_file_is_written lamax_file_is_written dnl_with_start_utc laeq_needs_a_window
  lden_needs_a_clock dnl_needs_a_clock band_frames_need_spectra observer_mixed_widths_rejected
  lamax_weighting_runs_on_a_finer_step lamax_weighting_refuses_a_coarse_history
  vehicle_weight_mismatch_warns matched_vehicle_is_silent shadow_phase_needs_resolved
  shadow_model_needs_stratified epnl_needs_its_bands procedure_check_refuses_wrong_vehicle
  csel_refuses_db_cutoff csel_refuses_ambient nax_needs_a_threshold
  screening_min_height_needs_screening provenance_run provenance_windowed
  provenance_spheres_silent_run provenance_spheres_match provenance_spheres_silent
  provenance_shares_held_states provenance_leaves_sel_alone static_hover_runs
)
# Tests that exist only with a GPU backend, which read extended runs.
set(EXTENDED_PLATFORM_TESTS
  urban_table_wave_gpu_build urban_table_wave_gpu urban_table_wave_gpu_matches_cpu
  buildings_gpu_run buildings_gpu_every_state_run buildings_gpu_matches_cpu buildings_gpu_every_state_matches_cpu
)
foreach(test IN LISTS EXTENDED_PLATFORM_TESTS)
  if(TEST ${test})
    list(APPEND EXTENDED_TESTS ${test})
  endif()
endforeach()

foreach(test IN LISTS EXTENDED_TESTS)
  if(NOT TEST ${test})
    message(FATAL_ERROR "tests/extended.cmake names ${test}, which is not a test")
  endif()
endforeach()
set_tests_properties(${EXTENDED_TESTS} PROPERTIES LABELS extended)

# CTest runs every setup of a fixture that a selected test requires, whatever its label, so an
# extended setup in a fixture a kept test requires would run anyway.  Such a fixture is split here:
# its extended setups set up <fixture>_extended instead, which its extended consumers require as well.
# A kept test must still find a kept setup for each of its fixtures.
get_property(all_tests DIRECTORY PROPERTY TESTS)
set(fixtures)
foreach(test IN LISTS all_tests)
  get_test_property(${test} FIXTURES_SETUP setups)
  if(setups)
    list(APPEND fixtures ${setups})
  endif()
endforeach()
list(REMOVE_DUPLICATES fixtures)
foreach(fixture IN LISTS fixtures)
  set(kept_setups)
  set(extended_setups)
  set(kept_users)
  set(extended_users)
  foreach(test IN LISTS all_tests)
    get_test_property(${test} FIXTURES_SETUP setups)
    get_test_property(${test} FIXTURES_REQUIRED required)
    if(test IN_LIST EXTENDED_TESTS)
      set(group extended)
    else()
      set(group kept)
    endif()
    if(setups AND fixture IN_LIST setups)
      list(APPEND ${group}_setups ${test})
    endif()
    if(required AND fixture IN_LIST required)
      list(APPEND ${group}_users ${test})
    endif()
  endforeach()
  if(kept_users AND NOT kept_setups)
    list(GET kept_users 0 user)
    message(FATAL_ERROR "${user} runs on CI but its fixture ${fixture} is set up only by extended "
            "tests (tests/extended.cmake): label it extended too, or keep a setup")
  endif()
  if(kept_users AND extended_setups)
    foreach(test IN LISTS extended_setups)
      get_test_property(${test} FIXTURES_SETUP setups)
      list(TRANSFORM setups REPLACE "^${fixture}$" "${fixture}_extended")
      set_tests_properties(${test} PROPERTIES FIXTURES_SETUP "${setups}")
    endforeach()
    foreach(test IN LISTS extended_users)
      get_test_property(${test} FIXTURES_REQUIRED required)
      list(APPEND required ${fixture}_extended)
      set_tests_properties(${test} PROPERTIES FIXTURES_REQUIRED "${required}")
    endforeach()
  endif()
endforeach()
