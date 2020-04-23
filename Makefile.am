AM_CPPFLAGS = -I${top_srcdir}/drifters -I${top_srcdir}/include -I${top_srcdir}/mosaic -I${top_srcdir}/mpp/include
AM_FCFLAGS = -I${top_srcdir}/include -I${top_srcdir}/mpp/include -I${top_srcdir}/fms -I${top_srcdir}/fms2_io/include

lib_LTLIBRARIES = libfmsc.la libFMS.la
libfmsc_la_SOURCES = ${top_srcdir}/affinity/affinity.c \
                     ${top_srcdir}/mosaic/create_xgrid.c \
                     ${top_srcdir}/mosaic/gradient_c2l.c \
                     ${top_srcdir}/mosaic/interp.c \
                     ${top_srcdir}/mosaic/mosaic_util.c \
                     ${top_srcdir}/mosaic/read_mosaic.c \
                     ${top_srcdir}/mpp/mpp_memuse.c \
                     ${top_srcdir}/mpp/nsclock.c \
                     ${top_srcdir}/mpp/threadloc.c

libFMS_la_SOURCES = ${top_srcdir}/affinity/fms_affinity.F90 \
    ${top_srcdir}/amip_interp/amip_interp.F90 \
    ${top_srcdir}/astronomy/astronomy.F90 \
    ${top_srcdir}/axis_utils/axis_utils.F90 \
    ${top_srcdir}/axis_utils/axis_utils2.F90 \
    ${top_srcdir}/block_control/block_control.F90 \
    ${top_srcdir}/column_diagnostics/column_diagnostics.F90 \
    ${top_srcdir}/constants/constants.F90 \
    ${top_srcdir}/coupler/atmos_ocean_fluxes.F90 \
    ${top_srcdir}/coupler/coupler_types.F90 \
    ${top_srcdir}/coupler/ensemble_manager.F90 \
    ${top_srcdir}/data_override/data_override.F90 \
    ${top_srcdir}/diag_integral/diag_integral.F90 \
    ${top_srcdir}/diag_manager/diag_axis.F90 \
    ${top_srcdir}/diag_manager/diag_data.F90 \
    ${top_srcdir}/diag_manager/diag_grid.F90 \
    ${top_srcdir}/diag_manager/diag_manager.F90 \
    ${top_srcdir}/diag_manager/diag_output.F90 \
    ${top_srcdir}/diag_manager/diag_table.F90 \
    ${top_srcdir}/diag_manager/diag_util.F90 \
    ${top_srcdir}/drifters/cloud_interpolator.F90 \
    ${top_srcdir}/drifters/drifters.F90 \
    ${top_srcdir}/drifters/drifters_comm.F90 \
    ${top_srcdir}/drifters/drifters_core.F90 \
    ${top_srcdir}/drifters/drifters_input.F90 \
    ${top_srcdir}/drifters/drifters_io.F90 \
    ${top_srcdir}/drifters/quicksort.F90 \
    ${top_srcdir}/exchange/stock_constants.F90 \
    ${top_srcdir}/exchange/xgrid.F90 \
    ${top_srcdir}/field_manager/field_manager.F90 \
    ${top_srcdir}/field_manager/fm_util.F90 \
    ${top_srcdir}/fms/fms.F90 \
    ${top_srcdir}/fms/fms_io.F90 \
    ${top_srcdir}/fms2_io/blackboxio.F90 \
    ${top_srcdir}/fms2_io/fms2_io.F90 \
    ${top_srcdir}/fms2_io/fms_io_utils.F90 \
    ${top_srcdir}/fms2_io/fms_netcdf_domain_io.F90 \
    ${top_srcdir}/fms2_io/fms_netcdf_unstructured_domain_io.F90 \
    ${top_srcdir}/fms2_io/netcdf_io.F90 \
    ${top_srcdir}/horiz_interp/horiz_interp.F90 \
    ${top_srcdir}/horiz_interp/horiz_interp_bicubic.F90 \
    ${top_srcdir}/horiz_interp/horiz_interp_bilinear.F90 \
    ${top_srcdir}/horiz_interp/horiz_interp_conserve.F90 \
    ${top_srcdir}/horiz_interp/horiz_interp_spherical.F90 \
    ${top_srcdir}/horiz_interp/horiz_interp_type.F90 \
    ${top_srcdir}/interpolator/interpolator.F90 \
    ${top_srcdir}/memutils/memutils.F90 \
    ${top_srcdir}/monin_obukhov/monin_obukhov.F90 \
    ${top_srcdir}/monin_obukhov/monin_obukhov_kernel.F90 \
    ${top_srcdir}/mosaic/gradient.F90 \
    ${top_srcdir}/mosaic/grid.F90 \
    ${top_srcdir}/mosaic/mosaic.F90 \
    ${top_srcdir}/mosaic/mosaic2.F90 \
    ${top_srcdir}/mpp/mpp.F90 \
    ${top_srcdir}/mpp/mpp_data.F90 \
    ${top_srcdir}/mpp/mpp_domains.F90 \
    ${top_srcdir}/mpp/mpp_efp.F90 \
    ${top_srcdir}/mpp/mpp_io.F90 \
    ${top_srcdir}/mpp/mpp_memutils.F90 \
    ${top_srcdir}/mpp/mpp_parameter.F90 \
    ${top_srcdir}/mpp/mpp_pset.F90 \
    ${top_srcdir}/mpp/mpp_utilities.F90 \
    ${top_srcdir}/platform/platform.F90 \
    ${top_srcdir}/random_numbers/MersenneTwister.F90 \
    ${top_srcdir}/random_numbers/random_numbers.F90 \
    ${top_srcdir}/sat_vapor_pres/sat_vapor_pres.F90 \
    ${top_srcdir}/sat_vapor_pres/sat_vapor_pres_k.F90 \
    ${top_srcdir}/time_interp/time_interp.F90 \
    ${top_srcdir}/time_interp/time_interp_external.F90 \
    ${top_srcdir}/time_interp/time_interp_external2.F90 \
    ${top_srcdir}/time_manager/get_cal_time.F90 \
    ${top_srcdir}/time_manager/time_manager.F90 \
    ${top_srcdir}/topography/gaussian_topog.F90 \
    ${top_srcdir}/topography/topography.F90 \
    ${top_srcdir}/tracer_manager/tracer_manager.F90 \
    ${top_srcdir}/tridiagonal/tridiagonal.F90
libFMS_la_LIBADD = libfmsc.la



${top_builddir}/affinity/fms_affinity.lo: ${top_srcdir}/affinity/fms_affinity.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/fms/fms.lo
${top_builddir}/amip_interp/amip_interp.lo: ${top_srcdir}/amip_interp/amip_interp.F90 ${top_srcdir}/time_interp/time_interp.lo ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/time_manager/get_cal_time.lo ${top_srcdir}/horiz_interp/horiz_interp.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/constants/constants.lo ${top_srcdir}/platform/platform.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/fms2_io/fms2_io.lo
${top_builddir}/astronomy/astronomy.lo: ${top_srcdir}/astronomy/astronomy.F90 ${top_srcdir}/fms/fms.lo ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/constants/constants.lo ${top_srcdir}/mpp/mpp.lo
${top_builddir}/axis_utils/axis_utils.lo: ${top_srcdir}/axis_utils/axis_utils.F90 ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/fms/fms.lo
${top_builddir}/axis_utils/axis_utils2.lo: ${top_srcdir}/axis_utils/axis_utils2.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/fms2_io/fms2_io.lo
${top_builddir}/block_control/block_control.lo: ${top_srcdir}/block_control/block_control.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo
${top_builddir}/column_diagnostics/column_diagnostics.lo: ${top_srcdir}/column_diagnostics/column_diagnostics.F90 ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/constants/constants.lo ${top_srcdir}/mpp/mpp.lo
${top_builddir}/constants/constants.lo: ${top_srcdir}/constants/constants.F90 ${top_srcdir}/platform/platform.lo
${top_builddir}/coupler/atmos_ocean_fluxes.lo: ${top_srcdir}/coupler/atmos_ocean_fluxes.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/coupler/coupler_types.lo ${top_srcdir}/coupler/coupler_types.lo ${top_srcdir}/coupler/coupler_types.lo ${top_srcdir}/coupler/coupler_types.lo ${top_srcdir}/coupler/coupler_types.lo ${top_srcdir}/coupler/coupler_types.lo ${top_srcdir}/field_manager/field_manager.lo ${top_srcdir}/field_manager/field_manager.lo ${top_srcdir}/field_manager/field_manager.lo ${top_srcdir}/field_manager/field_manager.lo ${top_srcdir}/field_manager/fm_util.lo ${top_srcdir}/field_manager/fm_util.lo ${top_srcdir}/field_manager/fm_util.lo ${top_srcdir}/field_manager/fm_util.lo ${top_srcdir}/field_manager/fm_util.lo ${top_srcdir}/field_manager/fm_util.lo ${top_srcdir}/field_manager/fm_util.lo ${top_srcdir}/field_manager/fm_util.lo ${top_srcdir}/field_manager/fm_util.lo
${top_builddir}/coupler/coupler_types.lo: ${top_srcdir}/coupler/coupler_types.F90 ${top_srcdir}/fms/fms.lo ${top_srcdir}/fms/fms_io.lo ${top_srcdir}/fms/fms_io.lo ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/diag_manager/diag_manager.lo ${top_srcdir}/data_override/data_override.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp.lo
${top_builddir}/coupler/ensemble_manager.lo: ${top_srcdir}/coupler/ensemble_manager.F90 ${top_srcdir}/fms/fms.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/fms/fms_io.lo
${top_builddir}/data_override/data_override.lo: ${top_srcdir}/data_override/data_override.F90 ${top_srcdir}/platform/platform.lo ${top_srcdir}/constants/constants.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/horiz_interp/horiz_interp.lo ${top_srcdir}/time_interp/time_interp_external2.lo ${top_srcdir}/fms/fms_io.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/axis_utils/axis_utils2.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/fms2_io/fms2_io.lo ${top_srcdir}/mosaic/mosaic2.lo
${top_builddir}/diag_integral/diag_integral.lo: ${top_srcdir}/diag_integral/diag_integral.F90 ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/fms2_io/fms2_io.lo ${top_srcdir}/constants/constants.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/coupler/ensemble_manager.lo
${top_builddir}/diag_manager/diag_axis.lo: ${top_srcdir}/diag_manager/diag_axis.F90 ${top_srcdir}/platform/platform.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/diag_manager/diag_data.lo
${top_builddir}/diag_manager/diag_data.lo: ${top_srcdir}/diag_manager/diag_data.F90 ${top_srcdir}/platform/platform.lo ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/fms2_io/fms2_io.lo
${top_builddir}/diag_manager/diag_grid.lo: ${top_srcdir}/diag_manager/diag_grid.F90 ${top_srcdir}/platform/platform.lo ${top_srcdir}/constants/constants.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo
${top_builddir}/diag_manager/diag_manager.lo: ${top_srcdir}/diag_manager/diag_manager.F90 ${top_srcdir}/platform/platform.lo ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/fms/fms_io.lo ${top_srcdir}/diag_manager/diag_axis.lo ${top_srcdir}/diag_manager/diag_util.lo ${top_srcdir}/diag_manager/diag_data.lo ${top_srcdir}/diag_manager/diag_table.lo ${top_srcdir}/diag_manager/diag_output.lo ${top_srcdir}/diag_manager/diag_grid.lo ${top_srcdir}/constants/constants.lo ${top_srcdir}/diag_manager/diag_axis.lo ${top_srcdir}/diag_manager/diag_axis.lo
${top_builddir}/diag_manager/diag_output.lo: ${top_srcdir}/diag_manager/diag_output.F90 ${top_srcdir}/platform/platform.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/diag_manager/diag_axis.lo ${top_srcdir}/diag_manager/diag_data.lo ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/fms2_io/fms2_io.lo ${top_srcdir}/axis_utils/axis_utils2.lo
${top_builddir}/diag_manager/diag_table.lo: ${top_srcdir}/diag_manager/diag_table.F90 ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/constants/constants.lo ${top_srcdir}/diag_manager/diag_data.lo ${top_srcdir}/diag_manager/diag_util.lo
${top_builddir}/diag_manager/diag_util.lo: ${top_srcdir}/diag_manager/diag_util.F90 ${top_srcdir}/platform/platform.lo ${top_srcdir}/diag_manager/diag_data.lo ${top_srcdir}/diag_manager/diag_axis.lo ${top_srcdir}/diag_manager/diag_output.lo ${top_srcdir}/diag_manager/diag_grid.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/fms/fms_io.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/fms/fms_io.lo ${top_srcdir}/constants/constants.lo ${top_srcdir}/fms2_io/fms2_io.lo
${top_builddir}/drifters/cloud_interpolator.lo: ${top_srcdir}/drifters/cloud_interpolator.F90
${top_builddir}/drifters/drifters.lo: ${top_srcdir}/drifters/drifters.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/drifters/drifters_core.lo ${top_srcdir}/drifters/drifters_input.lo ${top_srcdir}/drifters/drifters_io.lo ${top_srcdir}/drifters/drifters_comm.lo ${top_srcdir}/drifters/cloud_interpolator.lo
${top_builddir}/drifters/drifters_comm.lo: ${top_srcdir}/drifters/drifters_comm.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/drifters/drifters_core.lo ${top_srcdir}/drifters/drifters_input.lo
${top_builddir}/drifters/drifters_core.lo: ${top_srcdir}/drifters/drifters_core.F90
${top_builddir}/drifters/drifters_input.lo: ${top_srcdir}/drifters/drifters_input.F90
${top_builddir}/drifters/drifters_io.lo: ${top_srcdir}/drifters/drifters_io.F90
${top_builddir}/drifters/quicksort.lo: ${top_srcdir}/drifters/quicksort.F90
${top_builddir}/exchange/stock_constants.lo: ${top_srcdir}/exchange/stock_constants.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/diag_manager/diag_manager.lo
${top_builddir}/exchange/xgrid.lo: ${top_srcdir}/exchange/xgrid.F90 ${top_srcdir}/fms/fms.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/constants/constants.lo ${top_srcdir}/mosaic/mosaic2.lo ${top_srcdir}/exchange/stock_constants.lo ${top_srcdir}/mosaic/gradient.lo ${top_srcdir}/fms2_io/fms2_io.lo ${top_srcdir}/fms2_io/fms2_io.lo ${top_srcdir}/fms2_io/fms2_io.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/diag_manager/diag_manager.lo
${top_builddir}/field_manager/field_manager.lo: ${top_srcdir}/field_manager/field_manager.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/fms2_io/fms2_io.lo
${top_builddir}/field_manager/fm_util.lo: ${top_srcdir}/field_manager/fm_util.F90 ${top_srcdir}/field_manager/field_manager.lo ${top_srcdir}/field_manager/field_manager.lo ${top_srcdir}/field_manager/field_manager.lo ${top_srcdir}/field_manager/field_manager.lo ${top_srcdir}/field_manager/field_manager.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/mpp/mpp.lo
${top_builddir}/fms/fms.lo: ${top_srcdir}/fms/fms.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/fms/fms_io.lo ${top_srcdir}/fms2_io/fms2_io.lo ${top_srcdir}/memutils/memutils.lo ${top_srcdir}/constants/constants.lo
${top_builddir}/fms/fms_io.lo: ${top_srcdir}/fms/fms_io.F90 ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/platform/platform.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp_io.lo
${top_builddir}/fms2_io/blackboxio.lo: ${top_srcdir}/fms2_io/blackboxio.F90 ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/fms2_io/fms_io_utils.lo ${top_srcdir}/fms2_io/netcdf_io.lo ${top_srcdir}/fms2_io/fms_netcdf_domain_io.lo ${top_srcdir}/fms2_io/fms_netcdf_unstructured_domain_io.lo ${top_srcdir}/mpp/mpp.lo
${top_builddir}/fms2_io/fms2_io.lo: ${top_srcdir}/fms2_io/fms2_io.F90 ${top_srcdir}/fms2_io/fms_io_utils.lo ${top_srcdir}/fms2_io/netcdf_io.lo ${top_srcdir}/fms2_io/fms_netcdf_domain_io.lo ${top_srcdir}/fms2_io/fms_netcdf_unstructured_domain_io.lo ${top_srcdir}/fms2_io/blackboxio.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo
${top_builddir}/fms2_io/fms_io_utils.lo: ${top_srcdir}/fms2_io/fms_io_utils.F90 ${top_srcdir}/mpp/mpp.lo
${top_builddir}/fms2_io/fms_netcdf_domain_io.lo: ${top_srcdir}/fms2_io/fms_netcdf_domain_io.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/fms2_io/fms_io_utils.lo ${top_srcdir}/fms2_io/netcdf_io.lo
${top_builddir}/fms2_io/fms_netcdf_unstructured_domain_io.lo: ${top_srcdir}/fms2_io/fms_netcdf_unstructured_domain_io.F90 ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/fms2_io/fms_io_utils.lo ${top_srcdir}/fms2_io/netcdf_io.lo
${top_builddir}/fms2_io/netcdf_io.lo: ${top_srcdir}/fms2_io/netcdf_io.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/fms2_io/fms_io_utils.lo
${top_builddir}/horiz_interp/horiz_interp.lo: ${top_srcdir}/horiz_interp/horiz_interp.F90 ${top_srcdir}/fms/fms.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/constants/constants.lo ${top_srcdir}/horiz_interp/horiz_interp_type.lo ${top_srcdir}/horiz_interp/horiz_interp_type.lo ${top_srcdir}/horiz_interp/horiz_interp_conserve.lo ${top_srcdir}/horiz_interp/horiz_interp_conserve.lo ${top_srcdir}/horiz_interp/horiz_interp_bilinear.lo ${top_srcdir}/horiz_interp/horiz_interp_bilinear.lo ${top_srcdir}/horiz_interp/horiz_interp_bicubic.lo ${top_srcdir}/horiz_interp/horiz_interp_bicubic.lo ${top_srcdir}/horiz_interp/horiz_interp_spherical.lo ${top_srcdir}/horiz_interp/horiz_interp_spherical.lo
${top_builddir}/horiz_interp/horiz_interp_bicubic.lo: ${top_srcdir}/horiz_interp/horiz_interp_bicubic.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/horiz_interp/horiz_interp_type.lo ${top_srcdir}/constants/constants.lo
${top_builddir}/horiz_interp/horiz_interp_bilinear.lo: ${top_srcdir}/horiz_interp/horiz_interp_bilinear.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/constants/constants.lo ${top_srcdir}/horiz_interp/horiz_interp_type.lo
${top_builddir}/horiz_interp/horiz_interp_conserve.lo: ${top_srcdir}/horiz_interp/horiz_interp_conserve.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/fms/fms_io.lo ${top_srcdir}/constants/constants.lo ${top_srcdir}/horiz_interp/horiz_interp_type.lo
${top_builddir}/horiz_interp/horiz_interp_spherical.lo: ${top_srcdir}/horiz_interp/horiz_interp_spherical.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/constants/constants.lo ${top_srcdir}/horiz_interp/horiz_interp_type.lo
${top_builddir}/horiz_interp/horiz_interp_type.lo: ${top_srcdir}/horiz_interp/horiz_interp_type.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo
${top_builddir}/interpolator/interpolator.lo: ${top_srcdir}/interpolator/interpolator.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/diag_manager/diag_manager.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/fms2_io/fms2_io.lo ${top_srcdir}/horiz_interp/horiz_interp.lo ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/time_interp/time_interp.lo ${top_srcdir}/constants/constants.lo
${top_builddir}/memutils/memutils.lo: ${top_srcdir}/memutils/memutils.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_memutils.lo
${top_builddir}/monin_obukhov/monin_obukhov.lo: ${top_srcdir}/monin_obukhov/monin_obukhov.F90 ${top_srcdir}/constants/constants.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/monin_obukhov/monin_obukhov_kernel.lo
${top_builddir}/monin_obukhov/monin_obukhov_kernel.lo: ${top_srcdir}/monin_obukhov/monin_obukhov_kernel.F90
${top_builddir}/mosaic/gradient.lo: ${top_srcdir}/mosaic/gradient.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/constants/constants.lo
${top_builddir}/mosaic/grid.lo: ${top_srcdir}/mosaic/grid.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/constants/constants.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/fms/fms_io.lo ${top_srcdir}/mosaic/mosaic.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mosaic/mosaic.lo
${top_builddir}/mosaic/mosaic.lo: ${top_srcdir}/mosaic/mosaic.F90 ${top_srcdir}/fms/fms.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/fms/fms_io.lo ${top_srcdir}/constants/constants.lo
${top_builddir}/mosaic/mosaic2.lo: ${top_srcdir}/mosaic/mosaic2.F90 ${top_srcdir}/fms/fms.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/constants/constants.lo ${top_srcdir}/fms2_io/fms2_io.lo ${top_srcdir}/fms2_io/fms2_io.lo
${top_builddir}/mpp/mpp.lo: ${top_srcdir}/mpp/mpp.F90 ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_data.lo ${top_srcdir}/mpp/mpp_data.lo ${top_srcdir}/mpp/mpp_data.lo
${top_builddir}/mpp/mpp_data.lo: ${top_srcdir}/mpp/mpp_data.F90 ${top_srcdir}/mpp/mpp_parameter.lo
${top_builddir}/mpp/mpp_domains.lo: ${top_srcdir}/mpp/mpp_domains.F90 ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_data.lo ${top_srcdir}/mpp/mpp_data.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_memutils.lo ${top_srcdir}/mpp/mpp_pset.lo ${top_srcdir}/mpp/mpp_efp.lo
${top_builddir}/mpp/mpp_efp.lo: ${top_srcdir}/mpp/mpp_efp.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo
${top_builddir}/mpp/mpp_io.lo: ${top_srcdir}/mpp/mpp_io.F90 ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo
${top_builddir}/mpp/mpp_memutils.lo: ${top_srcdir}/mpp/mpp_memutils.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo
${top_builddir}/mpp/mpp_parameter.lo: ${top_srcdir}/mpp/mpp_parameter.F90
${top_builddir}/mpp/mpp_pset.lo: ${top_srcdir}/mpp/mpp_pset.F90 ${top_srcdir}/mpp/mpp.lo
${top_builddir}/mpp/mpp_utilities.lo: ${top_srcdir}/mpp/mpp_utilities.F90 ${top_srcdir}/mpp/mpp.lo
${top_builddir}/platform/platform.lo: ${top_srcdir}/platform/platform.F90
${top_builddir}/random_numbers/MersenneTwister.lo: ${top_srcdir}/random_numbers/MersenneTwister.F90
${top_builddir}/random_numbers/random_numbers.lo: ${top_srcdir}/random_numbers/random_numbers.F90 ${top_srcdir}/random_numbers/MersenneTwister.lo ${top_srcdir}/time_manager/time_manager.lo
${top_builddir}/sat_vapor_pres/sat_vapor_pres.lo: ${top_srcdir}/sat_vapor_pres/sat_vapor_pres.F90 ${top_srcdir}/constants/constants.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/sat_vapor_pres/sat_vapor_pres_k.lo
${top_builddir}/sat_vapor_pres/sat_vapor_pres_k.lo: ${top_srcdir}/sat_vapor_pres/sat_vapor_pres_k.F90
${top_builddir}/time_interp/time_interp.lo: ${top_srcdir}/time_interp/time_interp.F90 ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/mpp/mpp.lo
${top_builddir}/time_interp/time_interp_external.lo: ${top_srcdir}/time_interp/time_interp_external.F90 ${top_srcdir}/fms/fms.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/time_manager/get_cal_time.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/time_interp/time_interp.lo ${top_srcdir}/axis_utils/axis_utils.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/platform/platform.lo ${top_srcdir}/horiz_interp/horiz_interp.lo
${top_builddir}/time_interp/time_interp_external2.lo: ${top_srcdir}/time_interp/time_interp_external2.F90 ${top_srcdir}/fms/fms.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/time_manager/get_cal_time.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/time_interp/time_interp.lo ${top_srcdir}/axis_utils/axis_utils2.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/platform/platform.lo ${top_srcdir}/horiz_interp/horiz_interp.lo ${top_srcdir}/fms2_io/fms2_io.lo
${top_builddir}/time_manager/get_cal_time.lo: ${top_srcdir}/time_manager/get_cal_time.F90 ${top_srcdir}/fms/fms.lo ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/mpp/mpp.lo
${top_builddir}/time_manager/time_manager.lo: ${top_srcdir}/time_manager/time_manager.F90 ${top_srcdir}/constants/constants.lo ${top_srcdir}/fms/fms.lo
${top_builddir}/topography/gaussian_topog.lo: ${top_srcdir}/topography/gaussian_topog.F90 ${top_srcdir}/fms/fms.lo ${top_srcdir}/constants/constants.lo ${top_srcdir}/mpp/mpp.lo
${top_builddir}/topography/topography.lo: ${top_srcdir}/topography/topography.F90 ${top_srcdir}/topography/gaussian_topog.lo ${top_srcdir}/horiz_interp/horiz_interp.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/fms2_io/fms2_io.lo ${top_srcdir}/constants/constants.lo ${top_srcdir}/mpp/mpp.lo
${top_builddir}/tracer_manager/tracer_manager.lo: ${top_srcdir}/tracer_manager/tracer_manager.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/field_manager/field_manager.lo
${top_builddir}/tridiagonal/tridiagonal.lo: ${top_srcdir}/tridiagonal/tridiagonal.F90

CLEANFILES = *.mod




#    ${top_srcdir}/test_fms/affinity/test_affinity.F90 \
#    ${top_srcdir}/test_fms/axis_utils/test_axis_utils.F90 \
#    ${top_srcdir}/test_fms/data_override/test_data_override.F90 \
#    ${top_srcdir}/test_fms/diag_manager/test_diag_manager.F90 \
#    ${top_srcdir}/test_fms/drifters/test_cloud_interpolator.F90 \
#    ${top_srcdir}/test_fms/drifters/test_drifters.F90 \
#    ${top_srcdir}/test_fms/drifters/test_drifters_comm.F90 \
#    ${top_srcdir}/test_fms/drifters/test_drifters_core.F90 \
#    ${top_srcdir}/test_fms/drifters/test_drifters_input.F90 \
#    ${top_srcdir}/test_fms/drifters/test_drifters_io.F90 \
#    ${top_srcdir}/test_fms/drifters/test_quicksort.F90 \
#    ${top_srcdir}/test_fms/exchange/test_xgrid.F90 \
#    ${top_srcdir}/test_fms/field_manager/test_field_manager.F90 \
#    ${top_srcdir}/test_fms/fms/test_fms_io.F90 \
#    ${top_srcdir}/test_fms/fms/test_unstructured_fms_io.F90 \
#    ${top_srcdir}/test_fms/fms2_io/argparse.F90 \
#    ${top_srcdir}/test_fms/fms2_io/setup.F90 \
#    ${top_srcdir}/test_fms/fms2_io/test_atmosphere_io.F90 \
#    ${top_srcdir}/test_fms/fms2_io/test_fms2_io.F90 \
#    ${top_srcdir}/test_fms/horiz_interp/test2_horiz_interp.F90 \
#    ${top_srcdir}/test_fms/horiz_interp/test_horiz_interp.F90 \
#    ${top_srcdir}/test_fms/interpolator/test_interpolator.F90 \
#    ${top_srcdir}/test_fms/monin_obukhov/test_monin_obukhov.F90 \
#    ${top_srcdir}/test_fms/mosaic/test_mosaic.F90 \
#    ${top_srcdir}/test_fms/mpp/test_mpp.F90 \
#    ${top_srcdir}/test_fms/mpp/test_mpp_domains.F90 \
#    ${top_srcdir}/test_fms/mpp/test_mpp_mem_dump.F90 \
#    ${top_srcdir}/test_fms/mpp/test_mpp_memuse.F90 \
#    ${top_srcdir}/test_fms/mpp/test_mpp_memutils_begin_2x.F90 \
#    ${top_srcdir}/test_fms/mpp/test_mpp_memutils_begin_end.F90 \
#    ${top_srcdir}/test_fms/mpp/test_mpp_memutils_end_before_begin.F90 \
#    ${top_srcdir}/test_fms/mpp/test_mpp_print_memuse_stats_file.F90 \
#    ${top_srcdir}/test_fms/mpp/test_mpp_print_memuse_stats_stderr.F90 \
#    ${top_srcdir}/test_fms/mpp/test_mpp_pset.F90 \
#    ${top_srcdir}/test_fms/mpp_io/test_mpp_io.F90 \
#    ${top_srcdir}/test_fms/time_interp/test_time_interp.F90 \
#    ${top_srcdir}/test_fms/time_interp/test_time_interp_external.F90 \
#    ${top_srcdir}/test_fms/time_manager/test_time_manager.F90 \
#${top_builddir}/test_fms/affinity/test_affinity.lo: ${top_srcdir}/test_fms/affinity/test_affinity.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/affinity/fms_affinity.lo
#${top_builddir}/test_fms/axis_utils/test_axis_utils.lo: ${top_srcdir}/test_fms/axis_utils/test_axis_utils.F90 ${top_srcdir}/fms/fms.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/axis_utils/axis_utils.lo
#${top_builddir}/test_fms/data_override/test_data_override.lo: ${top_srcdir}/test_fms/data_override/test_data_override.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/affinity/fms_affinity.lo ${top_srcdir}/fms/fms_io.lo ${top_srcdir}/constants/constants.lo ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/diag_manager/diag_manager.lo ${top_srcdir}/diag_manager/diag_manager.lo ${top_srcdir}/data_override/data_override.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_memutils.lo
#${top_builddir}/test_fms/diag_manager/test_diag_manager.lo: ${top_srcdir}/test_fms/diag_manager/test_diag_manager.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/fms/fms_io.lo ${top_srcdir}/fms/fms_io.lo ${top_srcdir}/constants/constants.lo ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/diag_manager/diag_manager.lo ${top_srcdir}/diag_manager/diag_manager.lo ${top_srcdir}/diag_manager/diag_manager.lo ${top_srcdir}/diag_manager/diag_manager.lo ${top_srcdir}/diag_manager/diag_manager.lo ${top_srcdir}/diag_manager/diag_axis.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/diag_manager/diag_axis.lo ${top_srcdir}/diag_manager/diag_manager.lo ${top_srcdir}/time_manager/time_manager.lo
#${top_builddir}/test_fms/drifters/test_cloud_interpolator.lo: ${top_srcdir}/test_fms/drifters/test_cloud_interpolator.F90 ${top_srcdir}/drifters/cloud_interpolator.lo ${top_srcdir}/mpp/mpp.lo
#${top_builddir}/test_fms/drifters/test_drifters.lo: ${top_srcdir}/test_fms/drifters/test_drifters.F90 ${top_srcdir}/drifters/drifters.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp.lo
#${top_builddir}/test_fms/drifters/test_drifters_comm.lo: ${top_srcdir}/test_fms/drifters/test_drifters_comm.F90 ${top_srcdir}/drifters/drifters_core.lo ${top_srcdir}/drifters/drifters_comm.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo
#${top_builddir}/test_fms/drifters/test_drifters_core.lo: ${top_srcdir}/test_fms/drifters/test_drifters_core.F90 ${top_srcdir}/drifters/drifters_core.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/mpp/mpp.lo
#${top_builddir}/test_fms/drifters/test_drifters_input.lo: ${top_srcdir}/test_fms/drifters/test_drifters_input.F90 ${top_srcdir}/drifters/drifters_input.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/mpp/mpp.lo
#${top_builddir}/test_fms/drifters/test_drifters_io.lo: ${top_srcdir}/test_fms/drifters/test_drifters_io.F90 ${top_srcdir}/drifters/drifters_io.lo ${top_srcdir}/mpp/mpp.lo
#${top_builddir}/test_fms/drifters/test_quicksort.lo: ${top_srcdir}/test_fms/drifters/test_quicksort.F90
#${top_builddir}/test_fms/exchange/test_xgrid.lo: ${top_srcdir}/test_fms/exchange/test_xgrid.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/fms/fms_io.lo ${top_srcdir}/fms2_io/fms2_io.lo ${top_srcdir}/fms2_io/fms2_io.lo ${top_srcdir}/fms2_io/fms2_io.lo ${top_srcdir}/constants/constants.lo ${top_srcdir}/exchange/xgrid.lo ${top_srcdir}/exchange/xgrid.lo ${top_srcdir}/exchange/xgrid.lo ${top_srcdir}/exchange/xgrid.lo ${top_srcdir}/mosaic/mosaic2.lo ${top_srcdir}/mosaic/mosaic2.lo ${top_srcdir}/mosaic/grid.lo ${top_srcdir}/mosaic/gradient.lo ${top_srcdir}/coupler/ensemble_manager.lo ${top_srcdir}/coupler/ensemble_manager.lo
#${top_builddir}/test_fms/field_manager/test_field_manager.lo: ${top_srcdir}/test_fms/field_manager/test_field_manager.F90 ${top_srcdir}/field_manager/field_manager.lo ${top_srcdir}/mpp/mpp.lo
#${top_builddir}/test_fms/fms/test_fms_io.lo: ${top_srcdir}/test_fms/fms/test_fms_io.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/fms/fms_io.lo ${top_srcdir}/fms/fms_io.lo ${top_srcdir}/fms/fms_io.lo ${top_srcdir}/mpp/mpp_io.lo
#${top_builddir}/test_fms/fms/test_unstructured_fms_io.lo: ${top_srcdir}/test_fms/fms/test_unstructured_fms_io.F90 ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/fms/fms_io.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_parameter.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/fms/fms_io.lo
#${top_builddir}/test_fms/fms2_io/argparse.lo: ${top_srcdir}/test_fms/fms2_io/argparse.F90 ${top_srcdir}/test_fms/fms2_io/argparse.lo
#${top_builddir}/test_fms/fms2_io/setup.lo: ${top_srcdir}/test_fms/fms2_io/setup.F90 ${top_srcdir}/test_fms/fms2_io/argparse.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo
#${top_builddir}/test_fms/fms2_io/test_atmosphere_io.lo: ${top_srcdir}/test_fms/fms2_io/test_atmosphere_io.F90 ${top_srcdir}/fms2_io/fms2_io.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/test_fms/fms2_io/setup.lo
#${top_builddir}/test_fms/fms2_io/test_fms2_io.lo: ${top_srcdir}/test_fms/fms2_io/test_fms2_io.F90 ${top_srcdir}/test_fms/fms2_io/argparse.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/fms2_io/fms2_io.lo
#${top_builddir}/test_fms/horiz_interp/test2_horiz_interp.lo: ${top_srcdir}/test_fms/horiz_interp/test2_horiz_interp.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mosaic/mosaic.lo ${top_srcdir}/horiz_interp/horiz_interp.lo ${top_srcdir}/horiz_interp/horiz_interp.lo ${top_srcdir}/axis_utils/axis_utils.lo ${top_srcdir}/fms/fms_io.lo ${top_srcdir}/fms/fms_io.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/constants/constants.lo ${top_srcdir}/horiz_interp/horiz_interp.lo ${top_srcdir}/axis_utils/axis_utils.lo
#${top_builddir}/test_fms/horiz_interp/test_horiz_interp.lo: ${top_srcdir}/test_fms/horiz_interp/test_horiz_interp.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/horiz_interp/horiz_interp.lo ${top_srcdir}/horiz_interp/horiz_interp.lo ${top_srcdir}/constants/constants.lo
#${top_builddir}/test_fms/interpolator/test_interpolator.lo: ${top_srcdir}/test_fms/interpolator/test_interpolator.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/diag_manager/diag_manager.lo ${top_srcdir}/interpolator/interpolator.lo ${top_srcdir}/constants/constants.lo ${top_srcdir}/time_interp/time_interp.lo
#${top_builddir}/test_fms/monin_obukhov/test_monin_obukhov.lo: ${top_srcdir}/test_fms/monin_obukhov/test_monin_obukhov.F90 ${top_srcdir}/monin_obukhov/monin_obukhov_kernel.lo ${top_srcdir}/mpp/mpp.lo
#${top_builddir}/test_fms/mosaic/test_mosaic.lo: ${top_srcdir}/test_fms/mosaic/test_mosaic.F90 ${top_srcdir}/mosaic/mosaic.lo ${top_srcdir}/mosaic/mosaic.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/fms/fms_io.lo ${top_srcdir}/fms/fms.lo
#${top_builddir}/test_fms/mpp/test_mpp.lo: ${top_srcdir}/test_fms/mpp/test_mpp.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp.lo
#${top_builddir}/test_fms/mpp/test_mpp_domains.lo: ${top_srcdir}/test_fms/mpp/test_mpp_domains.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_memutils.lo ${top_srcdir}/affinity/fms_affinity.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo
#${top_builddir}/test_fms/mpp/test_mpp_mem_dump.lo: ${top_srcdir}/test_fms/mpp/test_mpp_mem_dump.F90 ${top_srcdir}/mpp/mpp_memutils.lo
#${top_builddir}/test_fms/mpp/test_mpp_memuse.lo: ${top_srcdir}/test_fms/mpp/test_mpp_memuse.F90
#${top_builddir}/test_fms/mpp/test_mpp_memutils_begin_2x.lo: ${top_srcdir}/test_fms/mpp/test_mpp_memutils_begin_2x.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_memutils.lo
#${top_builddir}/test_fms/mpp/test_mpp_memutils_begin_end.lo: ${top_srcdir}/test_fms/mpp/test_mpp_memutils_begin_end.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_memutils.lo
#${top_builddir}/test_fms/mpp/test_mpp_memutils_end_before_begin.lo: ${top_srcdir}/test_fms/mpp/test_mpp_memutils_end_before_begin.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_memutils.lo
#${top_builddir}/test_fms/mpp/test_mpp_print_memuse_stats_file.lo: ${top_srcdir}/test_fms/mpp/test_mpp_print_memuse_stats_file.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_memutils.lo ${top_srcdir}/mpp/mpp_memutils.lo
#${top_builddir}/test_fms/mpp/test_mpp_print_memuse_stats_stderr.lo: ${top_srcdir}/test_fms/mpp/test_mpp_print_memuse_stats_stderr.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_memutils.lo ${top_srcdir}/mpp/mpp_memutils.lo
#${top_builddir}/test_fms/mpp/test_mpp_pset.lo: ${top_srcdir}/test_fms/mpp/test_mpp_pset.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_pset.lo
#${top_builddir}/test_fms/mpp_io/test_mpp_io.lo: ${top_srcdir}/test_fms/mpp_io/test_mpp_io.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp.lo
#${top_builddir}/test_fms/time_interp/test_time_interp.lo: ${top_srcdir}/test_fms/time_interp/test_time_interp.F90 ${top_srcdir}/fms/fms.lo ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/time_interp/time_interp.lo
#${top_builddir}/test_fms/time_interp/test_time_interp_external.lo: ${top_srcdir}/test_fms/time_interp/test_time_interp_external.F90 ${top_srcdir}/constants/constants.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/mpp/mpp_io.lo ${top_srcdir}/mpp/mpp_domains.lo ${top_srcdir}/time_interp/time_interp_external.lo ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/horiz_interp/horiz_interp.lo ${top_srcdir}/axis_utils/axis_utils.lo
#${top_builddir}/test_fms/time_manager/test_time_manager.lo: ${top_srcdir}/test_fms/time_manager/test_time_manager.F90 ${top_srcdir}/mpp/mpp.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/fms/fms.lo ${top_srcdir}/constants/constants.lo ${top_srcdir}/fms/fms_io.lo ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/time_manager/time_manager.lo ${top_srcdir}/time_manager/time_manager.lo
