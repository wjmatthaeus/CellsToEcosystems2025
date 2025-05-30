# recursive makefile for Biome-BGC
# 1) BIOME-BGC core science library
# 2) pointbgc executable for single-point, single-biome BIOME-BGC simulations
#
# invoke by issuing command "make" from this directory
#

ROOTDIR=$(pwd)

CFLAGS_GENERIC = -I $ROOTDIR
LDFLAGS_GENERIC = -lm

CFLAGS = ${CFLAGS_GENERIC}
LDFLAGS = ${LDFLAGS_GENERIC}
CC = gcc

OBJS = bgc.o output_map_init.o make_zero_flux_struct.o atm_pres.o \
	prephenology.o restart_io.o firstday.o zero_srcsnk.o daymet.o \
	dayphen.o phenology.o radtrans.o prcp_route.o snowmelt.o baresoil_evap.o \
	soilpsi.o maint_resp.o canopy_et.o photosynthesis.o outflow.o decomp.o \
	daily_allocation.o annual_rates.o growth_resp.o state_update.o \
	nleaching.o mortality.o check_balance.o summary.o smooth.o \
	precision_control.o bgc_io.o output_ascii.o get_co2.o get_ndep.o \
	pointbgc.o met_init.o restart_init.o time_init.o scc_init.o co2_init.o \
	sitec_init.o epc_init.o state_init.o output_init.o metarr_init.o \
	presim_state_init.o ramp_ndep_init.o output_ctrl.o ndep_init.o \
	end_init.o ini.o state_init.o

INCLUDE = bgc_struct.h bgc_func.h bgc_constants.h ini.h bgc_struct.h \
	pointbgc_struct.h pointbgc_func.h misc_func.h bgc_io.h

%.o: %.c $(INCLUDE)
	$(CC) -c -o $@ $< $(CFLAGS)


bgc : ${OBJS}
	${CC} -o $@ $^ ${CFLAGS} -o $@ ${LDFLAGS}


