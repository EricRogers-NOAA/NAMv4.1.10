SHELL=/bin/bash
# All libraries should be set in the environment prior to compiling

SUBDIRS0 = nam_pcpbudget.fd nam_sndp.fd nam_firewx_gridspecs.fd \
           nam_lspasno24h.fd nam_makeprecip.fd nam_merge2n4.fd \
           nam_pcp_conus.fd nam_pcpprep.fd nam_post0.fd nam_icwf.fd \
           nam_ref2nemsio.fd nam_stage4_acc.fd nam_stnmlist.fd \
           nam_anl_step_b.fd nam_create_trak_fnl.fd nam_gettrk.fd \
           nam_inter_4to6.fd nam_merge_final.fd nam_merge_nest_b.fd \
           nam_nems_bin_io.fd nam_rem_specified_prepbufr.fd \
           nam_split_b.fd nam_supvit.fd nam_latlon_to_ij_h_bgrid_tc.fd

SUBDIRS1 = nam_get1236hpcp.fd

all: $(SUBDIRS0) $(SUBDIRS1)
	@for dir in $(SUBDIRS0); do \
           ( cd $$dir; echo "Making $@ in `pwd`" ; make ) ; \
        done
	@for dir in $(SUBDIRS1); do \
           ( cd $$dir; echo "Making $@ in `pwd`" ; ./compile ) ; \
        done
clean: $(SUBDIRS0)
	@for dir in $(SUBDIRS0); do \
          ( cd $$dir; echo "Cleaning `pwd`" ; make $@ ) ; \
	done
