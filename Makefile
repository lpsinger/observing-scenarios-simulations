RUNS = O3 O4 O5
POPS = bns_astro nsbh_astro bbh_astro
FILENAMES = events.xml.gz events.sqlite injections.dat coincs.dat
all: $(foreach run,$(RUNS),runs/$(run)/psds.xml) $(foreach run,$(RUNS),$(foreach pop,$(POPS),$(foreach filename,$(FILENAMES),runs/$(run)/$(pop)/$(filename))))


#
# PSD files to use from LIGO-T2000012.
#

O3-psds = \
	--H1 aligo_O3actual_H1.txt \
	--L1 aligo_O3actual_L1.txt \
	--V1 avirgo_O3actual.txt
O4-psds = \
	--H1 aligo_O4high.txt \
	--L1 aligo_O4high.txt \
	--V1 avirgo_O4high_NEW.txt \
	--K1 kagra_25Mpc.txt
O5-psds = \
	--H1 AplusDesign.txt \
	--L1 AplusDesign.txt \
	--V1 avirgo_O5high_NEW.txt \
	--K1 kagra_128Mpc.txt


#
# Download official PSD files from the LIGO DCC.
#

%.txt:
	curl -L https://dcc.ligo.org/LIGO-T2000012/public/$(@F) > $@


#
# Pack PSDs into XML files for input to BAYESTAR.
#

.SECONDEXPANSION:
psd_files = $(sort $(filter-out --%,$(value $(1)-psds)))
runs/%/psds.xml: $$(call psd_files,%)
	mkdir -p $(@D) && ./pack-psds.py -o $@ $(value $*-psds)


#
# Generate astrophysical distribution.
#

runs/%/injections.xml: $$(dir $$(@D))psds.xml
	mkdir -p $(@D) && cd $(@D) && bayestar-inject -l error --seed 1 -o $(@F) -j \
	--distribution $(notdir $(@D)) --reference-psd ../psds.xml


#
# Generate and detect simulated signals.
#

runs/%/events.xml.gz: runs/%/injections.xml $$(dir $$(@D))psds.xml
	mkdir -p $(@D) && cd $(@D) && bayestar-realize-coincs \
	--seed 1 -j -l error -o $(@F) $(<F) \
	--reference-psd ../psds.xml \
	--duty-cycle 0.7 --keep-subthreshold --measurement-error gaussian-noise \
	--detector $(subst --,,$(filter --%,$(value $(word 1,$(subst /, ,$*))-psds)))


#
# Convert simulated signals to SQLite format.
#

runs/%/events.sqlite: runs/%/events.xml.gz
	ligolw_sqlite -p -r -d $@ $<


#
# Generate convenient ASCII tables.
#

space := $(empty) $(empty)
tab := $(empty)	$(empty)
injections_dat_columns := simulation_id longitude latitude inclination distance mass1 mass2 spin1z spin2z
coincs_dat_columns := coinc_event_id ifos snr

%/injections.dat: %/events.xml.gz
	echo "$(subst $(space),$(tab),$(injections_dat_columns))" > $@ && \
	ligolw_print -t sim_inspiral $(injections_dat_columns:%=-c %) -d "$(tab)" $< >> $@

%/coincs.dat: %/events.xml.gz
	echo "$(subst $(space),$(tab),$(coincs_dat_columns))" > $@ && \
	ligolw_print -t coinc_inspiral $(coincs_dat_columns:%=-c %) -d "$(tab)" $< >> $@
