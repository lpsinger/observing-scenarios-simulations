RUNS = O3 O4a O4 O5a O5b O5c
POPS = fullpop4
FILENAMES = events events.xml.gz events.sqlite injections.dat coincs.dat

all: psds injections public-alerts.dat

psds: $(foreach run,$(RUNS),runs/$(run)/psds.xml)

injections: $(foreach run,$(RUNS),$(foreach pop,$(POPS),$(foreach filename,$(FILENAMES),runs/$(run)/$(pop)/$(filename))))

.PHONY: all psds injections

#
# Tabulate O3 public alerts.
#

public-alerts.dat:
	scripts/get-public-alerts.py


#
# PSD files to use.
#

O3-psds = \
	--H1 aligo_O3actual_H1.txt \
	--L1 aligo_O3actual_L1.txt \
	--V1 avirgo_O3actual.txt
O4a-psds = \
	--H1 aligo_O4low.txt \
	--L1 aligo_O4low.txt
O4-psds = \
	--H1 aligo_O4low.txt \
	--L1 aligo_O4low.txt \
	--V1 O4_Virgo_78.txt \
	--K1 kagra_3Mpc.txt
O5a-psds = \
	--H1 O5StrainCurves_freqabc.txt --H1-column O5aStrain \
	--L1 O5StrainCurves_freqabc.txt --L1-column O5aStrain \
	--V1 avirgo_O5low_NEW.txt \
	--K1 kagra_128Mpc.txt
O5b-psds = \
	--H1 O5StrainCurves_freqabc.txt --H1-column O5bStrain \
	--L1 O5StrainCurves_freqabc.txt --L1-column O5bStrain \
	--V1 avirgo_O5low_NEW.txt \
	--K1 kagra_128Mpc.txt
O5c-psds = \
	--H1 O5StrainCurves_freqabc.txt --H1-column O5cStrain \
	--L1 O5StrainCurves_freqabc.txt --L1-column O5cStrain \
	--V1 avirgo_O5low_NEW.txt \
	--K1 kagra_128Mpc.txt
O6-psds = \
	--H1 AplusDesign.txt \
	--L1 AplusDesign.txt \
	--I1 AplusDesign.txt \
	--V1 avirgo_O5high_NEW.txt \
	--K1 kagra_128Mpc.txt


#
# Download official PSD files from the LIGO DCC.
#

%.txt:
	curl -OL https://dcc.ligo.org/LIGO-T2200043-v3/public/$(@F)

O5StrainCurves_freqabc.txt:
	curl -OL https://dcc.ligo.org/LIGO-T2500310-v2/public/$(@F)


#
# Pack PSDs into XML files for input to BAYESTAR.
#

.SECONDEXPANSION:
psd_files = $(sort $(filter %.txt,$(value $(1)-psds)))
runs/%/psds.xml: $$(call psd_files,%)
	mkdir -p $(@D) && scripts/pack-psds.py -o $@ $(value $*-psds)


#
# Download population parameters from the GWTC-4 Rates & Populations paper
# and draw intrinsic samples suitable for bayestar-inject.
#

analyses_AllCBC.tar:
	curl -OL https://zenodo.org/records/16911563/files/$(@F)

AllCBC_FullPop.h5: analyses_AllCBC.tar
	tar -xf $< $@

#
# Convert the Farah samples to the format needed by bayestar-inject.
#

fullpop4.h5: scripts/fullpop4.py AllCBC_FullPop.h5
	$^ $@


#
# Generate astrophysical distribution.
#

runs/%/fullpop4/injections.xml: $$(dir $$(@D))psds.xml fullpop4.h5
	mkdir -p $(@D) && cd $(@D) && bayestar-inject -l error --seed 1 -o $(@F) -j \
	--snr-threshold 1 --distribution-samples ../../../fullpop4.h5 --reference-psd ../psds.xml \
	--min-triggers 1 --nsamples 1000000

runs/%/injections.xml: $$(dir $$(@D))psds.xml
	mkdir -p $(@D) && cd $(@D) && bayestar-inject -l error --seed 1 -o $(@F) -j \
	--snr-threshold 1 --distribution $(notdir $(@D)) --reference-psd ../psds.xml \
	--min-triggers 1 --nsamples 1000000


#
# Generate and detect simulated signals.
#

runs/%/events.xml.gz: runs/%/injections.xml $$(dir $$(@D))psds.xml
	mkdir -p $(@D) && cd $(@D) && bayestar-realize-coincs \
	--seed 1 -j -l error -o $(@F) $(<F) \
	--reference-psd ../psds.xml \
	--snr-threshold 1 \
	--net-snr-threshold $(if $(filter bbh_astro,$(word 2,$(subst /, ,$*))),9,8) \
	--min-triggers 1 \
	--duty-cycle 0.7 --keep-subthreshold --measurement-error gaussian-noise \
	--detector $(subst --,,$(filter --%,$(value $(firstword $(subst /, ,$*))-psds)))


#
# Split detected signals file into individual files, one file per event.
#

runs/%/events: runs/%/events.xml.gz
	scripts/split-events.py $< $@

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
