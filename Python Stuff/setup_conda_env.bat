ECHO "Setting Up Conda Environment For qutau"

conda env create -f conda_env.yml

ECHO "To run:"
ECHO "Activate conda with 'conda activate qutau'"
ECHO "Then run 'python'"
ECHO "In python 'import qutau_measurement'"
ECHO "Then run with 'qutau_measurement.run_experiment(runtime, channel_trig, channel_spcm, max_time, bin_width, fake_data=False, absolute_plotting=False)'"
PAUSE