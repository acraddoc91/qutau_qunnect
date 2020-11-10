import ctypes
import _ctypes
from matplotlib import pyplot as plt
import numpy as np
import time
import os
import tkinter as tk
import tkinter.filedialog
import h5py

tagger_resolution = 82.3e-12*2

working_directory = os.getcwd()
lib_name = 'tdcbase.dll'
extra_lib_name = 'qutau_dll.dll'

## Accumulate data for defined amount of time on the trigger and spcm channels
def take_data(runtime, channel_trig, channel_spcm):
    
    base = ctypes.CDLL(working_directory + '/' + lib_name)

    BUFFER_SIZE = 1000000
    
    # We're going to use these in a bit
    timestamps = np.array([])
    channels = np.array([])
    
    # For the ctypes stuff
    timestamps_dummy = (ctypes.c_longlong*BUFFER_SIZE)()
    channels_dummy = (ctypes.c_byte*BUFFER_SIZE)()
    tags_valid = ctypes.c_int(0)
    
    # Needs to be run to initialise the TDC, will pick up any TDC
    base.TDC_init(-1)
    
    # Figure out what channels need enabling and tell the TDC to do that thing
    channels_to_enable = (1 << (channel_trig-1)) + (1 << (channel_spcm-1))
    base.TDC_enableChannels(channels_to_enable)
    
    # Set the size of the buffer for the TDC, set to max here
    base.TDC_setTimestampBufferSize( BUFFER_SIZE )
    
    # Loop to accumulate tags
    start = time.time()
    while time.time() - start < runtime:
        # Wait to accumualte tags
        time.sleep(0.05)
        
        # Grab tags and stick them into the arrays we've made
        base.TDC_getLastTimestamps( 1, timestamps_dummy, channels_dummy, ctypes.byref(tags_valid))
        timestamps = np.append(timestamps,timestamps_dummy[0:tags_valid.value])
        channels = np.append(channels,channels_dummy[0:tags_valid.value])
    
    # Release TDC
    base.TDC_deInit()

    # Free the library, do this or windows loses it's shit
    _ctypes.FreeLibrary(base._handle)

    # Return tags and channels
    return timestamps, channels

# For testing purposes, to generate data
def take_fake_data(runtime, channel_trig, channel_spcm):

    base = ctypes.CDLL(working_directory + '/' + lib_name)
    extra = ctypes.CDLL(working_directory + '/' + extra_lib_name)

    BUFFER_SIZE = 1000000

    # We're going to use these in a bit
    timestamps = np.array([])
    channels = np.array([])

    # For the ctypes stuff
    timestamps_dummy = (ctypes.c_longlong*BUFFER_SIZE)()
    channels_dummy = (ctypes.c_byte*BUFFER_SIZE)()
    tags_valid = ctypes.c_int(0)

    # Needs to be run to initialise the TDC, will pick up any TDC
    base.TDC_init(-1)
    
    # Figure out what channels need enabling and tell the TDC to do that thing
    channels_to_enable = (1 << (channel_trig-1)) + (1 << (channel_spcm-1))
    base.TDC_enableChannels(channels_to_enable)

    # Set the size of the buffer for the TDC, set to max here
    base.TDC_setTimestampBufferSize( BUFFER_SIZE )

    # This is the simulated data distribution
    simPar = (ctypes.c_double*2)()
    # Let's go with a mean value of 5us between events
    simPar[0] = 32000
    simPar[1] = 32000

    # Figure out how many tags gets us to around the runtime
    num_loops = int(np.round(runtime/(5e-6*1000)))
    for _ in range(0,num_loops):
        # Generate tags
        extra.TDC_generateTimestamps_flat(simPar,1000)

        # Grab tags and stick them into the arrays we've made
        base.TDC_getLastTimestamps( 1, timestamps_dummy, channels_dummy, ctypes.byref(tags_valid))
        timestamps = np.append(timestamps,timestamps_dummy[0:tags_valid.value])
        channels = np.append(channels,channels_dummy[0:tags_valid.value])
    
    # Release TDC
    base.TDC_deInit()

    # Free the library, do this or windows loses it's shit
    _ctypes.FreeLibrary(base._handle)
    _ctypes.FreeLibrary(extra._handle)

    # Return tags and channels
    return timestamps, channels

## Get SPCM profile relative to trigger
def get_profile(max_time, bin_width, tags_trig, tags_spcm):
    
    # Make sure bin_width is an integer multiple of the tagger bin width
    # tagger_width = 81e-12
    tagger_width = 88.18e-12

    # Factor of 2 because there's some stupid issue with TDC chip with asymmetric bins
    act_bin_width_int = float(2*np.round(bin_width/(2*tagger_width)))
    max_time_int = int(np.round(max_time/(act_bin_width_int*tagger_width)))
    
    
    # Rescale the tags
    rescaled_trig = np.round(tags_trig.astype(float)/act_bin_width_int).astype(int)
    rescaled_spcm = np.round(tags_spcm.astype(float)/act_bin_width_int).astype(int)
    
    # Vectors to hold coincidences
    coinc = np.zeros(max_time_int+1)
    tau = np.array(range(len(coinc)))*act_bin_width_int*tagger_width
    
    # Define a pointers that indicate elements from the SPCM vector which are in range
    lower_pointer = 0
    upper_pointer = 0
    
    # Loop through trigger tags
    for k in range(len(rescaled_trig)):
        
        # Find the first tag from the SPCM vector which is equal or larger than the current trigger tag
        going = True
        j = lower_pointer
        while going:
            if rescaled_trig[k] > rescaled_spcm[j]:
                j += 1
                lower_pointer = j
            else:
                going = False
                lower_pointer = j
            
            # Make sure we don't go out of range on the spcm vector
            if j >= len(rescaled_spcm):
                going = False
                lower_pointer = j
            
        # Find the last tag in the SPCM vector which is smaller than the trigger plus the max time 
        going = True
        j = upper_pointer
        # Make sure that j can't be negative
        if j < 0:
            j = 0
        while going:
            if rescaled_trig[k] + max_time_int >= rescaled_spcm[j]:
                j += 1
                upper_pointer = j
            else:
                going = False
                upper_pointer = j - 1
            
            # Make sure we don't go out of range on the spcm vector. Minus one here so that upper_pointer < lower_pointer
            # for out of range elements
            if j >= len(rescaled_spcm):
                going = False
                upper_pointer = len(rescaled_spcm)-1
            
        # If we reach the end of the vector on the lower pointer we might as well stop
        if lower_pointer >= len(rescaled_spcm):
            break
        
        # Figure out when the coincidences were
        if upper_pointer>=lower_pointer:
            for l in range(lower_pointer,upper_pointer + 1):
                coinc_index = int(rescaled_spcm[l] - rescaled_trig[k])
                coinc[coinc_index] += 1
        
    return tau, coinc

## This acts like a wrapper to the profile stuff, but breaks stuff out into channels
def get_histogram(timestamps, channels, channel_trig, channel_spcm, max_time, bin_width):
    
    # Sort the tags so we know what came from what
    tags_trig = np.array(timestamps)[np.array(channels) == (channel_trig-1)]
    tags_spcm = np.array(timestamps)[np.array(channels) == (channel_spcm-1)]
    
    tau, coinc = get_profile(max_time, bin_width, tags_trig, tags_spcm)
    return tau, coinc

def yes_or_no(question):
    while "the answer is invalid":
        reply = str(input(question+' (y/n): ')).lower().strip()
        if reply[:1] == 'y':
            return True
        if reply[:1] == 'n':
            return False

# Dialog for getting notes to save to file
class MyDialog:

    def __init__(self, parent):

        top = self.top = tk.Toplevel(parent)

        #tk.Label(top, text="Value").pack()

        self.e = tk.Text(top)
        self.e.pack(padx=5)

        b = tk.Button(top, text="Save Notes", command=self.ok)
        b.pack(pady=5)

    def ok(self):

        self.result = self.e.get("1.0",'end-1c')
        self.top.destroy()

## For Jupyter
def run_experiment_jupyter(runtime, channel_trig, channel_spcm, max_time, bin_width, fake_data=False):
    # First let's grab the data
    if fake_data:
        timestamps, channels = take_fake_data(runtime, channel_trig, channel_spcm)
    else:
        timestamps, channels = take_data(runtime, channel_trig, channel_spcm)

    print("Collected data")
    # Then get coincidence profile
    tau, coinc = get_histogram(timestamps, channels, channel_trig, channel_spcm, max_time, bin_width)
    print("Processed coincs")

    return tau,coinc
    
## Runs the experiment
def run_experiment(runtime, channel_trig, channel_spcm, max_time, bin_width, fake_data=False, absolute_plotting=False):

    # First let's grab the data
    if fake_data:
        timestamps, channels = take_fake_data(runtime, channel_trig, channel_spcm)
    else:
        timestamps, channels = take_data(runtime, channel_trig, channel_spcm)
    
    # Then get coincidence profile
    tau, coinc = get_histogram(timestamps, channels, channel_trig, channel_spcm, max_time, bin_width)
    
    # Plot that shit
    if absolute_plotting:
        
        plt.plot(tau*1e6,coinc)
        plt.xlabel('t (us)')
        plt.ylabel('Absolute Coincidences')
        plt.show()
    else:
        # Rescale by number of triggers
        num_triggers = len(np.array(timestamps)[np.array(channels) == (channel_trig-1)])
        plt.plot(tau*1e6,coinc.astype(float)/float(num_triggers))
        plt.xlabel('t (us)')
        plt.ylabel('Coincidences Per Trigger')
        plt.show()
    
    # Ask if it wants saving to file
    save_to_file = yes_or_no("Save to file?")

    # Get filename and notes and save to file
    if save_to_file:

        # Run this or weird shit happens
        tk_window = tk.Tk()
        tk_window.withdraw()

        # Get filename and any notes
        filename = tkinter.filedialog.asksaveasfilename(filetypes=(("HDF5","*.h5"),))
        d = MyDialog(tk_window)
        tk_window.wait_window(d.top)
        notes = d.result

        # Make sure filename ends with the HDF5 type
        if filename[-3:] != ".h5":
            filename = filename + ".h5"

        # Save to HDF5
        h5f = h5py.File(filename, "w")
        h5f.create_group("Raw_Data")
        h5f.create_group("Information")
        h5f.create_group("Processed_Data")
        # Saving raw data in case we wanna re-process
        h5f.create_dataset("Raw_Data/tags",data=timestamps,dtype='i8')
        h5f.create_dataset("Raw_Data/channels",data=channels,dtype='i1')
        h5f.create_dataset("Information/channel_trig",data=channel_trig,dtype='i1')
        h5f.create_dataset("Information/channel_spcm",data=channel_spcm,dtype='i1')
        h5f.create_dataset("Information/Notes",data=notes)
        # Save histogram too, may not be necessary
        h5f.create_dataset("Processed_Data/time",data=tau,dtype='f8')
        h5f.create_dataset("Processed_Data/coincidences",data=coinc,dtype='i8')
        h5f.close()