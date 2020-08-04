# Volk Profile for Android

This is an Android app that profiles the [Volk library](https://www.libvolk.org/) and selects the best implementations for the architecture, similar to `volk_profile` that comes as part of Volk. It is configured to run in *benchmark* mode, i.e., it will profile all kernels, even the ones without alternative implementations. Running the application generates two files:
- `volk_results.json` with benchmark results
- `volk_config` with the best implementations for the architecture

Files are stored on *External Storage* (which might be an SD card or some internal storage). The directory can be accessed when the phone is connected to a PC and is available to all applications on the phone that have the corresponding permission. This way apps that use Volk (like GNU Radio Android apps) can use it to select the implementations that perform best.

Volk comes with a [basic evaluation script](https://github.com/gnuradio/volk/blob/master/apps/plot_best_vs_generic.py) for the JSON file.

## Installation

Building the app requires the [GNU Radio Android toolchain](https://github.com/bastibl/gnuradio-android/). Please see this repository for further instructions on how to build the toolchain and apps that use it.


## Running the App

The app is minimalistic. There is only a button to start profiling and a textbox that indicates when it is done. You can use the `logcat` output to monitor progress.


## Exemplary Results

The following plot shows the speedup of unaligned implementations on a OnePlus 5T smartphone.

![Volk Unaligned Implementations](doc/volk_unaligned.png)

## Publication

If you use this project, we would appreciate a reference to:

<ul>
<li>
<a href="http://dx.doi.org/10.1145/3411276.3412184"><img src="https://www.bastibl.net/bib/icons/ACM-logo.gif" title="ACM" alt=""></a> <a class="bibauthorlink" href="https://www.bastibl.net/">Bastian Bloessl</a>, Lars Baumgärtner and Matthias Hollick, “<strong>Hardware-Accelerated Real-Time Stream Data Processing on Android with GNU Radio</strong>,” Proceedings of 14th International Workshop on Wireless Network Testbeds, Experimental evaluation &amp; Characterization (WiNTECH’20), London, UK, September 2020.
 <small>[<a href="http://dx.doi.org/10.1145/3411276.3412184">DOI</a>, <a href="https://www.bastibl.net/bib/bloessl2020hardware/bloessl2020hardware.bib">BibTeX</a>, <a href="https://www.bastibl.net/bib/bloessl2020hardware/">PDF and Details…</a>]</small></p>
</li>
</ul>
