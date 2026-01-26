# conda activate your kilosort4 environment before running this scrip
# this code is adapted from the kilosort4 repository

# accept input argument as neuralRawDataPath
if __name__ == "__main__":
    import argparse
    import numpy as np
    import re
    import sys
    from pathlib import Path
    from kilosort.utils import download_probes
    from kilosort import run_kilosort, DEFAULT_SETTINGS
    download_probes() 
    
    settings = DEFAULT_SETTINGS.copy()
    
    # print(f"Results will be saved to: {settings['results_dir']}")
    settings['n_chan_bin'] = 64  # neuropixels1 has 385 channels (384+1). replace with your number of channels
    settings['batch_size'] = 64 * 1024  # default is 60000
    settings['nt'] = 41 # default. number of time samples per spike waveform. 61 as default for 30khz sampling rate, 41 for 20khz and 81 for 40khz
    settings['nblocks'] = 1 # default. recommended for single shank neuropixels probe.
    settings['tmin'] = 0  # default. start time in seconds
    settings['tmax'] = np.inf  # default. end time in seconds. by default all data is used but you can set a specific time range if there's noise in the beginning or end of your recording
    settings['probe_path'] = r"C:\\Users\scushing6\Desktop\TempKilosort\CN_doublesided_P2D_KSchmap.mat"
    settings['Th_universal'] = 8 # default is 9
    settings['Th_learned'] = 7 # default is 8,  reducing to handle cells appearing/disappearing
    
    parser = argparse.ArgumentParser(description='Run Kilosort4 on Neuropixels data.')
    parser.add_argument('--data_dir', type=str, required=True, help='Path to the directory containing the neural data.')
    parser.add_argument('--sess_str', type=str, default='', help='Sesssion string, e.g. NeuroPixels_X43_250719')
    parser.add_argument('--combined_ap_file', type=str, default='', help='Name of the combined ap.bin file if available (e.g. all_ap.bin)')
    args = parser.parse_args()
    # resolve the path to the data directory
    
    neuralRawDataPath = args.data_dir
    base_dir = Path(neuralRawDataPath).resolve()
    print(f"Using data directory: {base_dir}")
    settings['data_dir'] = neuralRawDataPath  # replace with your data path
    
    
    sess_str = args.sess_str
    combined_ap_file = args.combined_ap_file
    if sess_str == '' and combined_ap_file == '':
        # throw error
        raise ValueError("Either session string or combined .ap.bin path must be provided.")
    elif combined_ap_file != '' and sess_str != '':
        raise ValueError("Provide either session string or combined .ap.bin path, not both.")
    elif sess_str != '':
        print(f"Using session string: {sess_str}")
        # Pattern to match session folders
        session_pattern = re.compile(rf"{re.escape(sess_str)}_(\d+)_g\d+")

        # Dictionary to hold session number → .ap.bin Path
        ap_bin_files = {}

        for folder in base_dir.glob(f"{sess_str}_*_g*/"):
            print(f"Checking folder: {folder}")  # Debug: show which folders are found
            match = session_pattern.search(folder.name)
            if not match:
                print("No regex match")
                continue
            session_num = int(match.group(1))
            print(f"Found session number: {session_num}")  # Debug: show matched session number
            
            imec0_folder = folder / f"{folder.name}_imec0"
            if not imec0_folder.exists():
                print(f"imec0 folder does not exist: {imec0_folder}")
                continue
                
            # Find .ap.bin file in imec0 folder
            for ap_file in imec0_folder.glob("*.ap.bin"):
                ap_bin_files[session_num] = ap_file
                break  # Assume only one .ap.bin file per session

        # Get ordered list of Path objects by session
        ordered_ap_paths = [ap_bin_files[key] for key in sorted(ap_bin_files)]

        for path in ordered_ap_paths:
            print(path)
        # commented out the next line because it is easier to visualize raw waveform in phy2
        # settings['results_dir'] = base_dir / 'kilosort4'
        ops, st, clu, tF, Wall, similar_templates, is_ref, est_contam_rate, kept_spikes = \
        run_kilosort(settings, filename=ordered_ap_paths) 
    elif combined_ap_file != '':
        print(f"Using combined .ap.bin file: {combined_ap_file}")
        # commented out the next line because it is easier to visualize raw waveform in phy2. Not affected if combined ap.bin file is available!
        settings['results_dir'] = base_dir / 'kilosort4'
        # use base_dir and combined_ap_file to get combined ap.bin file path
        combined_ap_path = base_dir / combined_ap_file
        ops, st, clu, tF, Wall, similar_templates, is_ref, est_contam_rate, kept_spikes = \
        run_kilosort(settings, data_dir=base_dir, filename=combined_ap_path)

    # plot outpus of kilosort

    import numpy as np
    import pandas as pd
    from pathlib import Path

    # outputs saved to results_dir
    results_dir = Path(settings['data_dir']).joinpath('kilosort4_plots')
    # create results directory if it doesn't exist
    results_dir.mkdir(parents=True, exist_ok=True)
    ops = np.load(results_dir / 'ops.npy', allow_pickle=True).item()
    camps = pd.read_csv(results_dir / 'cluster_Amplitude.tsv', sep='\t')['Amplitude'].values
    contam_pct = pd.read_csv(results_dir / 'cluster_ContamPct.tsv', sep='\t')['ContamPct'].values
    chan_map =  np.load(results_dir / 'channel_map.npy')
    templates =  np.load(results_dir / 'templates.npy')
    chan_best = (templates**2).sum(axis=1).argmax(axis=-1)
    chan_best = chan_map[chan_best]
    amplitudes = np.load(results_dir / 'amplitudes.npy')
    st = np.load(results_dir / 'spike_times.npy')
    clu = np.load(results_dir / 'spike_clusters.npy')
    firing_rates = np.unique(clu, return_counts=True)[1] * 30000 / st.max()
    dshift = ops['dshift']


    import matplotlib.pyplot as plt
    from matplotlib import gridspec, rcParams
    rcParams['axes.spines.top'] = False
    rcParams['axes.spines.right'] = False
    gray = .5 * np.ones(3)

    fig = plt.figure(figsize=(10,10), dpi=100)
    grid = gridspec.GridSpec(3, 3, figure=fig, hspace=0.5, wspace=0.5)

    ax = fig.add_subplot(grid[0,0])
    ax.plot(np.arange(0, ops['Nbatches'])*2, dshift);
    ax.set_xlabel('time (sec.)')
    ax.set_ylabel('drift (um)')

    ax = fig.add_subplot(grid[0,1:])
    t0 = 0
    t1 = np.nonzero(st > ops['fs']*5)[0][0]
    ax.scatter(st[t0:t1]/30000., chan_best[clu[t0:t1]], s=0.5, color='k', alpha=0.25)
    ax.set_xlim([0, 5])
    ax.set_ylim([chan_map.max(), 0])
    ax.set_xlabel('time (sec.)')
    ax.set_ylabel('channel')
    ax.set_title('spikes from units')

    ax = fig.add_subplot(grid[1,0])
    nb=ax.hist(firing_rates, 20, color=gray)
    ax.set_xlabel('firing rate (Hz)')
    ax.set_ylabel('# of units')

    ax = fig.add_subplot(grid[1,1])
    nb=ax.hist(camps, 20, color=gray)
    ax.set_xlabel('amplitude')
    ax.set_ylabel('# of units')

    ax = fig.add_subplot(grid[1,2])
    nb=ax.hist(np.minimum(100, contam_pct), np.arange(0,105,5), color=gray)
    ax.plot([10, 10], [0, nb[0].max()], 'k--')
    ax.set_xlabel('% contamination')
    ax.set_ylabel('# of units')
    ax.set_title('< 10% = good units')

    for k in range(2):
        ax = fig.add_subplot(grid[2,k])
        is_ref = contam_pct<10.
        ax.scatter(firing_rates[~is_ref], camps[~is_ref], s=3, color='r', label='mua', alpha=0.25)
        ax.scatter(firing_rates[is_ref], camps[is_ref], s=3, color='b', label='good', alpha=0.25)
        ax.set_ylabel('amplitude (a.u.)')
        ax.set_xlabel('firing rate (Hz)')
        ax.legend()
        if k==1:
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_title('loglog')
    # save the figure
    fig.savefig(results_dir / 'kilosort4_stats.png', dpi=300, bbox_inches='tight')

    probe = ops['probe']
    # x and y position of probe sites
    xc, yc = probe['xc'], probe['yc']
    nc = 16 # number of channels to show
    good_units = np.nonzero(contam_pct <= 0.1)[0]
    mua_units = np.nonzero(contam_pct > 0.1)[0]


    gstr = ['good', 'mua']
    for j in range(2):
        print(f'~~~~~~~~~~~~~~ {gstr[j]} units ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        print('title = number of spikes from each unit')
        units = good_units if j==0 else mua_units
        fig = plt.figure(figsize=(12,3), dpi=150)
        grid = gridspec.GridSpec(2,20, figure=fig, hspace=0.25, wspace=0.5)

        for k in range(40):
            wi = units[np.random.randint(len(units))]
            wv = templates[wi].copy()
            cb = chan_best[wi]
            nsp = (clu==wi).sum()

            ax = fig.add_subplot(grid[k//20, k%20])
            n_chan = wv.shape[-1]
            ic0 = max(0, cb-nc//2)
            ic1 = min(n_chan, cb+nc//2)
            wv = wv[:, ic0:ic1]
            x0, y0 = xc[ic0:ic1], yc[ic0:ic1]

            amp = 4
            for ii, (xi,yi) in enumerate(zip(x0,y0)):
                t = np.arange(-wv.shape[0]//2,wv.shape[0]//2,1,'float32')
                t /= wv.shape[0] / 20
                ax.plot(xi + t, yi + wv[:,ii]*amp, lw=0.5, color='k')

            ax.set_title(f'{nsp}', fontsize='small')
            ax.axis('off')
        plt.show()
        # save the figure for good and mua units
        fig.savefig(results_dir / f'kilosort4_{gstr[j]}_units.png', dpi=300, bbox_inches='tight')
