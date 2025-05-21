import importlib
import numpy as np
import tmd

import os
from typing import List, Dict


import gtda
from gtda.diagrams import Silhouette

def check_and_load_neurons(
    main_path: str,
    folders_to_treat: List[str]
) -> Dict[str, List]:
    grouped_neurons = {}

    for folder in folders_to_treat:
        data_path = os.path.join(main_path, folder)
        neurons = []
        loaded_files = []
        all_asc_files = []
        log_lines = []

        if not os.path.exists(data_path):
            msg = f"Folder '{data_path}' does not exist.\n"
            print(msg)
            log_lines.append(msg)
        else:
            for filename in os.listdir(data_path):
                if filename.endswith(".asc"):
                    all_asc_files.append(filename)
                    filepath = os.path.join(data_path, filename)
                    try:
                        neu = tmd.io.load_neuron_from_morphio(filepath)
                        neurons.append(neu)
                        loaded_files.append(filename)
                    except Exception as e:
                        log_lines.append(f"Failed to load {filename}: {e}\n")

            unloaded_files = set(all_asc_files) - set(loaded_files)
            total_file_count = sum(
                os.path.isfile(os.path.join(data_path, f)) for f in os.listdir(data_path)
            )

            log_lines.append(f"Folder: {folder}")
            log_lines.append(f"Neurons loaded: {len(neurons)}")
            log_lines.append(f"Total files in folder: {total_file_count}")
            log_lines.append(f".asc files in folder: {len(all_asc_files)}")
            log_lines.append(f"Failed to load: {len(unloaded_files)} files")
            log_lines.append(f"Unloaded files: {sorted(unloaded_files)}")

        summary = "\n".join(log_lines) + "\n" + "=" * 40 + "\n"

        log_file_path = os.path.join(main_path, f"log_{folder}.txt")
        with open(log_file_path, "w") as f:
            f.write(summary)

        print(summary)

        # Store the loaded neurons for this folder
        grouped_neurons[folder] = neurons

    return grouped_neurons

def compute_apical_diagrams(
    neurons: List,
    group_name: str,
    neuron_part: str
) -> List:
    """
    Compute the persistence diagrams for a specified list of neurons

    Arguments:
        - neurons: list of the neurons we want to compute the persistence diagram of
        - group_name: name of sample of neurons
        - neuron_part: specify the type of neurite we are interested in. 
            Possible arguments: apical_dendrite, basal_dendrite, axon

    Return:
        - diagrams: list of the computed persistence diagrams
    """
    
    diagrams = []
    
    for idx, neu in enumerate(neurons):
        try:
            ph = tmd.methods.get_ph_neuron(neu, neurite_type=neuron_part)
            if ph is not None:
                diagrams.append(ph)
            else:
                print(f"[{group_name}] Neuron {idx} has no {neuron_part}.")
        except Exception as e:
            print(f"[{group_name}] Failed to compute {neuron_part} diagram for neuron {idx}: {e}")
    
    print(f"[{group_name}] Computed {len(diagrams)} persistence diagrams for {neuron_part}s.")
    
    return diagrams

def format_diagrams_for_gtda(diagrams, hom_dim=1):
    """
    Convert diagrams from (death, birth) or (birth, death) to (birth, death, dim),
    flipping coordinates only if birth > death.
    """
    formatted = []
    for diag in diagrams:
        arr = np.array(diag)
        corrected = []
        for point in arr:
            birth, death = sorted(point[:2])  # ensures birth ≤ death
            corrected.append([birth, death])
        corrected = np.array(corrected)
        dim_col = np.full((corrected.shape[0], 1), hom_dim)
        diagram_hom = np.hstack([corrected, dim_col])  # (n_points, 3)
        formatted.append(diagram_hom)
    return formatted


def pad_diagrams(diagrams, dim=1):
    max_len = max(len(d) for d in diagrams)
    padded = []
    for d in diagrams:
        n = len(d)
        if n < max_len:
            # use zeros instead of -inf
            padding = np.zeros((max_len - n, 3))
            padding[:, 2] = dim  # preserve homology dimension
            d = np.vstack([d, padding])
        padded.append(d)
    return np.stack(padded)


def bootstrap_sample_from_diagrams(
    diagrams: np.ndarray,
    B: int
):
    """
    Apply the bootstrap procedure using persistence diagrams to compute
    distances between bootstrap mean silhouettes and the original mean.

    Arguments:
        - diagrams : np.ndarray
            Original persistence diagrams
        - B : int
            Number of bootstrap repetitions

    Returns:
        - resampled_silhouettes_mouse : 
            Silhouette persistence from the resampled diagrams
        - bootstrap_distances : List[float]
            Sup-norm distances between bootstrap means and original mean.
        - bootstrap_means : List[np.ndarray]
            List of (1, n_bins) mean silhouette curves from each bootstrap.
        - bootstrap_diagrams : List[np.ndarray]
            List of (n, ?, 3) arrays of the diagrams used in each bootstrap round.
    """
    n = len(diagrams)
    bootstrap_distances = []
    bootstrap_means = []
    bootstrap_diagrams = []
    resampled_silhouettes = []
    sil = Silhouette()

    #Compute the original mean
    formated_og = format_diagrams_for_gtda(diagrams)
    padded_og= pad_diagrams(formated_og)
    silhouettes = sil.fit_transform(padded_og)
    silhouettes_mean = np.mean(silhouettes, axis=0)


    #Bootstrap procedure
    for _ in range(B):
        # Step 1: Resample diagrams with replacement
        indices = np.random.choice(n, size=n, replace=True)
        resampled_diagrams = [diagrams[i] for i in indices]

        # Step 2: Transform the resampled diagrams to fit the right format for the silhouette
        formated_diagrams = format_diagrams_for_gtda(resampled_diagrams)
        padded_diagrams = pad_diagrams(formated_diagrams)

        # Step 3: Compute silhouette summaries for resampled diagrams
        resampled_silhouette = sil.fit_transform(padded_diagrams)

        # Step 4: Compute bootstrap mean silhouette
        bootstrap_mean = np.mean(resampled_silhouette, axis=0)  
        
        # Step 5: Compute L∞ (sup) distance to original mean
        dist = np.max(np.abs(bootstrap_mean - silhouettes_mean))

        # Store results
        resampled_silhouettes.append(resampled_silhouette)
        bootstrap_distances.append(dist)
        bootstrap_means.append(bootstrap_mean)
        bootstrap_diagrams.append(resampled_diagrams)

    return silhouettes_mean, resampled_silhouettes, bootstrap_distances, bootstrap_means, bootstrap_diagrams