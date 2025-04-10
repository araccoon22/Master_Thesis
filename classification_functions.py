import importlib
import numpy as np
import tmd
import os
import sklearn.tree
import sklearn 

from typing import List, Dict


def load_neurons_from_folders(
    data_path:str,
    folders_to_treat: List
):
    all_neurons = {}

    for folder in folders_to_treat:
        folder_path = os.path.join(data_path, folder)
        neurons = []

        for filename in os.listdir(folder_path):
            if filename.endswith(".asc"):
                filepath = os.path.join(folder_path, filename)
                neu = tmd.io.load_neuron_from_morphio(filepath)
                neurons.append(neu)

        all_neurons[folder] = neurons
        print(f"Neurons loaded from {folder}: {len(neurons)}")

    return all_neurons

def train(
    mod, 
    classifier, 
    data, 
    labels, 
    **kwargs
):
    clas_mod = importlib.import_module("sklearn." + mod)
    clf = getattr(clas_mod, classifier)()

    # Apply regularization if needed
    if classifier == "LinearDiscriminantAnalysis":
        clf.set_params(solver="lsqr", shrinkage="auto")
    elif classifier == "QuadraticDiscriminantAnalysis":
        clf.set_params(reg_param=0.1)

    clf.set_params(**kwargs)
    clf.fit(data, labels)
    return clf


def predict(
    clf, 
    data
):
    """Predict label for data for the trained classifier clf.

    Returns the index of the predicted class for each datapoint in data.
    """
    predict_label = clf.predict([data])

    return predict_label[0]


def classify_cell_in_groups(
    all_neurons: Dict[str, List],
    types: List[str],
    cell_to_classify: str,
    neurite_type: str,
    number_of_trials: int,
    classifier_combinations: List[tuple],  # (module, classifier)
):
    """Classify the cell using multiple classifiers and return prediction percentages."""

    # ------------------------ Training dataset --------------------------------
    groups = [all_neurons[group] for group in types]
    labels = [i + 1 for i, group in enumerate(groups) for _ in group]

    pers_diagrams = [
        tmd.methods.get_ph_neuron(neuron, neurite_type=neurite_type)
        for group in groups
        for neuron in group
    ]
    xlim, ylim = tmd.vectorizations.get_limits(pers_diagrams)
    pers_images = [
        tmd.analysis.persistence_image_data(p, xlim=xlim, ylim=ylim)
        for p in pers_diagrams
    ]
    train_dataset = [img.flatten() for img in pers_images]

    # ------------------------ Test cell dataset -------------------------------
    neuron2test = tmd.io.load_neuron_from_morphio(cell_to_classify)
    pers2test = tmd.methods.get_ph_neuron(neuron2test, neurite_type=neurite_type)
    pers_image2test = tmd.analysis.persistence_image_data(pers2test, xlim=xlim, ylim=ylim)
    test_dataset = pers_image2test.flatten()

    # ------------------------ Try multiple classifiers ------------------------
    results = {}

    for mod, method in classifier_combinations:
        predict_labels = []

        for _ in range(number_of_trials):
            clf = train(mod, method, train_dataset, labels)
            predict_labels.append(predict(clf, test_dataset))

        percentages = {
            types[i - 1]: float(np.sum(np.array(predict_labels) == i)) / len(predict_labels)
            for i in np.unique(labels)
        }

        key = f"{mod}.{method}"
        results[key] = percentages

    return results