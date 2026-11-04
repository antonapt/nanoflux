from pathlib import Path
from pickle import load

import numpy as np
from sklearn.linear_model import LogisticRegression


def load_model(path: Path | str) -> LogisticRegression:
    path = Path(path)

    with open(path.with_suffix(".features"), "rb") as file:
        feature_names = load(file)

    with open(path.with_suffix(".onehot"), "rb") as file:
        onehot = load(file)

    model_path = path.with_suffix(".model")
    header = {}
    with model_path.open() as file:
        for line in file:
            line = line.rstrip()
            if line == "w":
                break
            key, *values = line.split()
            header[key] = values

    labels = np.asarray(header["label"], dtype=np.int64)
    n_features = int(header["nr_feature"][0])
    n_classes = int(header["nr_class"][0])
    bias = float(header["bias"][0])

    if bias != -1:
        raise ValueError("Models with bias are not supported by this loader yet.")
    if labels.size != n_classes or not np.array_equal(
        np.sort(labels), np.arange(1, n_classes + 1)
    ):
        raise ValueError("Invalid liblinear label metadata.")

    weights = np.loadtxt(model_path, skiprows=6).reshape(n_features, n_classes)
    order = np.argsort(labels)

    model = LogisticRegression(solver="liblinear", fit_intercept=False)
    model.classes_ = np.asarray(onehot.categories_[0])
    model.coef_ = weights[:, order].T
    model.intercept_ = np.zeros(n_classes, dtype=np.float64)
    model.n_features_in_ = n_features
    model.feature_names_in_ = np.asarray(feature_names, dtype=object)
    return model
