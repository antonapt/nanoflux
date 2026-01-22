from pathlib import Path

import numpy as np
import pandas as pd

from src.utils.log import logger


def load_pipeline_data(sample_path: Path) -> pd.DataFrame:
    logger.info("Reading sample.")

    df = pd.read_csv(sample_path, sep="\t", header=None)

    df = df.iloc[:, [0, 1, 2, 9, 13]]
    df.loc[:, 9] = df[9].apply(lambda x: x.split(" ")[1]).astype(float)
    df.columns = ["chr", "start", "end", "mod", "probe"]
    df.loc[:, "mod"] = df.loc[:, "mod"] * 0.01
    df.loc[df["mod"] >= 0.7, "mod"] = 1.0
    df.loc[df["mod"] < 0.3, "mod"] = -1.0
    df.loc[(df["mod"] >= 0.3) & (df["mod"] < 0.7), "mod"] = np.nan
    df = df.set_index("probe")
    df = df.loc[:, "mod"].to_frame()

    return df


def normalize_X(X: pd.DataFrame):
    """Normalize the input data matrix X by scaling each row based on the number of measured features."""
    measured_features = np.sum(X != 0, axis=1)
    available_features = X.shape[1]
    scaling_factor = available_features / measured_features
    X = X.mul(scaling_factor, axis=0)
    return X


def prepare_features(X: pd.DataFrame, curr_features_used: list) -> pd.DataFrame:
    """Prepare features by ensuring all required features are present, filling missing ones with 0."""
    if isinstance(X, pd.DataFrame):
        if X.shape[1] > 1:
            logger.warning(
                "Input data has more than one column; expected a single-column DataFrame. Using the first column only."
            )
        X = X.iloc[:, 0]
    elif not isinstance(X, pd.Series):
        raise ValueError("Input X must be a pandas DataFrame or Series.")

    if not np.all([val in [-1, 0, 1] or np.isnan(val) for val in sorted(X.unique())]):
        raise ValueError("Input data contains values other than -1, 0, 1, or NaN.")

    logger.info("Low-variance site filtering.")
    X_prepared = pd.Series(0.0, index=curr_features_used, name=X.name)
    common_rows = X.index.intersection(X_prepared.index)
    X_prepared.loc[common_rows] = X.loc[common_rows].astype(np.float32)
    X_prepared[X_prepared.isna()] = 0
    X_prepared = X_prepared.to_frame().T
    return X_prepared


def prepare_inference_data(input_file: Path, curr_features_used: list) -> pd.DataFrame:
    X = load_pipeline_data(input_file)
    logger.info("Loaded preprocessed data.")
    X_prepared = prepare_features(X, curr_features_used)
    logger.info("Prepared features for inference.")
    X_prepared = normalize_X(X_prepared)
    logger.info("Normalized features for inference.")
    return X_prepared
