import logging
from importlib.resources import files
from pathlib import Path

import pandas as pd

from src.model.ensemble_model import EnsembleWrapper
from src.utils.data_utils import prepare_inference_data
from src.utils.filehandling import prepare_location
from src.utils.log import logger


def load_features_list(model_dir) -> list[str]:
    """Load the list of features to be used for inference.
    TODO: Not sure yet if this is the correct place to put this function.
    """
    lsf_remaining_feats = model_dir / "features_min_sds_above_25.csv"
    if not lsf_remaining_feats.exists():  # type: ignore
        raise FileNotFoundError(
            f"Features list file not found at {lsf_remaining_feats}"
        )
    lsf_remaining_feats = pd.read_csv(lsf_remaining_feats)  # type: ignore
    lsf_remaining_feats = lsf_remaining_feats.iloc[:, 0].tolist()
    return lsf_remaining_feats


def main(args):
    if args.debug:
        logger.setLevel(logging.DEBUG)

    input_file = Path(args.input).resolve()
    output_dir = Path(args.output).resolve()
    # Load and prepare data
    lsf_remaining_feats = load_features_list(files("data") / "features" / args.model)
    X = prepare_inference_data(input_file, lsf_remaining_feats)
    # Load ensemble model
    ensemble_model = EnsembleWrapper(
        models_dir=files("data") / "models" / args.model,  # type: ignore
        n_threads=args.threads,
    )
    # Make predictions
    probas, preds = ensemble_model.predict(X)
    # Save results
    prepare_location(output_dir / "predicted_probabilities.csv", args.create_dir)
    probas.T.sort_values(by="mod", ascending=False).to_csv(
        output_dir / "predicted_probabilities.csv"
    )
    preds.to_json(output_dir / "predicted_labels.json")

    logger.info(f"Saved results to {output_dir}")
