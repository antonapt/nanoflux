from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.linear_model import LogisticRegression

from nanoflux.model.linear_model import load_model
from nanoflux.utils.download import download_models
from nanoflux.utils.log import logger


class EnsembleModel:
    def __init__(self, models_dir: Path | str, n_threads: int = -1):
        self.models_dir = Path(models_dir)
        self.n_threads = n_threads
        self.models: list[LogisticRegression] = []
        self.classes_: np.ndarray | None = None
        if not self.models_dir.exists():
            raise FileNotFoundError(
                f"Models directory {self.models_dir} does not exist."
            )
        self._load_ensemble_model()

    def _load_ensemble_model(self):
        """
        Finds all .model files in the directory and loads
        the corresponding LiblinearWrapper instances.
        """
        logger.info(f"Loading ensemble models from {self.models_dir}...")
        model_files = list(self.models_dir.glob("*.model"))

        if not model_files:
            try:
                download_models()
            except Exception:
                raise FileNotFoundError(f"No .model files found in {self.models_dir}")

        logger.info(f"Loading {len(model_files)} models...")

        self.models = Parallel(n_jobs=self.n_threads, backend="loky")(
            delayed(load_model)(p) for p in model_files
        )

        self._models_check()
        self.classes_ = self.models[0].classes_
        logger.info("Ensemble model loaded successfully.")

    def _models_check(self):
        """
        Ensures all models have the same categories.
        """
        if not self.models:
            raise ValueError("No models loaded in the ensemble.")

        reference_categories = self.models[0].classes_

        for i, model in enumerate(self.models[1:], start=1):
            if not np.array_equal(reference_categories, model.classes_):
                raise ValueError(f"Categories mismatch between model 0 and model {i}.")

    def _run_single_prediction(self, model, X, logits):
        """Helper function that lives outside the class or as a static method"""
        X = X.reindex(model.feature_names_in_, axis=1)
        if logits:
            return model.predict_logits(X)
        return model.predict_proba(X)

    def _predict_full(self, X: pd.DataFrame, logits: bool = False) -> np.ndarray:
        """
        Predicts the average probabilities or logits across all models in the ensemble.
        """
        results = Parallel(n_jobs=self.n_threads, backend="loky")(
            delayed(self._run_single_prediction)(m, X, logits) for m in self.models
        )
        return np.mean(results, axis=0)

    def predict_logits(self, X: pd.DataFrame) -> np.ndarray:
        """Predicts the average logits across all models in the ensemble."""
        return self._predict_full(X, logits=True)

    def predict_proba(self, X: pd.DataFrame) -> np.ndarray:
        """Predicts the average probabilities across all models in the ensemble."""
        return self._predict_full(X, logits=False)

    def predict(self, X: pd.DataFrame) -> np.ndarray:
        """
        Predicts the class with the highest average probability.
        """
        avg_probas = self.predict_proba(X)
        avg_probas = pd.DataFrame(avg_probas, index=X.index, columns=self.classes_)

        avg_preds = pd.DataFrame(
            dict(preds=avg_probas.idxmax(axis=1), confidence=avg_probas.max(axis=1)),
            index=X.index,
        )
        return avg_probas, avg_preds
