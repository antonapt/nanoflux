from pathlib import Path
from pickle import load

import liblinear.liblinearutil as ll
import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator


class LiblinearWrapper(BaseEstimator):
    def __init__(self):
        self.model = None
        self.onehot = None
        self.is_fitted = False
        self.feature_names_in_: list | None = None

    def encode_y(self, y):
        if self.onehot is None:
            raise ValueError(
                "OneHotEncoder is not initialized, please make sure for each `.model` there exists a `.onehot` file."
            )
        y_enc = self.onehot.fit_transform(y.to_numpy().reshape(-1, 1))
        y_enc = np.argmax(y_enc, axis=1) + 1  # liblinear needs labels starting at 1
        return y_enc

    def _to_numpy(self, X: pd.DataFrame):
        return np.ascontiguousarray(X.to_numpy(dtype=np.float32))

    def _prepare_X(self, X: pd.DataFrame):
        if len(X.columns) != len(self.feature_names_in_) or not all(
            X.columns == self.feature_names_in_
        ):
            X = X.reindex(
                columns=self.feature_names_in_,
            )
        return X

    def predict(self, X: pd.DataFrame):
        if (
            self.onehot is not None
            and self.is_fitted
            and self.feature_names_in_ is not None
        ):
            X = self._prepare_X(X)

            X_np = self._to_numpy(X)
            labels, _, _ = ll.predict([], X_np, self.model, "-q")
            labels = (
                np.array(labels).astype(int) - 1
            )  # convert liblinear output to zero-based
            decoded = self.onehot.categories_[0][labels]
            return decoded
        raise ValueError(
            "Model was not correctly loaded. Please make sure to have a .model, .features, and .onehot file."
        )

    def _predict_full(self, X: pd.DataFrame, cmd: str = "-b 1 -q"):
        if (
            self.onehot is not None
            and self.is_fitted
            and self.feature_names_in_ is not None
        ):
            X = self._prepare_X(X)

            X_np = self._to_numpy(X)
            _, _, probs = ll.predict([], X_np, self.model, cmd)

            return np.array(probs)[:, np.argsort(self.model.get_labels())]
        raise ValueError(
            "Model was not correctly loaded. Please make sure to have a .model, .features, and .onehot file."
        )

    def predict_proba(self, X: pd.DataFrame):
        return self._predict_full(X, cmd="-b 1 -q")

    def predict_logits(self, X: pd.DataFrame):
        return self._predict_full(X, cmd="-b 0 -q")

    def load_model(self, path: Path | str):
        path = Path(path)

        if self.model is None:
            self.model = ll.load_model(str(path.with_suffix(".model")))

            with open(path.with_suffix(".features"), "rb") as pkl:
                self.feature_names_in_ = load(pkl)

            with open(path.with_suffix(".onehot"), "rb") as pkl:
                self.onehot = load(pkl)

            self.is_fitted = True
        else:
            raise ValueError(
                "Model does already exist."
            )  # TODO: not sure if we should raise this error or just reload the model
