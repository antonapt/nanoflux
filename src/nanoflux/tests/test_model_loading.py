from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from nanoflux.model.linear_model import LiblinearWrapper, load_model

MODELS_DIR = Path(__file__).parents[3] / "data" / "models"
MODEL_PATHS = sorted(MODELS_DIR.glob("*.model"))


@pytest.mark.parametrize("model_path", MODEL_PATHS, ids=lambda path: path.stem)
def test_loaders_produce_same_predictions(model_path: Path) -> None:
    sklearn_model = load_model(model_path)

    wrapper = LiblinearWrapper()
    wrapper.load_model(model_path)

    rng = np.random.default_rng(0)
    features = sklearn_model.feature_names_in_
    X = pd.DataFrame(
        rng.integers(0, 2, size=(32, len(features))),
        columns=features,
    )

    np.testing.assert_allclose(
        sklearn_model.predict_proba(X),
        wrapper.predict_proba(X),
        rtol=1e-5,
        atol=1e-7,
    )
    breakpoint()
    np.testing.assert_array_equal(sklearn_model.predict(X), wrapper.predict(X))
