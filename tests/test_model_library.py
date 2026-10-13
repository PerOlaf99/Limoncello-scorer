"""A fitted standard survives the disk.

Auto-genotyping fits a plate's internal-standard shape from the operator's
marks and then propagates it across the plate.  That fit is written to the
library so a later plate of the same assay can be seeded with it, and so a
restart never costs it.  These tests pin the two halves of that promise: the
plain-data round-trip that makes a model storable, and the store that keeps it.
"""
import json
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import genotyping  # noqa: E402
import mark_library  # noqa: E402


def _model():
    """A four-band model with every primed reference set filled in."""
    base = 2600.0
    model = genotyping.PlateISModel.from_marks(
        [[base, base + 81, base + 145, base + 155] for _ in range(4)])
    model.std_snr = [12.0, 20.0, 15.0, 18.0]
    model.std_heights = [900.0, 950.0, 880.0, 910.0]
    model.sample_snr = [40.0, 41.0, 39.0, 42.0]
    model.co_mig = [[-10.0, -9.0, -11.0], [], [3.0, 4.0], []]
    return model


def test_geometry_roundtrips_including_merged():
    merged = genotyping.ISGeometry(81, 64, None, center_tol=200.0)
    back = genotyping.ISGeometry.from_dict(merged.to_dict())
    assert back.d1 == merged.d1 and back.d2 == merged.d2
    assert back.d3 is None and back.merged is True
    assert back.tol1 == merged.tol1 and back.center_tol == merged.center_tol

    split = genotyping.ISGeometry(81, 64, 10)
    back2 = genotyping.ISGeometry.from_dict(split.to_dict())
    assert back2.d3 == 10 and back2.merged is False
    assert back2.tol3 == split.tol3


def test_model_roundtrips_every_field():
    model = _model()
    back = genotyping.PlateISModel.from_dict(model.to_dict())
    assert back is not None
    assert back.geometry.to_dict() == model.geometry.to_dict()
    assert back.merged_geometry is None
    assert back.centers == model.centers
    assert back.n_marked == model.n_marked
    assert back.merged_seen == model.merged_seen
    assert back.std_snr == model.std_snr
    assert back.std_heights == model.std_heights
    assert back.sample_snr == model.sample_snr
    assert back.co_mig == model.co_mig
    assert back.scaled == model.scaled


def test_roundtripped_model_still_finds_the_standard():
    model = _model()
    rng = np.random.default_rng(3)
    trace = rng.normal(0.0, 1.0, 3400)
    for s in (2600, 2681, 2745, 2755):
        for i in range(s - 5, s + 6):
            trace[i] += 800.0 * np.exp(-0.5 * ((i - s) / 2.0) ** 2)
    found = model.find(trace)
    assert found is not None and list(found[0]) == [2600, 2681, 2745, 2755]
    # The loaded model must read the plate exactly as the fitted one did.
    back = genotyping.PlateISModel.from_dict(model.to_dict())
    assert back.find(trace) == found


def test_model_dict_is_json_safe():
    json.dumps(_model().to_dict())        # must not raise on numpy scalars


def test_merged_model_roundtrips():
    base = 2600.0
    geom = genotyping.ISGeometry(81, 64, None)
    model = genotyping.PlateISModel(geom, [base, base + 5], n_marked=2,
                                    merged_seen=True, scaled=1.6)
    back = genotyping.PlateISModel.from_dict(model.to_dict())
    assert back.geometry.merged is True
    assert back.merged_seen is True
    assert back.scaled == 1.6


def test_from_dict_refuses_a_record_it_cannot_place_with():
    assert genotyping.PlateISModel.from_dict(None) is None
    assert genotyping.PlateISModel.from_dict({}) is None
    # Geometry but no centre: nothing to anchor a search to.
    assert genotyping.PlateISModel.from_dict(
        {"geometry": {"d1": 81, "d2": 64}}) is None
    assert genotyping.ISGeometry.from_dict("nope") is None
    assert genotyping.ISGeometry.from_dict({}) is None


def test_library_stores_and_lists_a_model(tmp_path, monkeypatch):
    monkeypatch.setattr(mark_library, "LIBRARY_DIR", tmp_path)

    assert mark_library.load_model("runA") == {}
    assert mark_library.saved_models() == []

    model = _model()
    mark_library.save_model("runA", model.to_dict())

    payload = mark_library.load_model("runA")
    assert payload["centers"] == model.centers
    stored = mark_library.saved_models()
    assert [run for run, _updated, _data in stored] == ["runA"]
    assert genotyping.PlateISModel.from_dict(stored[0][2]) is not None


def test_model_files_do_not_look_like_mark_runs(tmp_path, monkeypatch):
    monkeypatch.setattr(mark_library, "LIBRARY_DIR", tmp_path)
    mark_library.save_model("runA", _model().to_dict())
    # The run scan reads library/*.json, not the models/ subdirectory.
    assert mark_library.runs() == []
    assert mark_library.summary() == []


def test_save_model_overwrites_and_is_json_safe(tmp_path, monkeypatch):
    monkeypatch.setattr(mark_library, "LIBRARY_DIR", tmp_path)
    mark_library.save_model("runA", _model().to_dict())
    other = _model()
    other.centers = [1.0, 2.0]
    mark_library.save_model("runA", other.to_dict())
    assert mark_library.load_model("runA")["centers"] == [1.0, 2.0]
    assert len(mark_library.saved_models()) == 1
    stored = (tmp_path / "models" / "runA.json").read_text(encoding="utf-8")
    json.loads(stored)
