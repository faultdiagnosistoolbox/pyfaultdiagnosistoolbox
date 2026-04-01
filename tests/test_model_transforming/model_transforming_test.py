import time
from unittest.mock import patch
import faultdiagnosistoolbox as fdt
import faultdiagnosistoolbox.ModelTransforming as mt
from test_models import *


def test_subexpression_identification():
    """Test that subexpressions in model input are correctly identified"""
    analyzed_model, reduced = mt.analyze_modelDef(subexpr_id_modelDef)
    assert analyzed_model == subexpr_id_correct_reduction


def test_logging_timing():
    """Test the logging timing. Requirement #13"""
    events = []

    def timed_print(*args, **kwargs):
        events.append((time.monotonic(), args))

    start_time = time.monotonic()
    # Temporarly replace built in print function with modified one
    with patch("builtins.print", timed_print):
        mt.get_optimized_modelDef(bad_subexpr_modelDef)

    end_time = time.monotonic()

    assert len(events) > 0

    # Get list of timestamps from events
    timestamps = [x[0] for x in events]

    # Check if process ran for longer than 5s
    if (end_time - start_time) > 5:
        # Check that the distance in time for every timestep is under 1s
        for i, timestamp in enumerate(timestamps[1:], start=1):
            assert timestamp - timestamps[i - 1] < 1


def test_identical_result():
    """Test that same input multiple times produces identical result"""
    results = []
    for i in range(10):
        result = mt.get_optimized_modelDef(bad_subexpr_modelDef)
        if result not in results:
            results.append(result)
    assert len(results) == 1


def test_timed_model_transformation():
    """Test that the model transformation for a model with 10 equations takes less than 20s"""
    start_time = time.monotonic()
    mt.get_optimized_modelDef(twelve_eq_modelDef)
    end_time = time.monotonic()
    assert (end_time - start_time) < 20


def test_find_and_replace():
    """Test that given a model, the identified reoccuring subexpressions and their replacements, a new model
    is returned with the replacements made. Requirement #10, #14"""
    reduced_modelDef = mt.find_and_replace(
        subexpr_id_modelDef,
        subexpr_id_correct_reduction,
    )

    assert reduced_modelDef == subexpr_id_reduced_modelDef
    result = 0
    offset = len(reduced_modelDef) - len(subexpr_id_modelDef["rels"])
    # Iterate backwards to get equivalent equations from both
    for i in range(len(subexpr_id_modelDef["rels"]) - 1, 0, -1):
        result += abs(
            reduced_modelDef["rels"][i + offset] - subexpr_id_modelDef["rels"][i]
        )
    assert result == 0


def test_view_reduced_subexprs(capsys):
    """Test that user can view manipulated expressions"""
    mt.view_reduced_subexprs(bad_subexpr_modelDef)
    captured = capsys.readouterr()
    replacements, reduced = mt.analyze_modelDef(bad_subexpr_modelDef)
    assert str(replacements) in captured.out


def test_model_copy():
    """Test that copy of model is properly created"""
    model_copy = mt.create_model_copy(twelve_eq_modelDef)
    assert twelve_eq_modelDef == model_copy
    assert id(twelve_eq_modelDef) != id(model_copy)


def test_model_transformation(capsys):
    """Test that model is transformed to representation that satisfies toolbox limitations"""
    test_modelDef = mt.get_optimized_modelDef(bad_subexpr_modelDef)
    test = fdt.DiagnosisModel(test_modelDef)
    test.Lint()
    captured = capsys.readouterr()
    assert "validation finished with 0 errors" in captured.out
