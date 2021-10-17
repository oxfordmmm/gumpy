# conftest.py

import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--run_slow", action="store_true", default=False, help="run slow tests"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as involving M. tuberculosis and therefore slow to run")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--run_slow"):
        # --runslow given in cli: do not skip slow tests
        return
    skip_slow = pytest.mark.skip(reason="need --run_slow option to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)
