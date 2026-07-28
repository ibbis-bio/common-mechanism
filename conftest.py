"""
Additional Configurations for Pytest
"""


def pytest_addoption(parser):
    """Adds unique argument to pytest for commec database example outputs generation."""
    print("Test Configuration loaded!")
    parser.addoption(
        "--gen-examples",
        action="store_true",
        default=False,
        help="Generate exemplar output files instead of testing against them.",
    )
    parser.addoption(
        "--commec-db-dir", action="store", default=None,
        help=(
            "Path to a real commec reference database directory. Required for the"
            " end-to-end detection tests in test_evaluation.py; those tests are skipped"
            " when it is not provided. May also be set via the COMMEC_DB_DIR env var."
        )
    )
