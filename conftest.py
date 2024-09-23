"""
Additional Configurations for Pytest
"""
def pytest_addoption(parser):
    """ Adds unique argument to pytest for commec database example outputs generation."""
    print("Test Configuration loaded!")
    parser.addoption(
        "--gen-examples", action="store_true", default=False,
        help="Generate exemplar output files instead of testing against them."
    )
