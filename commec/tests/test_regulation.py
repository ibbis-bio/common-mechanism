import logging
import pytest

from commec.regulation.containers import (
    add_regulated_list,
    RegulationList
)
from commec.utils.logger import setup_console_logging


@pytest.mark.parametrize("new_list,expected_outcome", [
    pytest.param(*case) for case in [
        (RegulationList("list01","L1","www.l1.com",["NZ"],1), True), # Duplicate
        (RegulationList("list01","L2","www.l1.com",["NZ"],0), True), # Non-Duplicate
        (RegulationList("list01","L1","www.l2.com",["NZ"],0), False), # Wrong URL
        (RegulationList("list01","L1","www.l1.com",["NZ", "AU"],0), False), # Wrong Regions
        (RegulationList("list01","L1","www.l1.com",["AU"],0), False), # Wrong Region
        (RegulationList("list02","L1","www.l1.com",["NZ"],0), False), # Wrong Name
    ]
])
def test_list_overwrite_protections(new_list, expected_outcome):
    """
    Tests what occurs when two lists share an acronym that should not be concatenated.
    """
    setup_console_logging(logging.DEBUG)
    existing_list = RegulationList("list01","L1","www.l1.com",["NZ"],0)
    add_regulated_list(existing_list)
    assert expected_outcome == add_regulated_list(new_list)
