# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science

"""
Responsible for the ingestion and interrogation of Control Lists,
sets of regionally relevant data that connect accessions to annotations
to inform on regulation, or export control.

Control lists databases used for commec are downloaded from the commec
database repository, but can also be trivially added with custom controls.

"""

from .cli import add_args, format_control_lists
from .containers import (
    Accession,
    AccessionFormat,
    ControlList,
    ControlListContext,
    ControlListOutput,
    ListMode,
    Region,
    derive_accession_format,
)
from .control_list import (
    DESCRIPTION,
    get_cluster_hash,
    get_control_lists,
    get_regulation,
    import_data,
    is_regulated,
    run,
    should_ignore,
)
from .region import get_regions_set

__all__ = [
    run,
    DESCRIPTION,
    get_control_lists,
    get_regulation,
    is_regulated,
    import_data,
    get_cluster_hash,
    should_ignore,
    ListMode,
    ControlList,
    ControlListOutput,
    ControlListContext,
    Region,
    Accession,
    AccessionFormat,
    derive_accession_format,
    get_regions_set,
    add_args,
    format_control_lists,
]
