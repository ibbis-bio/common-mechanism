# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science

"""
Responsible for the ingestion and interrogation of Control Lists, 
sets of regionally relevant data that connect accessions to annotations
to inform on regulation, or export control.

Control lists databases used for commec are downloaded from the commec
database repository, but can also be trivially added with custom controls.

"""
from .control_list import (
    run,
    DESCRIPTION,
    get_control_lists,
    get_regulation,
    is_regulated,
    import_data
)

from .cli import (
    add_args,
    format_control_lists
)

from .containers import (
    ListMode,
    ControlList,
    ControlListInfo,
    Region,
    Accession,
    AccessionFormat,
    derive_accession_format
)

from .region import (
    get_regions_set
)
