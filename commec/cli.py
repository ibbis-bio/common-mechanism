#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Command-line entrypoint for the package. Exposes the following subcommands:

    setup   Download or update the reference databases required for screening
    screen  Run Common Mechanism screening on an input FASTA.
    flag    Parse .screen or .json files in a directory and create CSVs of flags raised
    list    Display information on available annotated control lists

Command-line usage examples:
    - commec screen -y /path/to/config.yaml input.fasta
    - commec flag /path/to/directory/with/output.screen
    - commec -h, --help
    - commec -v, --version
"""

from commec.control_list import (
    DESCRIPTION as list_DESCRIPTION,
    add_args as list_add_args,
    run as list_run,
)
from commec.flag import (
    DESCRIPTION as flag_DESCRIPTION,
    add_args as flag_add_args,
    run as flag_run,
)
from commec.screen import (
    DESCRIPTION as screen_DESCRIPTION,
    add_args as screen_add_args,
    run as screen_run,
    ScreenArgumentParser,
)
from commec.setup import (
    DESCRIPTION as setup_DESCRIPTION,
    add_args as setup_add_args,
    run as setup_run,
)

from commec import __version__ as COMMEC_VERSION


def main():
    """
    Parse the command line arguments and call the relevant sub-command.
    """
    parser = ScreenArgumentParser(
        prog="commec", description="Command-line entrypoint for the Common Mechanism"
    )
    # Sub argument for version information
    parser.add_argument(
        "-v",
        "--version",
        dest="version",
        action="store_true",
        default=False,
        help="show version information and exit",
    )

    # Setup sub parsers:
    subparsers = parser.add_subparsers(dest="command")

    # Sub-command for "screen"
    screen_parser = subparsers.add_parser("screen", description=screen_DESCRIPTION)
    screen_add_args(screen_parser)

    # Sub-command for "flag"
    flag_parser = subparsers.add_parser("flag", description=flag_DESCRIPTION)
    flag_add_args(flag_parser)

    # Sub-command for "setup"
    setup_parser = subparsers.add_parser("setup", description=setup_DESCRIPTION)
    setup_add_args(setup_parser)

    # Sub-command for "list"
    list_parser = subparsers.add_parser("list", description=list_DESCRIPTION)
    list_add_args(list_parser)

    args = parser.parse_args()

    if args.command == "screen":
        screen_run(args)
    elif args.command == "flag":
        flag_run(args)
    elif args.command == "setup":
        setup_run(args)
    elif args.command == "list":
        list_run(args)
    elif args.version:
        print(
            "Commec  : The Common Mechanism\n"
            f"Version : {COMMEC_VERSION}\n"
            "Copyright IBBIS (c) 2021-2025\n"
            "International Biosecurity and Biosafety Initiative for Science"
        )
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
