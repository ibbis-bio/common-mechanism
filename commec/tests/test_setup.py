from commec.setup import CliSetup

def test_setup_tutorial_url_is_valid():
    """
    Ensure that the default tutorial database is a valid URL.
    """
    dont_automate : bool = False
    psuedo_setup = CliSetup(dont_automate)
    is_tutorial_url_valid : bool = psuedo_setup.check_url_exists(
        psuedo_setup.default_tutorial_download_url
        )
    assert is_tutorial_url_valid
