"""Tests for commec.control_list.containers — data structures and enums."""

import pytest

from commec.control_list.containers import (
    Accession,
    AccessionFormat,
    CategoryType,
    ControlList,
    ControlListContext,
    ControlListOutput,
    ListMode,
    ListUseAcronym,
    Region,
    derive_accession_format,
)


# ---------------------------------------------------------------------------
# Accession
# ---------------------------------------------------------------------------

class TestAccession:

    def test_from_valid_taxid(self):
        acc = Accession("12345")
        assert acc.code == "12345"
        assert acc.type == AccessionFormat.TAXID

    def test_from_none(self):
        acc = Accession(None)
        assert acc.code is None
        assert acc.type == AccessionFormat.UNKNOWN

    def test_from_empty_string(self):
        acc = Accession("")
        assert acc.code is None
        assert acc.type == AccessionFormat.UNKNOWN

    def test_hash_equality(self):
        """Two Accessions with the same code must hash identically (DataFrame index requirement)."""
        a = Accession("999")
        b = Accession("999")
        assert hash(a) == hash(b)
        assert a == b

    def test_hash_inequality(self):
        a = Accession("100")
        b = Accession("200")
        assert hash(a) != hash(b)

    def test_str(self):
        assert str(Accession("42")) == "42"

    def test_repr_contains_format_and_code(self):
        acc = Accession("42")
        r = repr(acc)
        assert "42" in r
        assert "TaxonomyID" in r

    def test_get_returns_tuple(self):
        code, fmt = Accession("555").get()
        assert code == "555"
        assert fmt == AccessionFormat.TAXID

    def test_get_format(self):
        assert Accession("1").get_format() == AccessionFormat.TAXID
        assert Accession(None).get_format() == AccessionFormat.UNKNOWN

    def test_usable_as_dict_key(self):
        d = {Accession("1"): "val"}
        assert d[Accession("1")] == "val"

    def test_usable_in_set(self):
        s = {Accession("1"), Accession("1"), Accession("2")}
        assert len(s) == 2


# ---------------------------------------------------------------------------
# derive_accession_format
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("input_val,expected", [
    ("12345", AccessionFormat.TAXID),
    ("0", AccessionFormat.TAXID),
    ("999999999", AccessionFormat.TAXID),
    ("abc", AccessionFormat.UNKNOWN),
    ("12.34", AccessionFormat.UNKNOWN),
    ("", AccessionFormat.UNKNOWN),
    ("12 34", AccessionFormat.UNKNOWN),
    ("ABC123", AccessionFormat.UNKNOWN),
])
def test_derive_accession_format(input_val, expected):
    assert derive_accession_format(input_val) == expected


# ---------------------------------------------------------------------------
# Region
# ---------------------------------------------------------------------------

class TestRegion:

    def test_str(self):
        assert str(Region("New Zealand", "NZ")) == "New Zealand"

    def test_repr(self):
        assert repr(Region("New Zealand", "NZ")) == "New Zealand [NZ]"

    def test_hash_same_acronym(self):
        a = Region("European Union", "EU")
        b = Region("European Union", "EU")
        assert hash(a) == hash(b)

    def test_hash_different_acronym(self):
        assert hash(Region("New Zealand", "NZ")) != hash(Region("Australia", "AU"))

    def test_equality(self):
        assert Region("New Zealand", "NZ") == Region("New Zealand", "NZ")

    def test_inequality_different_name(self):
        assert Region("New Zealand", "NZ") != Region("Aotearoa", "NZ")

    def test_defaults(self):
        r = Region()
        assert r.name == ""
        assert r.acronym == ""


# ---------------------------------------------------------------------------
# ListMode
# ---------------------------------------------------------------------------

def test_list_mode_values():
    assert ListMode.COMPLIANCE == "Compliance"
    assert ListMode.CONDITIONAL_NUM == "Conditional Compliance"
    assert ListMode.COMPLIANCE_WARN == "Comply with Warning"
    assert ListMode.IGNORE == "Ignored"


def test_list_mode_is_str_enum():
    for member in ListMode:
        assert isinstance(member, str)


# ---------------------------------------------------------------------------
# ListUseAcronym
# ---------------------------------------------------------------------------

def test_list_use_acronym_values():
    assert ListUseAcronym.EXPORTCONTROLS == "EXPORT"
    assert ListUseAcronym.LICENCED == "LICENCE"
    assert ListUseAcronym.OTHERPATHOGEN == "PATHGN"


# ---------------------------------------------------------------------------
# ControlList
# ---------------------------------------------------------------------------

class TestControlList:

    def _make(self, **overrides):
        defaults = dict(
            name="The List", name_translated="Listus", acronym="SCM", url="http://example.com",
            region=Region("European Union", "EU"),
            status=ListMode.COMPLIANCE, use=ListUseAcronym.EXPORTCONTROLS,
        )
        defaults.update(overrides)
        return ControlList(**defaults)

    def test_str_format(self):
        assert str(self._make()) == "EU_SCM_EXPORT"

    def test_description_contains_key_fields(self):
        desc = self._make().description()
        assert "The List" in desc
        assert "European Union" in desc
        assert "http://example.com" in desc
        assert "EU_SCM_EXPORT" in desc

    def test_equality_same_core_fields(self):
        """__eq__ checks name, url, region, acronym — not status or use."""
        a = self._make(status=ListMode.COMPLIANCE, use=ListUseAcronym.EXPORTCONTROLS)
        b = self._make(status=ListMode.IGNORE, use=ListUseAcronym.LICENCED)
        assert a == b

    def test_inequality_different_name(self):
        assert self._make(name="A") != self._make(name="B")

    def test_inequality_different_url(self):
        assert self._make(url="http://a.com") != self._make(url="http://b.com")

    def test_inequality_different_acronym(self):
        assert self._make(acronym="X") != self._make(acronym="Y")

    def test_inequality_different_region(self):
        assert (
            self._make(region=Region("NZ", "NZ"))
            != self._make(region=Region("AU", "AU"))
        )



# ---------------------------------------------------------------------------
# ControlListOutput
# ---------------------------------------------------------------------------

def test_control_list_output_defaults():
    out = ControlListOutput()
    assert out.name == ""
    assert out.category == ""
    assert out.list == ""
    assert out.species == ""
    assert out.genus == ""


def test_control_list_output_fields():
    out = ControlListOutput("Influenza A", "Viruses", "L1", "flu_species", "flu_genus")
    assert out.name == "Influenza A"
    assert out.category == "Viruses"
    assert out.list == "L1"
    assert out.species == "flu_species"
    assert out.genus == "flu_genus"


# ---------------------------------------------------------------------------
# ControlListContext
# ---------------------------------------------------------------------------

def test_control_list_context_defaults():
    ctx = ControlListContext()
    assert ctx.derived_from is None
    assert ctx.is_child is False


def test_control_list_context_fields():
    ctx = ControlListContext(derived_from="parent_organism", is_child=True)
    assert ctx.derived_from == "parent_organism"
    assert ctx.is_child is True
