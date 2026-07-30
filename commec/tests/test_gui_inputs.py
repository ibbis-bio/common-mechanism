"""Input handling for the Commec GUI."""

import json
from io import BytesIO

import openpyxl
import pytest
import xlwt
import yaml

from commec.gui import server


@pytest.fixture
def client(tmp_path, monkeypatch):
    runs_dir = tmp_path / "runs"
    runs_dir.mkdir()
    databases = tmp_path / "databases"
    control_lists = databases / "control_lists"
    control_lists.mkdir(parents=True)
    (control_lists / "region_definitions.json").write_text(
        json.dumps(
            [
                {
                    "name": "Australia Group",
                    "acronym": "AG",
                    "regions": ["AU", "CA", "US"],
                },
                {
                    "name": "United Kingdom",
                    "acronym": "UK",
                    "regions": ["GB"],
                },
            ]
        ),
        encoding="utf-8",
    )
    direct_list = control_lists / "direct"
    direct_list.mkdir()
    (direct_list / "list_info.csv").write_text(
        "region_code\nBR\n",
        encoding="utf-8",
    )

    monkeypatch.setitem(server.CFG, "password_hash", None)
    monkeypatch.setitem(server.CFG, "require_local_auth", False)
    monkeypatch.setitem(server.CFG, "runs_dir", runs_dir)
    monkeypatch.setitem(server.CFG, "default_databases", str(databases))
    monkeypatch.setitem(
        server.CFG,
        "presets",
        [{"id": "test", "label": "Test", "config": {}}],
    )
    monkeypatch.setattr(server, "_missing_databases", lambda _preset: [])
    monkeypatch.setattr(server.threading.Thread, "start", lambda _thread: None)
    server.JOBS.clear()
    server.app.config.update(TESTING=True)
    yield server.app.test_client()
    server.JOBS.clear()


def _submit(client, **overrides):
    data = {
        "preset": "test",
        "label": "input-test",
        "sequence_text": ">short\nACGT\n",
    }
    data.update(overrides)
    return client.post(
        "/screen",
        data=data,
        environ_overrides={"REMOTE_ADDR": "127.0.0.1"},
    )


def _xlsx_file():
    book = openpyxl.Workbook()
    samples = book.active
    samples.title = "Samples"
    samples.append(["sample_id", "nucleotide", "note"])
    samples.append(["alpha", "ACGTACGTACGT", "first"])
    samples.append(["beta", "TTTTCCCCAAAA", "second"])
    other = book.create_sheet("Other")
    other.append(["name", "sequence"])
    other.append(["gamma", "GGGGAAAACCCC"])
    data = BytesIO()
    book.save(data)
    data.seek(0)
    return data


def _xls_file():
    book = xlwt.Workbook()
    sheet = book.add_sheet("Legacy")
    sheet.write(0, 0, "name")
    sheet.write(0, 1, "sequence")
    sheet.write(1, 0, "legacy")
    sheet.write(1, 1, "ACACACACACAC")
    data = BytesIO()
    book.save(data)
    data.seek(0)
    return data


def _long_xlsx_file():
    book = openpyxl.Workbook()
    sheet = book.active
    sheet.title = "Long samples"
    sheet.append(["sample_id", "nucleotide"])
    for i in range(100):
        sheet.append([f"sample_{i + 1:03d}", "ACGT" * 25])
    data = BytesIO()
    book.save(data)
    data.seek(0)
    return data


def test_short_sequences_are_passed_to_commec(tmp_path, client):
    """The GUI leaves sequence-length handling to the screening pipeline."""
    response = _submit(
        client,
        # Older clients may still submit the removed field. It must not
        # restore the GUI's former length filter.
        skip_short="1",
    )

    assert response.status_code == 200
    run_dir = next(path for path in (tmp_path / "runs").iterdir() if path.is_dir())
    assert (run_dir / "input.fasta").read_text(encoding="utf-8") == ">short\nACGT\n"
    config = yaml.safe_load((run_dir / "config.used.yaml").read_text(encoding="utf-8"))
    assert config["databases"]["control_lists"]["regions"] == "all"


def test_selected_regions_are_recorded_in_run_config(tmp_path, client):
    response = _submit(client, regions="US,CA")

    assert response.status_code == 200
    run_dir = next(path for path in (tmp_path / "runs").iterdir() if path.is_dir())
    config = yaml.safe_load((run_dir / "config.used.yaml").read_text(encoding="utf-8"))
    assert config["databases"]["control_lists"]["regions"] == "US,CA"


def test_hidden_single_country_alias_is_accepted(tmp_path, client):
    response = _submit(client, regions="UK")

    assert response.status_code == 200
    run_dir = next(path for path in (tmp_path / "runs").iterdir() if path.is_dir())
    config = yaml.safe_load((run_dir / "config.used.yaml").read_text(encoding="utf-8"))
    assert config["databases"]["control_lists"]["regions"] == "UK"


@pytest.mark.parametrize("region", ["NOT_A_REGION", "AQ"])
def test_unsupported_region_is_rejected(tmp_path, client, region):
    response = _submit(client, regions=region)

    assert response.status_code == 400
    assert response.get_json() == {
        "error": f"Unknown regulatory jurisdiction: {region}."
    }
    assert list((tmp_path / "runs").iterdir()) == []


def test_config_exposes_supported_region_choices(client):
    response = client.get(
        "/config",
        environ_overrides={"REMOTE_ADDR": "127.0.0.1"},
    )

    assert response.status_code == 200
    regions = response.get_json()["regions"]
    assert {
        "code": "US",
        "name": "United States",
        "group": False,
        "memberships": ["Australia Group"],
    } in regions
    assert {
        "code": "BR",
        "name": "Brazil",
        "group": False,
        "memberships": [],
    } in regions
    assert not any(region["code"] == "AQ" for region in regions)
    assert [region for region in regions if region["code"] == "AG"] == [
        {"code": "AG", "name": "Australia Group", "group": True}
    ]
    assert not any(region["code"] == "UK" for region in regions)
    assert {
        "code": "GB",
        "name": "United Kingdom",
        "group": False,
        "memberships": [],
    } in regions


def test_gui_serves_bundled_report_font(client):
    response = client.get("/fonts/CrimsonPro.woff2")

    assert response.status_code == 200
    assert response.mimetype == "font/woff2"
    assert response.data[:4] == b"wOF2"


def test_xlsx_preview_lists_sheets_and_rows(client):
    response = client.post(
        "/spreadsheet-preview",
        data={"spreadsheet_file": (_xlsx_file(), "samples.xlsx")},
        content_type="multipart/form-data",
    )

    assert response.status_code == 200
    body = response.get_json()
    assert body["sheets"] == ["Samples", "Other"]
    assert body["sheet"] == "Samples"
    assert body["rows"][1][:2] == ["alpha", "ACGTACGTACGT"]


@pytest.mark.parametrize("filename,contents,sheet", [
    ("samples.xlsx", _xlsx_file, "Samples"),
    ("legacy.xls", _xls_file, "Legacy"),
    ("samples.csv", lambda: BytesIO(b"name,sequence\ncsv,CCCCAAAATTTT\n"), "Data"),
    ("samples.tsv", lambda: BytesIO(b"name\tsequence\ntsv\tAAAACCCCGGGG\n"), "Data"),
])
def test_tabular_upload_maps_selected_columns(tmp_path, client, filename, contents, sheet):
    response = _submit(
        client,
        sequence_text="",
        spreadsheet_file=(contents(), filename),
        spreadsheet_sheet=sheet,
        spreadsheet_nucleotide_column="1",
        spreadsheet_name_column="0",
        spreadsheet_skip_first_row="1",
    )

    assert response.status_code == 200
    run_dir = next(path for path in (tmp_path / "runs").iterdir() if path.is_dir())
    fasta = (run_dir / "input.fasta").read_text(encoding="utf-8")
    assert ">" in fasta
    assert "sequence" not in fasta


def test_tabular_upload_can_generate_names(client, tmp_path):
    response = _submit(
        client,
        sequence_text="",
        spreadsheet_file=(_xlsx_file(), "samples.xlsx"),
        spreadsheet_sheet="Other",
        spreadsheet_nucleotide_column="1",
        spreadsheet_name_column="",
        spreadsheet_skip_first_row="1",
    )

    assert response.status_code == 200
    run_dir = next(path for path in (tmp_path / "runs").iterdir() if path.is_dir())
    assert ">Other_row_2" in (run_dir / "input.fasta").read_text(encoding="utf-8")


def test_long_spreadsheet_preview_is_capped_but_all_rows_are_imported(client, tmp_path):
    preview = client.post(
        "/spreadsheet-preview",
        data={"spreadsheet_file": (_long_xlsx_file(), "long-samples.xlsx")},
        content_type="multipart/form-data",
    )

    assert preview.status_code == 200
    preview_data = preview.get_json()
    assert len(preview_data["rows"]) == server._TABULAR_PREVIEW_ROWS == 3
    assert preview_data["has_more_rows"] is True

    response = _submit(
        client,
        sequence_text="",
        spreadsheet_file=(_long_xlsx_file(), "long-samples.xlsx"),
        spreadsheet_sheet="Long samples",
        spreadsheet_nucleotide_column="1",
        spreadsheet_name_column="0",
        spreadsheet_skip_first_row="1",
    )

    assert response.status_code == 200
    run_dir = next(path for path in (tmp_path / "runs").iterdir() if path.is_dir())
    fasta = (run_dir / "input.fasta").read_text(encoding="utf-8")
    assert fasta.count(">sample_") == 100
