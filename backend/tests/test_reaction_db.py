import csv
import gzip

import pytest

from app.reaction_db import ReactionTemplateStore

FIELDS = [
    "TEMPLATE_ID", "SCORE", "TEMPLATE", "REACTIONS", "RADII", "REACTIONS_COUNT",
    "ECS", "ECS_COUNT", "RADIUS_MIN", "RADIUS_MAX", "VALID", "DATASETS",
]


def _write_dataset(path):
    rows = [
        ["RR:test-001", "", "[O;H1:1]>>[C;H3]-[O;H0:1]", "USPTOB_1;USPTOB_2", "0", "2", "", "0", "0", "0", "True", "uspto"],
        ["RR:test-002", "", "[C;H3]-[O;H0:1]>>[O;H1:1]", "USPTOB_3", "0", "1", "", "0", "0", "0", "False", "uspto"],
        ["RR:test-003", "", "[C:1]>>[C:1]", "MNX_1", "1", "1", "1.1.1.1", "1", "1", "1", "True", "metanetx;rhea"],
    ]
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(FIELDS)
        writer.writerows(rows)


def test_reaction_database_index_and_search(tmp_path):
    source = tmp_path / "templates.csv.gz"
    database = tmp_path / "index.sqlite3"
    _write_dataset(source)
    store = ReactionTemplateStore(source, database)

    before = store.status()
    assert before["dataset_exists"] is True
    assert before["indexed"] is False

    indexed = store.build_index()
    assert indexed["indexed"] is True
    assert indexed["template_count"] == "3"
    assert indexed["valid_count"] == "2"
    assert indexed["invalid_count"] == "1"
    assert {x["dataset"] for x in indexed["source_distribution"]} == {"uspto", "metanetx", "rhea"}

    results = store.search(dataset="uspto", valid=True, limit=10)
    assert results["total"] == 1
    assert results["results"][0]["template_id"] == "RR:test-001"

    by_id = store.get("RR:test-003", validate_smarts=False)
    assert by_id["datasets"] == "metanetx;rhea"
    assert by_id["valid"] is True

    unchanged = store.build_index()
    assert "already matches" in unchanged["message"]


def test_reaction_database_rejects_missing_columns(tmp_path):
    source = tmp_path / "bad.csv.gz"
    with gzip.open(source, "wt", encoding="utf-8", newline="") as handle:
        handle.write("TEMPLATE_ID,TEMPLATE\nRR:1,[C:1]>>[C:1]\n")
    store = ReactionTemplateStore(source, tmp_path / "index.sqlite3")
    with pytest.raises(ValueError, match="missing required columns"):
        store.build_index()


def test_reaction_database_search_requires_index(tmp_path):
    source = tmp_path / "missing.csv.gz"
    store = ReactionTemplateStore(source, tmp_path / "index.sqlite3")
    with pytest.raises(RuntimeError, match="has not been indexed"):
        store.search(q="amide")


def test_reaction_database_api_routes_are_registered():
    from app.main import app
    paths = app.openapi()["paths"]
    assert "/api/reactions/status" in paths
    assert "/api/reactions/index" in paths
    assert "/api/reactions/search" in paths
    assert "/api/reactions/template/{template_id}" in paths
