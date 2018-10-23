import pytest
from scripts import collapse_homopolymers

@pytest.mark.parametrize("data", [
    ("GGGGGGGCCCCCCCCCAAAAAAAAAATTTTTTTTTTTT", 1, "GCAT"),
    ("GGGGGGGCCCCCCCCCAAAAAAAAAATTTTTTTTTTTT", 2, "GGCCAATT"),
    ("GGGGGGGCCCCCCCCCAAAAAAAAAATTTTTTTTTTTT", 200, "GGGGGGGCCCCCCCCCAAAAAAAAAATTTTTTTTTTTT")])
def test_collapse_homopolymers_valid(data):
    assert collapse_homopolymers.collapse_homopolymers(data[0], max_length=data[1]) == data[2]

@pytest.mark.parametrize("data", [
    ("GGGGGGGCCCCCCCCCAAAAAAAAAATTTTTTTTTTTT", 0),
    ("GGGGGGGCCCCCCCCCAAAAAAAAAATTTTTTTTTTTT", -5)])
def test_collapse_homopolymers_invalid(data):
    with pytest.raises(Exception):
        collapse_homopolymers.collapse_homopolymers(data[0], max_length=data[1])
