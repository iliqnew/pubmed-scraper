import json


from Bio import Entrez
from typing import Union, Literal, Any


EMAIL = "iliqanew@email.com"


with open("person_profiles.json", "rb") as f:
    STORED_PROFILES = json.loads(f.read())


def clean_names(
    fore_name: Union[str, None], last_name: Union[str, None]
) -> dict[Literal["firstName", "lastName"], Union[str, None]]:
    # apply data cleaning logic here

    return {"firstName": fore_name, "lastName": last_name}


def clean_is_stored(author: dict[str, Any]) -> dict[Literal["isStored"], bool]:
    fn = author["firstName"]
    ln = author["lastName"]

    _first_name = fn + " " if fn else ""
    _last_name = ln if ln else ""

    is_stored = _first_name + _last_name in [
        profile["lastName"].split(",")[0]
        for profile in STORED_PROFILES
        if not profile["lastName"] is None
    ]

    return {"isStored": is_stored}


def search_query(query: str, retstart: int = 0, retmax: int = 9999):
    Entrez.email = EMAIL
    handle = Entrez.esearch(
        db="pubmed",
        sort="relevance",
        retmax=str(retmax),
        retstart=str(retstart),
        retmode="xml",
        term=query,
    )
    ids = Entrez.read(handle)
    return ids


def get_publications_by_ids(id_list):
    ids = ",".join(id_list)
    Entrez.email = EMAIL
    handle = Entrez.efetch(db="pubmed", retmode="xml", id=ids)
    raw_publications = Entrez.read(handle)
    return raw_publications
