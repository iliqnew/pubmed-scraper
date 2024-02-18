import json

from Bio import Entrez
from typing import Literal, Union, Any


QUERY = "Ulcerative Colitis"
EMAIL = "iliqanew@email.com"

with open("person_profiles.json", "rb") as f:
    STORED_PROFILES = json.loads(f.read())

AUTHOR_T = {
    # "_id": {"$oid": "5c6ad92e1fba1942ccb101c2"},
    "firstName": None,
    # "middleName": None,
    "lastName": None,  # "Tjalling W de Vries, dr.",
    "affiliation": None,
    # "degrees": None,
    # "phone": "0031582863390",
    # "phoneExt": None,
    # "email": "tjalling.de.vries@znb.nl",
    # "roleSpecified": True,
    # "role": 0,
    # "roles": [1],
    # "personType": 0,
    # "name": "Tjalling W de Vries, dr.",
    # "formattedName": "Tjalling De Vries",
    # "isMarkedAsDuplicate": False,
    # "lastDuplicationCheck": {"$date": "2022-04-26T01:44:11.340Z"},
    # "personTypeFixDate": {"$date": {"$numberLong": "-62135596800000"}},
    # "createdOn": {"$date": {"$numberLong": "-62135596800000"}},
    # "debarmentData": None,
    # "formattedAddress": "Leeuwarden, Friesland, Netherlands",
    # "website": None,
    # "placeId": None,
    # "location": None,
    # "contactsDownloadedFrom": None,
    # "googleSERPProcessed": False,
    # "notRealPerson": False,
    # "npi": None,
    # "medicalSchool": None,
    # "specialties": None,
    # "affiliateHospitals": None,
    # "medicareProcessed": True,
    # "medicareProcessedFix": True,
    # "phonesInfo": [
    #     {"trialId": "NCT03654508", "value": "0031582863390"},
    #     {"trialId": "NCT02197780", "value": "0031 582 866 666"},
    # ],
    # "emailsInfo": [
    #     {"trialId": "NCT03654508", "value": "tjalling.de.vries@znb.nl"},
    #     {"trialId": "NCT02197780", "value": "Tjalling.de.Vries@ZNB.nl"},
    # ],
    # "emailsInfoUpdated": True,
    # "fDAInspections": None,
    # "normalizedEmails": ["TJALLING.DE.VRIES@ZNB.NL"],
    # "normalizedPhones": ["0031582863390", "0031582866666"],
    # "trials": ["NCT03654508", "NCT02197780"],
    "isStored": None,
}


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


def parse_studies(raw_studies):
    for raw_study in raw_studies:
        raw_authors = raw_study["MedlineCitation"]["Article"].get("AuthorList") or []
        for raw_author in raw_authors:
            author = AUTHOR_T.copy()

            cleaned_names = clean_names(
                raw_author.get("ForeName"), raw_author.get("LastName")
            )
            first_name = cleaned_names["firstName"]
            last_name = cleaned_names["lastName"]
            if first_name is last_name is None:
                continue

            author.update({"firstName": first_name, "lastName": last_name})

            is_stored = clean_is_stored(author)

            author.update(is_stored)

            yield author


def search_publications(query: str, retstart: int = 0, retmax: int = 9999):
    Entrez.email = EMAIL
    handle = Entrez.esearch(
        db="pubmed",
        sort="relevance",
        retmax=str(retmax),
        retstart=str(retstart),
        retmode="json",
        term=query,
    )
    return json.loads(handle.read())


def get_publications_by_ids(id_list):
    ids = ",".join(id_list)
    Entrez.email = EMAIL
    handle = Entrez.efetch(db="pubmed", retmode="xml", id=ids)
    raw_publications = Entrez.read(handle)
    return raw_publications


def paginate_publicaitons_ids(query: str, batch_size: int = 9999):
    retstart = 0
    retmax = batch_size
    studies_ids = search_publications(query, retstart=retstart, retmax=retmax)[
        "esearchresult"
    ]
    studies_id_list = studies_ids["idlist"]

    id_count = int(studies_ids["count"])
    while retmax < id_count:
        raw_studies = get_publications_by_ids(
            studies_id_list,
        )
        for author in parse_studies(raw_studies["PubmedArticle"]):
            yield author

        retstart = retmax + 1
        retmax += batch_size

        studies_ids = search_publications(query, retstart=retstart, retmax=retmax)[
            "esearchresult"
        ]
        studies_id_list = studies_ids["idlist"]


for author in paginate_publicaitons_ids(QUERY, batch_size=10):
    print(author)
